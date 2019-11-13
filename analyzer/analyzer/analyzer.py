import json
import multiprocessing
import os
import pathlib
import subprocess
import sys

import numpy as np
import scipy
import scipy.stats
from scipy import interpolate
from numba import jit, guvectorize, float64, int32, prange, cuda
from numba.typed import List as TypedList

import settings

from . import ornstein_uhlenbeck


# inspired by scipy.special.logsumexp
@jit(nopython=True, fastmath=True)
def logsumexp(x):
    x = np.asarray(x)
    _, length = x.shape
    xmax = np.empty(length, dtype=x.dtype)
    for i in range(length):
        xmax[i] = max(x[:, i])

    tmp = np.exp(x - xmax)

    s = np.sum(tmp, axis=0)
    out = np.log(s)

    return out + xmax


def generate_signal(name, max_time=1000, step_size=0.01, x0=500, mean=500, correlation_time=10, diffusion=10):
    timestamps = np.arange(0, max_time, step_size)
    x = np.clip(ornstein_uhlenbeck.generate(
        timestamps, x0, correlation_time, diffusion, mean=mean), 0.0, None)
    return {'timestamps': timestamps, 'components': {name: x}}


@jit(nopython=True, fastmath=True)
def calculate_sum_of_reaction_propensities(components, reactions):
    """ Accelerated reaction propensity calculation

    Arguments:
        components {np.ndarray} -- an array of shape (num_components, num_trajectories, num_events_per_trajectory)
        reaction_k {np.ndarray} -- [description]
        reaction_reactants {np.ndarray} -- a two-dimensional array (num_reactions, num_reactants)

    Returns:
        [type] -- [description]
    """
    result = np.empty_like(components[0])

    for n_reaction in range(reactions.size):
        tmp = np.full_like(components[0], reactions.k[n_reaction])

        for j_reactant in reactions.reactants[n_reaction]:
            if j_reactant >= 0:
                tmp *= components[j_reactant]

        if n_reaction == 0:
            result = tmp
        else:
            result += tmp

    return result


@jit(nopython=True, fastmath=True)
def calculate_selected_reaction_propensities(components, reaction_events, reactions):
    """ Accelerated reaction propensity calculation

    Arguments:
        components {numba.typed.List} -- a list of components
        reaction_k {np.ndarray} -- [description]
        reaction_reactants {np.ndarray} -- a two-dimensional array (num_reactions, num_reactants)

    Returns:
        [type] -- [description]
    """
    length, = components[0].shape
    assert length == reaction_events.shape[-1]

    propensities = np.empty(
        reactions.k.shape + components[0].shape, dtype=components[0].dtype)
    for n_reaction in range(reactions.size):
        propensities[n_reaction] = reactions.k[n_reaction]

        for j_reactant in reactions.reactants[n_reaction]:
            if j_reactant >= 0:
                propensities[n_reaction] *= components[j_reactant]

    result = np.empty_like(components[0])
    for i in range(length):
        event = reaction_events[i]
        result[i] = propensities[event, i]

    return result


@cuda.jit
def gpu_selected_reaction_propensities(components, reaction_events, reactions, propensities):
    pos = cuda.grid(1)
    if pos < reactions.size:
        event = reaction_events[pos]
        propensities[pos] = reactions.k[event]
        for j_reactant in reactions.reactants[event]:
            if j_reactant >= 0:
                propensities[pos] *= components[j_reactant][pos]


@jit(nopython=True, fastmath=True)
def evaluate_trajectory_at(trajectory, old_timestamps, new_timestamps, out=None):
    """ Evaluate trajectory with events at `old_timestamps` at the times in `new_timestamps`.

    Note: This function assumes that both `old_timestamps` and `new_timestamps` are ordered.
    """
    if out is None:
        out = np.empty(new_timestamps.shape, trajectory.dtype)
    new_idx = 0
    for old_idx, ts in enumerate(old_timestamps):
        while new_timestamps[new_idx] < ts and new_idx < len(out):
            out[new_idx] = trajectory[max(old_idx - 1, 0)]
            new_idx += 1
        if new_idx >= len(out):
            break

    while new_idx < len(out):
        out[new_idx] = trajectory[-1]
        new_idx += 1
    return out


@jit(nopython=True, fastmath=True)
def time_average(trajectory, old_timestamps, new_timestamps, dtype=np.double, out=None, evaluated=None):
    """ Average of `trajectory` with `old_timestamp` in the time-intervals specified by `new_timestamps`.

    Note: This function assumes that both `old_timestamps` and `new_timestamps` are ordered.

    # Discussion

    This function is used for the calculation of the mean propensity over a timestep in the
    response trajectory. Since the propensity for a reaction is variable, the survival probability
    depends on the time-integrated value of the total propensity. Since the signal trajectory is
    assumed to be piecewise constant we can calculate the average signal between every
    pair of response timestamps. If we use this average signal to compute the total propensity,
    we don't require the calculation of the time integral anymore.

    This function computes the average of a trajectory given by the values `trajectory`
    at the timestamps `old_timestamps` where the averaging happens for the intervals between pairs
    of timestamps in `new_timestamps`.

    Returns a list of averages of size `len(new_timestamps) - 1`.

              |                                    |
              |        +---------------------------|
              |========|===========================| <== average[i]
              |        |                           |
              |--------+                           |
              |                                    |
              +------------------------------------+---> time
        old_timestamps[i]                  old_timestamps[i+1]
    """
    if out is None:
        out = np.empty(len(new_timestamps) - 1, dtype=dtype)

    old_idx = 0

    def iter_trajectory(old_idx):
        old_idx += 1
        if old_idx > len(old_timestamps) - 1:
            return np.inf, trajectory[-1]
        else:
            return old_timestamps[old_idx], trajectory[old_idx - 1]

    next_trajectory_change, trajectory_value = iter_trajectory(old_idx)
    old_idx += 1
    for idx, (low, high) in enumerate(zip(new_timestamps[:-1], new_timestamps[1:])):
        delta_t = high - low
        acc = 0.0
        while low < high:
            while next_trajectory_change <= low:
                next_trajectory_change, trajectory_value = iter_trajectory(
                    old_idx)
                old_idx += 1

            acc += trajectory_value * (min(high, next_trajectory_change) - low)
            low = next_trajectory_change

        out[idx] = acc / delta_t
        if evaluated is not None:
            evaluated[idx] = trajectory_value
    return out


@jit(nopython=True, fastmath=True)
def log_likelihood_inner(signal_components, signal_timestamps, response_components, response_timestamps, reaction_events, reactions, dtype=np.dtype(np.double), out=None):
    num_signal_comps, _ = signal_components.shape
    num_response_comps, length = response_components.shape

    if out is not None:
        dtype = out.dtype

    # resampled signal components
    rsc = np.empty((num_signal_comps, 2, length - 1), dtype=dtype)
    for i in range(num_signal_comps):
        time_average(signal_components[i], signal_timestamps,
                     response_timestamps, dtype=dtype, out=rsc[i, 0], evaluated=rsc[i, 1])

    components = TypedList()
    for i in range(num_signal_comps):
        components.append(rsc[i, 0])

    for i in range(num_response_comps):
        components.append(response_components[i, 1:].astype(dtype))

    averaged_rates = calculate_sum_of_reaction_propensities(
        components, reactions)

    for i in range(num_signal_comps):
        components[i] = rsc[i, 1]

    instantaneous_rates = calculate_selected_reaction_propensities(
        components, reaction_events, reactions)

    # return the logarithm of `np.cumprod(instantaneous_rates * np.exp(-averaged_rates * dt))`

    # perform the operations in-place (reuse the output array for the dt's)
    dt = out if out is not None else np.empty(length - 1, dtype=dtype)
    for i in range(length-1):
        dt[i] = response_timestamps[i+1] - response_timestamps[i]
    likelihoods = np.log(instantaneous_rates)
    likelihoods -= averaged_rates * dt

    result = dt
    accumulator = 0.0
    for i in range(length-1):
        accumulator += likelihoods[i]
        result[i] = accumulator

    return result


@jit(nopython=True, fastmath=True, parallel=True, cache=True)
def log_likelihood(traj_lengths, signal_components, signal_timestamps, response_components, response_timestamps, reaction_events, reactions, out=None):
    num_r, _, _ = response_components.shape
    num_s, _, _ = signal_components.shape
    length, = traj_lengths.shape

    assert num_r == num_s

    result = out if out is not None else np.empty(
        (num_r, length), dtype=np.single)
    for r in prange(num_r):
        rc = response_components[r]
        rt = response_timestamps[r]
        sc = signal_components[r]
        st = signal_timestamps[r]

        log_p = log_likelihood_inner(
            sc, st, rc, rt, reaction_events[r], reactions)
        indices = np.digitize(traj_lengths, rt)
        result[r] = log_p[indices]

    return result


@jit(nopython=True, fastmath=True, parallel=True, cache=True)
def log_averaged_likelihood(traj_lengths, signal_components, signal_timestamps, response_components, response_timestamps, reaction_events, reactions, p_zero, out=None):
    """
    Calculates the log likelihoods of the responses for various signals and averages over the signals.
    """
    num_r, _, _ = response_components.shape
    num_s, _, _ = signal_components.shape
    length, = traj_lengths.shape

    result = out if out is not None else np.empty(
        (num_r, length), dtype=np.single)
    for r in prange(num_r):
        rc = response_components[r]
        rt = response_timestamps[r]
        tmp = np.empty((num_s, length), dtype=np.single)
        indices = np.digitize(traj_lengths, rt)
        for s in range(num_s):
            sc = signal_components[s]
            st = signal_timestamps[s]

            log_p = log_likelihood_inner(
                sc, st, rc, rt, reaction_events[r], reactions, dtype=np.single)
            tmp[s] = log_p[indices] + p_zero[s, r]

        result[r] = logsumexp(tmp) - np.log(num_s)

    return result
