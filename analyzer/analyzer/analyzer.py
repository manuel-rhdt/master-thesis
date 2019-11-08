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


def generate_signal(name, max_time=1000, step_size=0.01, x0=500, mean=500, correlation_time=10, diffusion=10):
    timestamps = np.arange(0, max_time, step_size)
    x = np.clip(ornstein_uhlenbeck.generate(
        timestamps, x0, correlation_time, diffusion, mean=mean), 0.0, None)
    return {'timestamps': timestamps, 'components': {name: x}}


@jit(nopython=True, fastmath=True)
def calculate_sum_of_reaction_propensities(components, reaction_k, reaction_reactants):
    """ Accelerated reaction propensity calculation

    Arguments:
        components {np.ndarray} -- an array of shape (num_components, num_trajectories, num_events_per_trajectory)
        reaction_k {np.ndarray} -- [description]
        reaction_reactants {np.ndarray} -- a two-dimensional array (num_reactions, num_reactants)

    Returns:
        [type] -- [description]
    """
    result = np.empty_like(components[0])

    for n_reaction in range(len(reaction_k)):
        tmp = np.full_like(components[0], reaction_k[n_reaction])

        for j_reactant in reaction_reactants[n_reaction]:
            if j_reactant >= 0:
                tmp *= components[j_reactant]

        if n_reaction == 0:
            result = tmp
        else:
            result += tmp

    return result


@jit(nopython=True, fastmath=True)
def calculate_selected_reaction_propensities(components, reaction_k, reaction_reactants, reaction_events):
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

    propensities = np.empty(reaction_k.shape + components[0].shape)
    for n_reaction in range(len(reaction_k)):
        propensities[n_reaction] = reaction_k[n_reaction]

        for j_reactant in reaction_reactants[n_reaction]:
            if j_reactant >= 0:
                propensities[n_reaction] *= components[j_reactant]

    result = np.empty_like(components[0])
    for i in range(length):
        event = reaction_events[i]
        result[i] = propensities[event, i]

    return result


@cuda.jit
def gpu_selected_reaction_propensities(components, reaction_k, reaction_reactants, reaction_events, propensities):
    pos = cuda.grid(1)
    if pos < reaction_events.size:
        event = reaction_events[pos]
        propensities[pos] = reaction_k[event]
        for j_reactant in reaction_reactants[event]:
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
def time_average(trajectory, old_timestamps, new_timestamps, out=None):
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
              |========|===========================| <== average
              |        |                           |
              |--------+                           |
              |                                    |
              +------------------------------------+---> time
        old_timestamps[i]                  old_timestamps[i+1]
    """
    if out is None:
        out = np.empty(len(new_timestamps) - 1, trajectory.dtype)
    old_idx = 0
    old_ts = old_timestamps[0]
    trajectory_value = trajectory[0]
    for new_idx in range(len(new_timestamps) - 1):
        low = new_timestamps[new_idx]
        high = new_timestamps[new_idx + 1]
        delta_t = high - low
        acc = 0.0
        while low < high:
            while old_ts <= low:
                old_idx += 1
                if old_idx == len(old_timestamps):
                    old_ts = np.inf
                    trajectory_value = trajectory[-1]
                else:
                    old_ts = old_timestamps[old_idx]
                    trajectory_value = trajectory[old_idx - 1]

            acc += trajectory_value * (min(old_ts, high) - low)
            low = old_ts

        out[new_idx] = acc / delta_t
    return out


@jit(nopython=True, fastmath=True)
def log_likelihood_inner(signal_components, signal_timestamps, response_components, response_timestamps, reaction_k, reaction_reactants, reaction_events, out=None):
    num_signal_comps, _ = signal_components.shape
    num_response_comps, length = response_components.shape

    components = TypedList()
    for i in range(num_signal_comps):
        components.append(np.empty((length - 1)))
        time_average(signal_components[i], signal_timestamps,
                     response_timestamps, out=components[i])
    for i in range(num_response_comps):
        components.append(response_components[i, 1:])

    averaged_rates = calculate_sum_of_reaction_propensities(
        components, reaction_k, reaction_reactants)

    for i in range(num_signal_comps):
        # we don't evaluate the trajectory at the first timestamp since it only specifies the
        # initial value (no reaction occurecd)
        evaluate_trajectory_at(signal_components[i], signal_timestamps,
                               response_timestamps[1:], out=components[i])

    instantaneous_rates = calculate_selected_reaction_propensities(
        components, reaction_k, reaction_reactants, reaction_events)

    # return the logarithm of `np.cumprod(instantaneous_rates * np.exp(-averaged_rates * time_delta))`

    # perform the operations in-place
    dt = np.zeros(length-1)
    for i in range(length-1):
        dt[i] = response_timestamps[i+1] - response_timestamps[i]
    likelihoods = np.log(instantaneous_rates)
    likelihoods -= averaged_rates * dt

    result = out if out is not None else np.empty(length - 1)
    accumulator = 0.0
    for i in range(length-1):
        accumulator += likelihoods[i]
        result[i] = accumulator

    return result


@jit(nopython=True, fastmath=True, parallel=True)
def log_likelihood(signal_components, signal_timestamps, response_components, response_timestamps, reaction_k, reaction_reactants, reaction_events, out=None):
    num_r, _, length = response_components.shape
    num_s, _, _ = signal_components.shape

    assert num_r == num_s

    result = out if out is not None else np.empty((num_r, length - 1))
    for r in prange(num_r):
        rc = response_components[r]
        rt = response_timestamps[r]
        sc = signal_components[r]
        st = signal_timestamps[r]

        log_likelihood_inner(sc, st, rc, rt, reaction_k,
                             reaction_reactants, reaction_events[r], out=result[r])

    return result


@jit(nopython=True, fastmath=True, parallel=True)
def log_averaged_likelihood(signal_components, signal_timestamps, response_components, response_timestamps, reaction_k, reaction_reactants, reaction_events, out=None):
    num_r, _, length = response_components.shape
    num_s, _, _ = signal_components.shape

    result = out if out is not None else np.empty((num_r, length - 1))
    for r in prange(num_r):
        rc = response_components[r]
        rt = response_timestamps[r]
        tmp = np.empty_like(result[r])
        for s in range(num_s):
            sc = signal_components[s]
            st = signal_timestamps[s]

            log_p = log_likelihood_inner(
                sc, st, rc, rt, reaction_k, reaction_reactants, reaction_events[r])

            if s > 0:
                # We correctly average the likelihood over multiple signals if requested.
                #
                # Note that since we computed the log likelihoods we have to use the
                # logaddexp operation to correctly perform the averaging
                #
                # The next line of code performs the following computation:
                # result <- log(exp(result) + exp(log_p) / num_s)
                tmp = np.logaddexp(tmp, log_p - np.log(num_s))
            else:
                tmp = log_p - np.log(num_s)

        result[r] = tmp

    return result


def compile_reactions(response, signal):
    name_idx_map = {}

    signal_names = list(signal['components'].keys())
    response_names = list(response['components'].keys())

    for idx, name in enumerate(signal_names):
        name_idx_map[name] = idx

    for idx, name in enumerate(response_names):
        name_idx_map[name] = idx + len(signal_names)

    reactions = response['reactions']

    max_num_reactants = max(len(r['reactants']) for r in reactions)
    reaction_reactants = np.full(
        (len(reactions), max_num_reactants), -1, dtype=np.int32)
    for i, reaction in enumerate(reactions):
        for j, reactant in enumerate(reaction['reactants']):
            reaction_reactants[i, j] = name_idx_map[reactant]

    signal_components = TypedList()
    for name in signal_names:
        signal_components.append(signal['components'][name])

    response_components = TypedList()
    for name in response_names:
        response_components.append(response['components'][name])

    reaction_k = np.array([r['k'] for r in reactions])

    return (signal_components, response_components, reaction_k, reaction_reactants)


def log_likelihoods_given_signal(response, signal):
    """ Calculates the individual transition probabilities for a response trajectory given a signal.

    Returns an array of transition probabilities which when added together result in the log likelihood
    $log p(X|S)$ where X is the response trajectory and S is the signal trajectory.
    """

    signal_comps, response_comps, reaction_k, reaction_reactants = compile_reactions(
        response, signal)
    signal_ts = signal['timestamps']
    response_ts = response['timestamps']
    events = response['reaction_events']

    return log_averaged_likelihood(signal_comps, signal_ts, response_comps, response_ts, reaction_k, reaction_reactants, events)


def generate_input(name, num_components, num_reactions, num_blocks, num_steps):
    return '{}\n{}\n{}\n{} {} {}\n{}\n'.format(name, num_components, num_reactions, num_blocks, 0, num_steps, 100)


def empty_trajectory(components, observations, length, events=False):
    """Create an empty trajectory with preallocated arrays

    Arguments:
        components {[type]} -- [description]
        observations {[type]} -- [description]
        length {[type]} -- [description]

    Keyword Arguments:
        events {bool} -- [description] (default: {False})

    Returns:
        [type] -- [description]
    """
    obs = observations
    l = length
    trajectory = {
        'timestamps': np.zeros((obs, l)),
        'components': {}
    }
    for c in components:
        trajectory['components'][c] = np.zeros((obs, l))

    if events:
        trajectory['reaction_events'] = np.zeros((obs, l), dtype=np.int32)

    return trajectory
