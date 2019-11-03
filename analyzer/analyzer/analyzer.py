import json
import multiprocessing
import os
import pathlib
import subprocess
import sys

import msgpack
import numpy as np
import scipy
import scipy.stats
from scipy import interpolate
from numba import jit, guvectorize, float64, int32, prange
from numba.typed import List as TypedList

import settings

from . import ornstein_uhlenbeck


def generate_signal(name, max_time=1000, step_size=0.01, x0=500, mean=500, correlation_time=10, diffusion=10):
    timestamps = np.arange(0, max_time, step_size)
    x = np.clip(ornstein_uhlenbeck.generate(
        timestamps, x0, correlation_time, diffusion, mean=mean), 0.0, None)
    return {'timestamps': timestamps, 'components': {name: x}}


def save_signal(name, path, duration=1000, step_size=0.01, x0=500, mean=500, correlation_time=10, diffusion=10):
    with open(str(path), 'wb') as file:
        signal = generate_signal(name, duration, step_size, x0,
                                 mean, correlation_time, diffusion)
        signal['timestamps'] = signal['timestamps'].tolist()
        signal['components'][name] = signal['components'][name].tolist()
        msgpack.pack(signal, file)


# inspired by scipy.misc.logsumexp
@jit(nopython=True, fastmath=True)
def logsumexp(x, scale=1):
    x = np.asarray(x)
    x_max = np.amax(x)
    tmp = scale * np.exp(x - x_max)

    s = np.sum(tmp)
    out = np.log(s)

    return out + x_max


@jit(nopython=True, fastmath=True)
def _calculate_sum_of_reaction_propensities(components, reaction_k, reaction_reactants, result=None):
    """ Accelerated reaction propensity calculation

    Arguments:
        components {np.ndarray} -- an array of shape (num_components, num_trajectories, num_events_per_trajectory)
        reaction_k {np.ndarray} -- [description]
        reaction_reactants {np.ndarray} -- a two-dimensional array (num_reactions, num_reactants)

    Returns:
        [type] -- [description]
    """
    traj_length = components[0].shape[-1]

    if result is None:
        result = np.zeros_like(components[0])

    for i in prange(traj_length):
        for n_reaction in range(len(reaction_k)):
            tmp = reaction_k[n_reaction]

            for j_reactant in reaction_reactants[n_reaction]:
                if j_reactant >= 0:
                    tmp *= components[j_reactant][i]

            result[i] += tmp

    return result


@jit(nopython=True, fastmath=True)
def _calculate_selected_reaction_propensities(components, reaction_k, reaction_reactants, reaction_events):
    """ Accelerated reaction propensity calculation

    Arguments:
        components {np.ndarray} -- an array of shape (num_components, num_trajectories, num_events_per_trajectory)
        reaction_k {np.ndarray} -- [description]
        reaction_reactants {np.ndarray} -- a two-dimensional array (num_reactions, num_reactants)

    Returns:
        [type] -- [description]
    """
    result = np.zeros_like(components[0])

    traj_length = components[0].shape[-1]

    assert reaction_events.shape == components[0].shape

    for i in range(traj_length):
        event = reaction_events[i]
        result[i] = reaction_k[event]
        for j_reactant in reaction_reactants[event]:
            if j_reactant >= 0:
                result[i] *= components[j_reactant][i]

    return result


@jit(nopython=True, fastmath=True)
def resample_trajectory(trajectory, old_timestamps, new_timestamps, out=None):
    """ Resample trajectory with `old_timestamp` at the times in `new_timestamps`.

    This function assumes that `new_timestamps` is ordered.

    Arguments:
        trajectory {[type]} -- [description]
        old_timestamps {[type]} -- [description]
        new_timestamps {[type]} -- [description]

    Returns:
        [type] -- an array of values
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
def resample_averaged_trajectory(trajectory, old_timestamps, new_timestamps, out=None):
    """ Resample and averagte trajectory with `old_timestamp` at the times in `new_timestamps`.

    This function assumes that `new_timestamps` is ordered.

    Arguments:
        trajectory {[type]} -- [description]
        old_timestamps {[type]} -- [description]
        new_timestamps {[type]} -- [description]

    Returns:
        [type] -- an array of values
    """
    if out is None:
        out = np.empty(new_timestamps.shape, trajectory.dtype)
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
    out[-1] = trajectory[min(old_idx, len(old_timestamps)-1)]
    return out


@jit(nopython=True)
def log_likelihood(signal_components, signal_timestamps, response_components, response_timestamps, reaction_k, reaction_reactants, reaction_events):
    length, = response_components[0].shape

    components = TypedList()
    for comp in signal_components:
        resampled = np.zeros(length)
        resample_averaged_trajectory(
            comp, signal_timestamps, response_timestamps, resampled)
        components.append(resampled)
    for comp in response_components:
        components.append(comp)

    averaged_rates = _calculate_sum_of_reaction_propensities(
        components, reaction_k, reaction_reactants)

    for i, comp in enumerate(signal_components):
        resampled = np.zeros(length)
        resample_trajectory(comp, signal_timestamps,
                            response_timestamps, resampled)
        components[i] = resampled

    instantaneous_rates = _calculate_selected_reaction_propensities(
        components, reaction_k, reaction_reactants, reaction_events)

    # return the logarithm of `np.cumprod(instantaneous_rates * np.exp(-averaged_rates * time_delta))`
    result = np.zeros(length)
    for i in range(1, length):
        result[i] = result[i - 1] + np.log(instantaneous_rates[i]) - averaged_rates[i] * (
            response_timestamps[i] - response_timestamps[i-1])
    return result


@jit(nopython=True, parallel=False, fastmath=True, debug=True)
def log_averaged_likelihood(signal_components, signal_timestamps, response_components, response_timestamps, reaction_k, reaction_reactants, reaction_events):
    num_r, length = response_components[0].shape
    num_s, _ = signal_components[0].shape

    result = np.empty((num_r, length))
    for j in prange(num_r * num_s):
        r = j // num_s
        s = j % num_s

        sc = TypedList()
        for comp in signal_components:
            sc.append(comp[s])
        rc = TypedList()
        for comp in response_components:
            rc.append(comp[r])

        log_p = log_likelihood(sc, signal_timestamps[s], rc, response_timestamps[r],
                               reaction_k, reaction_reactants, reaction_events[r])

        if s > 0:
            # We correctly average the likelihood over multiple signals if requested.
            #
            # Note that since we computed the log likelihoods we have to use the
            # logaddexp operation to correctly perform the averaging
            result[r] = np.logaddexp(result[r], log_p - np.log(num_s))
        else:
            result[r] = log_p - np.log(num_s)

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


def simulate_trajectory(input_name, output_name=None, trajectories=None, cwd=settings.configuration['gillespie_cwd'], seed=4252):
    cmd = [os.path.abspath(
        str(settings.configuration['executable'])), input_name, '-s', str(seed)]
    if output_name is not None:
        cmd.extend(['-o', os.path.abspath(output_name)])

    if trajectories is not None:
        assert len(trajectories) >= 1
        for t in trajectories:
            cmd.extend(['-t', os.path.abspath(t)])

    subprocess.run(cmd, stdout=subprocess.DEVNULL, cwd=os.path.abspath(cwd),
                   timeout=10.0, check=True)


def load_trajectories(path, glob):
    trajectories = []
    for file_path in pathlib.Path(path).glob(glob):
        trajectories.append(load_trajectory(file_path))
    return trajectories


def load_trajectory(path):
    with path.open('rb') as f:
        trajectory = msgpack.unpack(f, raw=False)
        trajectory['timestamps'] = np.array(trajectory['timestamps'])
        for component in trajectory['components']:
            trajectory['components'][component] = np.array(
                trajectory['components'][component])

        if 'random_variates' in trajectory:
            trajectory['random_variates'] = np.array(
                trajectory['random_variates'])

        if 'reaction_events' in trajectory:
            trajectory['reaction_events'] = np.array(
                trajectory['reaction_events'])
    return trajectory


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
