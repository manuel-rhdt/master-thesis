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
from numba import jit, guvectorize, float64, int32
from numba.typed import List as TypedList

import settings

from . import ornstein_uhlenbeck


def generate_signal(name, max_time=1000, length=100000, x0=500, mean=500, correlation_time=10, diffusion=10):
    timestamps = np.linspace(0, max_time, length)
    x = np.clip(ornstein_uhlenbeck.generate(
        timestamps, x0, correlation_time, diffusion, mean=mean), 0.0, None)
    return {'timestamps': timestamps, 'components': {name: x}}


def save_signal(name, path, duration=1000, length=100000, x0=500, mean=500, correlation_time=10, diffusion=10):
    with open(str(path), 'wb') as file:
        signal = generate_signal(name, duration, length, x0,
                                 mean, correlation_time, diffusion)
        signal['timestamps'] = signal['timestamps'].tolist()
        signal['components'][name] = signal['components'][name].tolist()
        msgpack.pack(signal, file)


@guvectorize([(float64[:, :], float64[:], int32[:, :], float64[:, :])], '(c,n),(r),(r,l)->(r,n)', nopython=True)
def _calculate_reaction_propensities(components, reaction_k, reaction_reactants, result=None):
    """ Accelerated reaction propensity calculation

    Arguments:
        components {np.ndarray} -- an array of shape (num_components, num_trajectories, num_events_per_trajectory)
        reaction_k {np.ndarray} -- [description]
        reaction_reactants {np.ndarray} -- a two-dimensional array (num_reactions, num_reactants)

    Returns:
        [type] -- [description]
    """
    # iterate over all reactions
    for i_reaction in range(reaction_k.shape[0]):
        result[i_reaction] = reaction_k[i_reaction]
        # iterate over all reactants
        for j_reactant in reaction_reactants[i_reaction]:
            if j_reactant >= 0:
                result[i_reaction] *= components[j_reactant]


def calculate_reaction_propensities(reactions, components):
    """ Uses the provided information about reactions and the component counts at each time step to
    calculate an array of reaction propensities at each timestamp.

    The outermost dimension corresponds to the number of reactions. The second outermost dimension
    corresponds to the number of observations
    """
    reaction_k = np.array([reaction['k'] for reaction in reactions])

    # turn a dictionary into an array
    keys = list(components.keys())
    component_array = None
    br = np.broadcast(*components.values())
    for i, x in enumerate(keys):
        x = np.asanyarray(components[x])
        if component_array is None:
            component_array = np.zeros((len(keys),) + br.shape, dtype=x.dtype)
        component_array[i] = x

    # we have to move the new axis
    # (components, responses, signals, timestamps) -> (responses, signals, components, timestamps)
    component_array = np.moveaxis(component_array, 0, 2)

    max_num_reactants = max(len(r['reactants']) for r in reactions)
    reaction_reactants = np.full(
        (len(reactions), max_num_reactants), -1, dtype=np.int32)
    for i, reaction in enumerate(reactions):
        for j, reactant in enumerate(reaction['reactants']):
            reaction_reactants[i, j] = keys.index(reactant)

    # return value has shape
    # (responses, signals, reactions, propensities)
    return _calculate_reaction_propensities(
        component_array, reaction_k, reaction_reactants)


# @guvectorize([(float64[:], float64[:], float64[:], float64[:])], '(o),(o),(n)->(n)', nopython=True)
@jit(nopython=True)
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


# @guvectorize([(float64[:], float64[:], float64[:], float64[:])], '(o),(o),(n)->(n)', nopython=True)
@jit(nopython=True)
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

        out[new_idx + 1] = acc / delta_t

@jit(nopython=True)
def resample_trajectory_outer(trajectory, old_timestamps, new_timestamps, output=None, averaged=False):
    assert len(trajectory.shape) == 2
    t = trajectory
    ots = old_timestamps
    nts = new_timestamps

    if output is None:
        output = np.empty((nts.shape[0], ots.shape[0], nts.shape[1]))

    for i in range(ots.shape[0]):
        for j in range(nts.shape[0]):
            if averaged:
                resample_averaged_trajectory(t[i], ots[i], nts[j], output[j, i])
            else:
                resample_trajectory(t[i], ots[i], nts[j], output[j,i])
    return output


def reaction_rates_given_signal(response, signal, averaged=False):
    """Calculate the reaction rates at every timestamp of a response trajectory given a specified signal trajectory.

    Arguments:
        response {[type]} -- [description]
        signal {[type]} -- [description]

    Returns:
        [type] -- [description]
    """
    components = {}
    for name, comp in response['components'].items():
        # we want every component to have the shape (responses, signals, trajectory)
        components[name] = np.expand_dims(comp, axis=1)

    # resample the signal concentration and set the corresponding dictionary entries in components
    for name, straj in signal['components'].items():
        stime = signal['timestamps']
        rtime = response['timestamps']
        # shape (responses, signals, trajectories)
        components[name] = resample_trajectory_outer(
            straj, stime, rtime, averaged=averaged)

    return calculate_reaction_propensities(response['reactions'], components)


def log_likelihoods_given_signal(response, signal):
    """ Calculates the individual transition probabilities for a response trajectory given a signal.

    Returns an array of transition probabilities which when multiplied together result in the likelihood
    $p(X|S)$ where X is the response trajectory and S is the signal trajectory.
    """
    instantaneous_rates = reaction_rates_given_signal(response, signal)

    # since the chosen reaction depends on the response we swap the signal and response axes
    instantaneous_rates = np.swapaxes(instantaneous_rates, 0, 1)

    # do some very clever indexing
    num_responses = instantaneous_rates.shape[1]
    length = instantaneous_rates.shape[-1]
    idx = np.ix_(range(num_responses), range(length))
    rate_of_chosen_reaction = instantaneous_rates[:, idx[0], response['reaction_events'], idx[1]]

    time_delta = np.diff(response['timestamps'], prepend=0.0)
    averaged_rates = reaction_rates_given_signal(
        response, signal, averaged=True)
    total_averaged_rate = np.swapaxes(np.sum(averaged_rates, axis=-2), 0, 1)

    # shape (signals, responses, log_likelihoods)
    assert len(total_averaged_rate.shape) == 3

    # return the logarithm of `rate_of_chosen_reaction * np.exp(-total_averaged_rate * time_delta)`
    # also swap the axes back so that the shape is (responses, signals, trajectory)
    return np.swapaxes(np.log(rate_of_chosen_reaction) - total_averaged_rate * time_delta, 0, 1)


def generate_input(name, num_components, num_reactions, num_blocks, num_steps):
    return '{}\n{}\n{}\n{} {} {}\n{}\n'.format(name, num_components, num_reactions, num_blocks, 0, num_steps, 100)


def simulate_trajectory(input_name, output_name=None, trajectories=None, seed=4252):
    cmd = [os.path.abspath(str(settings.configuration['executable'])), input_name, '-s', str(seed)]
    if output_name is not None:
        cmd.extend(['-o', os.path.abspath(output_name)])

    if trajectories is not None:
        assert len(trajectories) >= 1
        for t in trajectories:
            cmd.extend(['-t', os.path.abspath(t)])

    print(' '.join(cmd), file=sys.stderr)
    subprocess.run(cmd, stdout=subprocess.DEVNULL, cwd=os.path.abspath(settings.configuration['gillespie_cwd']),
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


def simulation1():
    num_signals = 100
    num_responses_per_signal = 100

    for i_signal in range(num_signals):
        signal_name = "/data/signal/sig{}.traj".format(i_signal)
        save_signal("S", signal_name)
        print("FINISHED SIGNAL {}".format(signal_name))
        for i_response in range(num_responses_per_signal):
            simulate_trajectory('response', '/data/response/res{}-{}.traj'.format(i_signal, i_response), [signal_name],
                                seed=i_response)
