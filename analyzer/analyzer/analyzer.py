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
from numba import jit

EXECUTABLE = os.getenv("GILLESPIE")


def ornstein_uhlenbeck_path(x0, t, mean_rev_speed, mean_rev_level, vola):
    """ Simulates a sample path for an Ornstein-Uhlenbeck process."""
    assert len(t) > 1
    x = scipy.stats.norm.rvs(size=len(t))
    x[0] = x0
    dt = np.diff(t)
    scale = std(dt, mean_rev_speed, vola)
    x[1:] = x[1:] * scale
    for i in range(1, len(x)):
        x[i] += mean(x[i - 1], dt[i - 1], mean_rev_speed, mean_rev_level)
    return np.abs(x)


def std(t, mean_rev_speed, vola):
    return np.sqrt(variance(t, mean_rev_speed, vola))


def variance(t, mean_rev_speed, vola):
    assert mean_rev_speed >= 0
    assert vola >= 0
    return vola * vola * (1.0 - np.exp(- 2.0 * mean_rev_speed * t)) / (2 * mean_rev_speed)


def mean(x0, t, mean_rev_speed, mean_rev_level):
    assert mean_rev_speed >= 0
    return x0 * np.exp(-mean_rev_speed * t) + (1.0 - np.exp(- mean_rev_speed * t)) * mean_rev_level


def generate_signal(name):
    timestamps = np.linspace(0, 100, 100000)
    x = ornstein_uhlenbeck_path(10000, timestamps, 0.001, 10000, 900)
    return {'timestamps': timestamps.tolist(), 'components': {name: x.tolist()}}


def save_signal(name, path):
    with open(str(path), 'wb') as file:
        msgpack.pack(generate_signal(name), file)


def calculate_reaction_propensities(reactions, components):
    """ Uses the provided information about reactions and the component counts at each time step to
    calculate an array of reaction propensities at each timestamp."""
    reaction_propensities = []
    for reaction in reactions:
        # multiply reaction constant by the concentrations of the reactants
        propensity = reaction['k'] * np.prod([components[x]
                                              for x in reaction['reactants']], axis=0)
        reaction_propensities.append(propensity)

    return reaction_propensities


@jit(nopython=True)
def resample_trajectory(trajectory, old_timestamps, new_timestamps):
    """ Resample trajectory with `old_timestamp` at the times in `new_timestamps`.

    This function assumes that `new_timestamps` is ordered.

    Arguments:
        trajectory {[type]} -- [description]
        old_timestamps {[type]} -- [description]
        new_timestamps {[type]} -- [description]

    Returns:
        [type] -- an array of values
    """
    trajectory = np.asarray(trajectory)

    resampled_trajectory = np.zeros_like(
        new_timestamps, dtype=trajectory.dtype)

    new_idx = 0
    for old_idx, ts in enumerate(old_timestamps):
        while new_timestamps[new_idx] < ts and new_idx < len(resampled_trajectory):
            resampled_trajectory[new_idx] = trajectory[max(old_idx - 1, 0)]
            new_idx += 1
        if new_idx >= len(resampled_trajectory):
            break

    while new_idx < len(resampled_trajectory):
        resampled_trajectory[new_idx] = trajectory[-1]
        new_idx += 1

    return resampled_trajectory


def reaction_rates_given_signal(response, signal):
    """Calculate the reaction rates at every timestamp of a response trajectory given a specified signal trajectory.

    Arguments:
        response {[type]} -- [description]
        signal {[type]} -- [description]

    Returns:
        [type] -- [description]
    """
    components = response['components']

    # resample the signal concentration and set the corresponding dictionary entries in components
    for comp, values in signal['components'].items():
        components[comp] = resample_trajectory(
            values, signal['timestamps'], response['timestamps'])

    return calculate_reaction_propensities(response['reactions'], components)


@jit(nopython=True)
def resample_averaged_trajectory(trajectory, old_timestamps, new_timestamps):
    """ Resample and averagte trajectory with `old_timestamp` at the times in `new_timestamps`.

    This function assumes that `new_timestamps` is ordered.

    Arguments:
        trajectory {[type]} -- [description]
        old_timestamps {[type]} -- [description]
        new_timestamps {[type]} -- [description]

    Returns:
        [type] -- an array of values
    """
    trajectory = np.asarray(trajectory)

    resampled_trajectory = np.zeros_like(
        new_timestamps, dtype=trajectory.dtype)

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

        resampled_trajectory[new_idx + 1] = acc / delta_t

    return resampled_trajectory


def mean_reaction_rates_given_signal(response, signal):
    """Similar to the function above. However here we don't calculate the immediate rates at every timestamp but
    rather the averaged reaction rates between two timestamps.

    Arguments:
        response {[type]} -- [description]
        signal {[type]} -- [description]
    """
    components = response['components']

    for comp, values in signal['components'].items():
        components[comp] = resample_averaged_trajectory(
            values, signal['timestamps'], response['timestamps'])

    return calculate_reaction_propensities(response['reactions'], components)


def likelihoods_given_signal(response, signal):
    """ Calculates the individual transition probabilities for a response trajectory given a signal.

    Returns an array of transition probabilities which when multiplied together result in the likelihood
    $p(X|S)$ where X is the response trajectory and S is the signal trajectory.
    """
    instantaneous_rates = reaction_rates_given_signal(response, signal)
    rate_of_chosen_reaction = np.choose(
        response['reaction_events'], instantaneous_rates)

    time_deltas = np.diff(response['timestamps'], prepend=0.0)
    averaged_rates = mean_reaction_rates_given_signal(response, signal)
    total_averaged_rate = np.sum(averaged_rates, axis=0)

    return rate_of_chosen_reaction * np.exp(-total_averaged_rate * time_deltas)


def generate_input(name, num_components, num_reactions, num_blocks, num_steps):
    return '{}\n{}\n{}\n{} {} {}\n{}\n'.format(name, num_components, num_reactions, num_blocks, 0, num_steps, 100)


def simulate_trajectory(input_name, output_name=None, trajectories=None, seed=4252):
    cmd = [str(EXECUTABLE), input_name + '.inp', '-s', str(seed)]
    if output_name is not None:
        cmd.extend(['-o', output_name])

    if trajectories is not None:
        assert len(trajectories) >= 1
        for t in trajectories:
            cmd.extend(['-t', t])

    print(' '.join(cmd), file=sys.stderr)
    subprocess.run(cmd, stdout=subprocess.DEVNULL,
                   timeout=10.0).check_returncode()


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


def simulate():
    likelihoods = []

    count = multiprocessing.cpu_count()
    pool = multiprocessing.Pool(processes=count)

    results = pool.map(calculate, range(50))

    print("save to likelihoods.txt")
    np.savetxt('likelihoods.txt', np.array(results))
    print(results)


def simulation1():
    num_signals = 10
    num_responses_per_signal = 10

    for i_signal in range(num_signals):
        signal_name = "/data/signal/sig{}.traj".format(i_signal)
        save_signal("S", signal_name)
        for i_response in range(num_responses_per_signal):
            simulate_trajectory('response', '/data/response/res{}-{}.traj'.format(i_signal, i_response), [signal_name],
                                seed=i_response)
