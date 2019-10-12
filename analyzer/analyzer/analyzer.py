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

EXECUTABLE = pathlib.Path(os.getenv("GILLESPIE"))


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


def calculate_reaction_propensities(reaction_events, reactions, components):
    """ Uses the provided information about reaction_events, reactions and the component counts at each time step to
    calculate an array of reaction propensities at each timestamp."""
    reaction_propensities = []
    for reaction in reactions:
        # multiply reaction constant by the concentrations of the reactants
        propensity = reaction['k'] * np.prod([components[x]
                                              for x in reaction['reactants']], axis=0)
        reaction_propensities.append(propensity)

    return np.choose(reaction_events, reaction_propensities)


# This function currently assumes that the first component of response is the original signal
#
# In the future we probably want to separate the signal from the responses
def likelihoods_of_reaction_events(signal, response):
    components = response['components']
    time_deltas = np.diff(response['timestamps'], prepend=0.0)

    # resample the signal concentration
    for comp, values in signal['components'].items():
        interpolation = interpolate.interp1d(signal['timestamps'], values, kind='previous',
                                             fill_value=(values[0], values[-1]), bounds_error=False, assume_sorted=True)
        components[comp] = interpolation(response['timestamps'])

    propensities = calculate_reaction_propensities(
        response['reaction_events'], response['reactions'], components)

    # to calculate the survival probability for each time delta we need to integrate the signal and then interpolate
    # it.
    integrated_signal = np.cumsum(
        signal['components']['S'] * np.diff(signal['timestamps'], prepend=0.0))
    integrated_signal = np.interp(
        response['timestamps'], signal['timestamps'], integrated_signal) / time_deltas
    components['S'] = integrated_signal

    total_integrated_propensities = np.zeros_like(
        response['reaction_events'], dtype='f8')
    for reaction in response['reactions']:
        # multiply reaction constant by the concentrations of the reactants
        propensity = reaction['k'] * np.prod([components[x]
                                              for x in reaction['reactants']], axis=0)
        total_integrated_propensities += propensity

    return propensities * np.exp(-total_integrated_propensities * time_deltas)


def single_trajectory_piecewise_likelihoods(trajectory):
    reaction_events = np.array(trajectory['reaction_events'])
    concentrations = np.array(trajectory['components'])

    assert reaction_events.shape[0] == concentrations.shape[1]

    chosen_reaction_propensity = calculate_reaction_propensities(reaction_events, trajectory['reactions'],
                                                                 concentrations)

    # since the random variates are sampled according to the survival probability the following association is correct
    survival_probabilities = np.array(trajectory['random_variates'])

    piecewise_probabilities = chosen_reaction_propensity * survival_probabilities

    return piecewise_probabilities


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


def calculate(trajectory_num):
    trajectory = simulate_trajectory('response', output_name='/data/response' + str(trajectory_num),
                                     seed=trajectory_num)
    return log_likelihood(trajectory)


def load_trajectories(path, glob):
    trajectories = []
    for file_path in pathlib.Path(path).glob(glob):
        trajectories.append(load_trajectory(file_path))
    return trajectories


def load_trajectory(path):
    print('loading ' + str(path))
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


def get_histogram(path, glob, equil_time):
    """
    DELETE THHIS
    :param path:
    :param glob:
    :param equil_time:
    :return:
    """
    datapoints = []
    for file in path.glob(glob):
        print('processing ' + str(file))
        with file.open() as f:
            trajectory = json.load(f)
            cutoff_index = np.searchsorted(
                trajectory['timestamps'], equil_time)
            for k, comp in enumerate(trajectory['components']):
                if len(datapoints) == k:
                    datapoints.append([])
                datapoints[k].append(comp[cutoff_index])
    return np.array(datapoints)


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
