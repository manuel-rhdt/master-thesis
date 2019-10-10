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
    return x


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
        propensity = reaction['k'] * np.prod([components[x] for x in reaction['reactants']], axis=0)
        reaction_propensities.append(propensity)

    return np.choose(reaction_events, reaction_propensities)


# This function currently assumes that the first component of response is the original signal
#
# In the future we probably want to separate the signal from the responses
def likelihoods_of_reaction_events(signal, response):
    reaction_events = np.array(response['reaction_events'])
    concentrations = np.array(response['components'])
    timestamps = np.array(response['timestamps'])

    # integrate the signal trajectory and evaluate (interpolate) it at the timestamps of the response events
    signal_concentrations = np.array(signal['components'][0])
    signal_timestamps = np.array(signal['timestamps'])
    interpolated_signal = interpolate.interp1d(signal_timestamps, signal_concentrations, kind='previous',
                                               assume_sorted=True)[timestamps]
    assert len(signal_timestamps) == len(timestamps)
    integrated_signal = np.cumsum(signal_concentrations * np.diff(signal_timestamps, prepend=0.0))
    integrated_signal = np.interp(timestamps, signal_timestamps, integrated_signal)

    # for every reaction create an array that contains the propensity at every event
    reaction_propensities = []
    total_integrated_propensities = np.zeros_like(reaction_events, dtype='f8')
    for reaction in response['reactions']:
        # multiply reaction constant by the concentrations of the reactants
        propensity = reaction['k'] * np.prod([(concentrations[x] if x != 0 else interpolated_signal)
                                              for x in reaction['reactants']], axis=0)
        reaction_propensities.append(propensity)

        if 0 in reaction['reactants']:
            total_integrated_propensities += reaction['k'] * np.prod(
                [(concentrations[x] if x != 0 else integrated_signal)
                 for x in reaction['reactants']], axis=0)
        else:
            total_integrated_propensities += reaction['k'] * np.diff(timestamps, prepend=0.0) * np.prod(
                [concentrations[x]
                 for x in reaction['reactants']], axis=0)

    assert len(reaction_propensities) == len(response['reactions'])

    survival_probabilities = np.choose(reaction_events, reaction_propensities) * np.exp(-total_integrated_propensities)


def log_likelihood(trajectory):
    reaction_events = np.array(trajectory['reaction_events'])
    concentrations = np.array(trajectory['components'])

    assert reaction_events.shape[0] == concentrations.shape[1]

    # for every reaction create an array that contains the propensity at every event
    reaction_propensities = []
    for reaction in trajectory['reactions']:
        # multiply reaction constant by the concentrations of the reactants
        propensity = reaction['k'] * np.prod([concentrations[x]
                                              for x in reaction['reactants']], axis=0)
        reaction_propensities.append(propensity)

    chosen_reaction_propensity = np.choose(
        reaction_events, reaction_propensities)

    # since the random variates are sampled according to the survival probability the following association is correct
    survival_probabilities = np.array(trajectory['random_variates'])

    piecewise_probabilities = chosen_reaction_propensity * survival_probabilities

    return np.sum(np.log(piecewise_probabilities))


def generate_input(name, num_components, num_reactions, num_blocks, num_steps):
    return '{}\n{}\n{}\n{} {} {}\n{}\n'.format(name, num_components, num_reactions, num_blocks, 0, num_steps, 100)


def simulate_trajectory(input_name, output_name=None, seed=4252):
    cmd = [str(EXECUTABLE), input_name + '.inp', '-s', str(seed)]
    trajectory_path = input_name + '.traj'
    if output_name is not None:
        trajectory_path = output_name + '.traj'
        cmd.extend(['-o', trajectory_path])
    print(' '.join(cmd), file=sys.stderr)
    subprocess.run(cmd, stdout=subprocess.DEVNULL)
    trajectory = load_trajectory(pathlib.Path(trajectory_path))

    return trajectory


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
            trajectory['components'][component] = np.array(trajectory['components'][component])

        if 'random_variates' in trajectory:
            trajectory['random_variates'] = np.array(trajectory['random_variates'])

        if 'reaction_events' in trajectory:
            trajectory['reaction_events'] = np.array(trajectory['reaction_events'])
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
            pass


simulation = {
    "command": "repeat",
    "times": 50,
    "foreach": {
        "command": save_signal,
        "output": "/data/signal/sig{}.traj",
        "after": {

        }
    }
}
