#!/usr/bin/env python3

import multiprocessing
import settings
from analyzer import analyzer, stochastic_sim, ornstein_uhlenbeck
import pathlib
import numpy as np
import os
import queue
import math
from tqdm import tqdm


CONFIGURATION = settings.configuration['mutual_information']
mean = CONFIGURATION['signal_mean']
duration = CONFIGURATION['signal_duration']
corr_time = CONFIGURATION['signal_correlation_time']
diffusion = CONFIGURATION['signal_diffusion']
resolution = CONFIGURATION['signal_resolution']

num_signals = CONFIGURATION['num_signals']
num_responses = CONFIGURATION['num_responses']

kappa = 20.0
lamda = 0.005
rho = 0.005
mu = 0.02
mean_s = kappa / lamda
mean_x = mean_s * rho / mu

reaction_k = np.array([rho, mu])
reaction_reactants = np.array([[0], [1]])
reaction_products = np.array([[0, 1], [-1, -1]])


def generate_signals_numerical(n_signals):
    stimestamps = np.arange(0, duration, 1/resolution)
    signal_c = np.empty((n_signals, 1, len(stimestamps)))
    for i in range(n_signals):
        signal_c[i][0] = np.clip(ornstein_uhlenbeck.generate(
            stimestamps, x0=mean, correlation_time=corr_time, diffusion_constant=diffusion, mean=mean), 0.0, None)
    return {
        'timestamps': np.broadcast_to(stimestamps, (n_signals, len(stimestamps))),
        'components': signal_c
    }


def generate_signals_sim(n_signals, length=100000):
    k = np.array([kappa, lamda])
    reactants = np.array([[-1], [0]])
    products = np.array([[0], [-1]])

    timestamps = np.zeros((n_signals, length))
    trajectory = np.zeros((n_signals, 1, length))
    reaction_events = np.zeros((n_signals, length - 1), dtype='i1')

    # initial values
    for s in range(n_signals):
        trajectory[s, 0, 0] = np.random.normal(
            loc=mean_s, scale=np.sqrt(mean_s))

    stochastic_sim.simulate(timestamps, trajectory, ext_components=None, ext_timestamps=None,
                            reaction_k=k, reaction_reactants=reactants, reaction_products=products, reaction_events=reaction_events)

    return {
        'timestamps': timestamps,
        'components': trajectory,
        'reaction_events': reaction_events
    }


def generate_responses(num_responses, signal_timestamps, signal_comps, length=100000):
    timestamps = np.zeros((num_responses, length))
    trajectory = np.zeros((num_responses, 1, length))
    reaction_events = np.zeros((num_responses, length - 1), dtype='i1')

    # initial values
    for r in range(num_responses):
        trajectory[r, 0, 0] = np.random.normal(
            loc=mean_x, scale=np.sqrt(mean_x * (1.0 + rho/(lamda + mu))))

    stochastic_sim.simulate(timestamps, trajectory, ext_components=signal_comps, ext_timestamps=signal_timestamps,
                            reaction_k=reaction_k, reaction_reactants=reaction_reactants, reaction_products=reaction_products, reaction_events=reaction_events)

    return {
        'components': trajectory,
        'timestamps': timestamps,
        'reaction_events': reaction_events,
    }


def calculate(i, num_responses, combined_signal):
    if num_responses == 0:
        return

    # generate responses from signals
    sig = generate_signals_sim(num_responses, length=50000)
    combined_response = generate_responses(num_responses,
                                           sig['timestamps'], sig['components'], length=50000)
    response_len = combined_response['timestamps'][0].shape[-1]

    result_size = (response_len - 1)

    # we create an array with the following dimensions to hold the results of our calculations
    # dimension1 = 2: mutual_information[0] holds the trajectory length times
    #                 mutual_information[1] holds the mutual information values
    # dimension2: arrays of responses
    # dimension3: arrays of trajectories
    mutual_information = np.empty(
        (2, num_responses, result_size))

    # store the trajectory lengths for which the mutual information is computed
    # we don't store the initial timestamp since this does not contain any information
    mutual_information[0] = combined_response['timestamps'][:, 1:]

    response_components = combined_response['components']
    response_timestamps = combined_response['timestamps']
    reaction_events = combined_response['reaction_events']

    analyzer.log_likelihood(
        sig['components'], sig['timestamps'], response_components, response_timestamps, reaction_k, reaction_reactants, reaction_events, out=mutual_information[1])

    signal_components = combined_signal['components']
    signal_timestamps = combined_signal['timestamps']
    mutual_information[1] -= analyzer.log_averaged_likelihood(
        signal_components, signal_timestamps, response_components, response_timestamps, reaction_k, reaction_reactants, reaction_events)

    name = CONFIGURATION['output']

    np.save(os.path.expandvars('{}.{}'.format(name, i)),
            np.swapaxes(mutual_information, 0, 1))


def main():
    num_signals = CONFIGURATION['num_signals']

    pathlib.Path(os.path.expandvars(
        CONFIGURATION['output'])).parent.mkdir(exist_ok=True)

    combined_signal = generate_signals_sim(num_signals)

    pbar = tqdm(total=num_responses)
    response_batch = multiprocessing.cpu_count()

    for i in range(num_responses // response_batch):
        calculate(i, response_batch, combined_signal)
        pbar.update(response_batch)
    i = num_responses // response_batch + 1
    remaining_responses = num_responses % response_batch
    calculate(i, remaining_responses, combined_signal)
    pbar.update(remaining_responses)


if __name__ == '__main__':
    main()
