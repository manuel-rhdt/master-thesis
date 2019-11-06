#!/usr/bin/env python3

import settings
from analyzer import analyzer, stochastic_sim
import pathlib
import numpy as np
import os
import multiprocessing
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


reaction_k = np.array([0.005, 0.02])
reaction_reactants = np.array([[0], [1]])
reaction_products = np.array([[0, 1], [-1, -1]])


length = 100000
timestamps = np.zeros((num_responses, length))
trajectory = np.zeros((num_responses, 1, length))
reaction_events = np.zeros((num_responses, length - 1), dtype='i4')


def generate_responses(s, num_responses, combined_signal):
    trajectory[..., 0, 0] = 1000.0

    ext_components = combined_signal['components']['S'][[s]]
    ext_timestamps = combined_signal['timestamps'][s]

    stochastic_sim.simulate(timestamps, trajectory, ext_components, ext_timestamps,
                            reaction_k, reaction_reactants, reaction_products, reaction_events)

    return {
        'components': trajectory,
        'timestamps': timestamps,
        'reaction_events': reaction_events,
    }


def calculate(s, combined_signal):
    combined_response = generate_responses(s, num_responses, combined_signal)
    response_len = combined_response['timestamps'][0].shape[-1]

    granularity = 100
    result_size = -(-(response_len - 1) // granularity)

    mutual_information = np.empty(
        (2, num_responses, result_size))

    mutual_information[0] = combined_response['timestamps'][:, 1::granularity]

    signal_components = np.expand_dims(
        combined_signal['components']['S'][[s]], axis=1)
    signal_timestamps = combined_signal['timestamps'][[s]]
    response_components = combined_response['components']
    response_timestamps = combined_response['timestamps']
    reaction_events = combined_response['reaction_events']

    mutual_information[1] = analyzer.log_averaged_likelihood(
        signal_components, signal_timestamps, response_components, response_timestamps, reaction_k, reaction_reactants, reaction_events, granularity=granularity)

    signal_components = np.expand_dims(np.delete(
        combined_signal['components']['S'], s, axis=0), axis=1)
    signal_timestamps = np.delete(combined_signal['timestamps'], s, axis=0)
    mutual_information[1] -= analyzer.log_averaged_likelihood(
        signal_components, signal_timestamps, response_components, response_timestamps, reaction_k, reaction_reactants, reaction_events, granularity=granularity)

    name = CONFIGURATION['output']
    np.save(os.path.expandvars('{}.{}'.format(name, s)),
            np.swapaxes(mutual_information, 0, 1))


def generate_signals():
    num_signals = CONFIGURATION['num_signals']

    combined_signal = None
    for s in range(num_signals):
        signal = analyzer.generate_signal(
            'S', duration, 1/resolution, x0=mean, mean=mean, correlation_time=corr_time, diffusion=diffusion)
        if combined_signal is None:
            length = len(signal['timestamps'])
            names = ['S']
            combined_signal = analyzer.empty_trajectory(
                names, num_signals, length)

        combined_signal['components']['S'][s] = signal['components']['S']
        combined_signal['timestamps'][s] = signal['timestamps']

    return combined_signal


def main():
    num_signals = CONFIGURATION['num_signals']

    pathlib.Path(os.path.expandvars(
        CONFIGURATION['output'])).parent.mkdir(exist_ok=True)

    combined_signal = generate_signals()

    pbar = tqdm(total=num_signals)
    for s in range(num_signals):
        calculate(s, combined_signal)
        pbar.update()


def profile():
    num_signals = CONFIGURATION['num_signals']
    pathlib.Path(os.path.expandvars(
        CONFIGURATION['output'])).parent.mkdir(exist_ok=True)
    combined_signal = generate_signals()

    for s in range(1):
        calculate(s, combined_signal)


if __name__ == '__main__':
    main()
