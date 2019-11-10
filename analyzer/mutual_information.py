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

OUT_PATH = os.path.expandvars(CONFIGURATION['output'])

mean = CONFIGURATION['signal_mean']
duration = CONFIGURATION['signal_duration']
corr_time = CONFIGURATION['signal_correlation_time']
diffusion = CONFIGURATION['signal_diffusion']
resolution = CONFIGURATION['signal_resolution']

num_signals = CONFIGURATION['num_signals']

kappa = 20.0
lamda = 0.005
rho = 0.005
mu = 0.02
mean_s = kappa / lamda
mean_x = mean_s * rho / mu

reactions = stochastic_sim.ReactionNetwork(2)
reactions.k = np.array([rho, mu])
reactions.reactants = np.array([[0], [1]], dtype=np.int32)
reactions.products = np.array([[0, 1], [-1, -1]], dtype=np.int32)


def generate_signals_numerical(count):
    stimestamps = np.arange(0, duration, 1/resolution)
    signal_c = np.empty((count, 1, len(stimestamps)))
    for i in range(count):
        signal_c[i][0] = np.clip(ornstein_uhlenbeck.generate(
            stimestamps, x0=mean, correlation_time=corr_time, diffusion_constant=diffusion, mean=mean), 0.0, None)
    return {
        'timestamps': np.broadcast_to(stimestamps, (count, len(stimestamps))),
        'components': signal_c
    }


def generate_signals_sim(count, length=100000, initial_values=None):
    sig_reactions = stochastic_sim.ReactionNetwork(2)
    sig_reactions.k = np.array([kappa, lamda])
    sig_reactions.reactants = np.array([[-1], [0]], dtype=np.int32)
    sig_reactions.products = np.array([[0], [-1]], dtype=np.int32)

    timestamps = np.zeros((count, length))
    trajectory = np.zeros((count, 1, length))
    reaction_events = np.zeros((count, length - 1), dtype='i1')

    # initial values
    if initial_values is not None:
        trajectory[:, 0, 0] = initial_values
    else:
        for s in range(count):
            trajectory[s, 0, 0] = np.random.normal(
                loc=mean_s, scale=np.sqrt(mean_s))

    stochastic_sim.simulate(
        timestamps, trajectory, reaction_events=reaction_events, reactions=sig_reactions)

    return {
        'timestamps': timestamps,
        'components': trajectory,
        'reaction_events': reaction_events
    }


def generate_responses(count, signal_timestamps, signal_comps, length=100000):
    timestamps = np.zeros((count, length))
    trajectory = np.zeros((count, 1, length))
    reaction_events = np.zeros((count, length - 1), dtype='i1')

    # initial values
    for r in range(count):
        trajectory[r, 0, 0] = np.random.normal(
            loc=mean_x, scale=np.sqrt(mean_x * (1.0 + rho/(lamda + mu))))

    stochastic_sim.simulate(timestamps, trajectory, reaction_events=reaction_events,
                            reactions=reactions, ext_components=signal_comps, ext_timestamps=signal_timestamps)

    return {
        'components': trajectory,
        'timestamps': timestamps,
        'reaction_events': reaction_events,
    }


def generate_histogram(signal_start, size=100, traj_length=5000):
    signals = generate_signals_sim(
        size, length=traj_length, initial_values=signal_start)

    # reverse the signals
    sig_t = -signals['timestamps'][..., ::-1]
    sig_t += sig_t[..., [0]]
    sig_comp = signals['components'][..., ::-1]

    stochastic_sim.simulate_until()
    responses = generate_responses(size, sig_t, sig_comp, length=traj_length)


def calculate(i, num_responses, averaging_signals):
    """ Calculates and stores the mutual information for `num_responses` respones.

    This function does the following:
    1. generate signals
    2. generate responses following the signals
    3. calculate the log likelihoods of the responses
    4. calculate the log marginal probabilities of the responses
    5. use both quantities to estimate the mutual information
    """
    if num_responses == 0:
        return

    response_len = 50000

    # generate responses from signals
    sig = generate_signals_sim(num_responses, length=response_len)
    responses = generate_responses(
        num_responses, sig['timestamps'], sig['components'], length=response_len)

    result_size = (response_len - 1)

    # we create an array with the following dimensions to hold the results of our calculations
    # dimension1 = 2: mutual_information[0] holds the trajectory length times
    #                 mutual_information[1] holds the mutual information values
    # dimension2: arrays of responses
    # dimension3: arrays of trajectories
    mutual_information = np.empty(
        (2, num_responses, result_size))

    # store the trajectory lengths for which the mutual information is computed
    mutual_information[0] = responses['timestamps'][:, 1:] -\
        responses['timestamps'][:, [0]]

    response_components = responses['components']
    response_timestamps = responses['timestamps']
    reaction_events = responses['reaction_events']

    analyzer.log_likelihood(sig['components'], sig['timestamps'], response_components,
                            response_timestamps, reaction_events, reactions, out=mutual_information[1])

    signal_components = averaging_signals['components']
    signal_timestamps = averaging_signals['timestamps']
    mutual_information[1] -= analyzer.log_averaged_likelihood(
        signal_components, signal_timestamps, response_components, response_timestamps, reaction_events, reactions)

    np.save(os.path.join(OUT_PATH, 'mi.{}'.format(i)),
            np.swapaxes(mutual_information, 0, 1))


def main():
    pathlib.Path(OUT_PATH).mkdir(exist_ok=True)

    combined_signal = generate_signals_sim(num_signals, length=50000)

    num_responses = int(CONFIGURATION['num_responses'])
    pbar = tqdm(total=num_responses, smoothing=0.9, desc='simulated responses')
    response_batch = multiprocessing.cpu_count()

    for i in range(num_responses // response_batch):
        calculate(i, response_batch, combined_signal)
        pbar.update(response_batch)
    i = num_responses // response_batch
    remaining_responses = num_responses % response_batch
    calculate(i, remaining_responses, combined_signal)
    pbar.update(remaining_responses)


if __name__ == '__main__':
    main()
