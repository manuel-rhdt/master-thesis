#!/usr/bin/env python3

import multiprocessing
from analyzer import analyzer, stochastic_sim, ornstein_uhlenbeck
from analyzer import configuration
import pathlib
import numpy as np
import os
import queue
import math
from tqdm import tqdm
import scipy
from scipy.stats import gaussian_kde
from datetime import datetime, timezone
import toml


CONFIGURATION = configuration.get()['mutual_information']

OUT_PATH = CONFIGURATION['output']
num_signals = CONFIGURATION['num_signals']

kappa = 20.0
lamda = 0.005
rho = 0.005
mu = 0.02
mean_s = kappa / lamda
mean_x = mean_s * rho / mu

reactions = stochastic_sim.ReactionNetwork(2)
reactions.k = np.array([rho, mu], dtype=np.single)
reactions.reactants = np.array([[0], [1]], dtype=np.int32)
reactions.products = np.array([[0, 1], [-1, -1]], dtype=np.int32)


# def generate_signals_numerical(count):
#     stimestamps = np.arange(0, duration, 1/resolution)
#     signal_c = np.empty((count, 1, len(stimestamps)))
#     for i in range(count):
#         signal_c[i][0] = np.clip(ornstein_uhlenbeck.generate(
#             stimestamps, x0=mean, correlation_time=corr_time, diffusion_constant=diffusion, mean=mean), 0.0, None)
#     return {
#         'timestamps': np.broadcast_to(stimestamps, (count, len(stimestamps))),
#         'components': signal_c
#     }


def generate_signals_sim(count, length=100000, initial_values=None):
    sig_reactions = stochastic_sim.ReactionNetwork(2)
    sig_reactions.k = np.array([kappa, lamda], dtype=np.single)
    sig_reactions.reactants = np.array([[-1], [0]], dtype=np.int32)
    sig_reactions.products = np.array([[0], [-1]], dtype=np.int32)

    timestamps = np.zeros((count, length))
    trajectory = np.zeros((count, 1, length), dtype=np.int16)
    reaction_events = np.zeros((count, length - 1), dtype=np.uint8)

    # initial values
    if initial_values is not None:
        trajectory[:, 0, 0] = initial_values
    else:
        trajectory[:, 0, 0] = np.random.normal(
            size=count, loc=mean_s, scale=np.sqrt(mean_s))

    stochastic_sim.simulate(
        timestamps, trajectory, reaction_events=reaction_events, reactions=sig_reactions)

    return {
        'timestamps': timestamps,
        'components': trajectory,
        'reaction_events': reaction_events
    }


def generate_responses(count, signal_timestamps, signal_comps, length=100000, initial_values=None):
    timestamps = np.zeros((count, length))
    trajectory = np.zeros((count, 1, length), dtype=np.int16)
    reaction_events = np.zeros((count, length - 1), dtype=np.uint8)

    # initial values
    if initial_values is not None:
        trajectory[:, 0, 0] = initial_values
    else:
        trajectory[:, 0, 0] = np.random.normal(
            size=count, loc=mean_x, scale=np.sqrt(mean_x * (1.0 + rho/(lamda + mu))))

    stochastic_sim.simulate(timestamps, trajectory, reaction_events=reaction_events,
                            reactions=reactions, ext_components=signal_comps, ext_timestamps=signal_timestamps)

    return {
        'components': trajectory,
        'timestamps': timestamps,
        'reaction_events': reaction_events,
    }


def calculate(i, num_responses, averaging_signals, kde_estimate):
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

    signal_distr = scipy.stats.norm(loc=mean_s, scale=np.sqrt(mean_s))

    # generate responses from signals
    initial_comps = kde_estimate.resample(num_responses)
    sig = generate_signals_sim(
        num_responses, length=response_len, initial_values=initial_comps[0])
    responses = generate_responses(
        num_responses, sig['timestamps'], sig['components'], length=response_len, initial_values=initial_comps[1])

    log_p_x_zero_this = kde_estimate.logpdf(
        initial_comps) - signal_distr.logpdf(initial_comps[0])

    values = np.empty((2, num_signals, num_responses))
    values[0, :, :] = averaging_signals['components'][:, 0, 0, np.newaxis]
    values[1, :, :] = responses['components'][:, 0, 0]
    values = np.reshape(values, (2, -1))

    log_p_x_zero = kde_estimate.logpdf(
        values) - signal_distr.logpdf(values[0])
    log_p_x_zero = np.reshape(log_p_x_zero, (num_signals, num_responses))

    result_size = 5000
    traj_lengths = np.geomspace(0.01, 1000, num=result_size, dtype=np.single)
    # we create an array with the following dimensions to hold the results of our calculations
    # dimension1 = 2: mutual_information[0] holds the trajectory length times
    #                 mutual_information[1] holds the mutual information values
    # dimension2: arrays of responses
    # dimension3: arrays of trajectories
    mutual_information = np.empty(
        (2, num_responses, result_size), dtype=np.single)

    # store the trajectory lengths for which the mutual information is computed
    mutual_information[0] = traj_lengths

    response_components = responses['components']
    response_timestamps = responses['timestamps']
    reaction_events = responses['reaction_events']

    analyzer.log_likelihood(traj_lengths, sig['components'], sig['timestamps'], response_components,
                            response_timestamps, reaction_events, reactions, out=mutual_information[1])
    mutual_information[1] += log_p_x_zero_this[:, np.newaxis]

    signal_components = averaging_signals['components']
    signal_timestamps = averaging_signals['timestamps']
    mutual_information[1] -= analyzer.log_averaged_likelihood(
        traj_lengths, signal_components, signal_timestamps, response_components, response_timestamps, reaction_events, reactions, log_p_x_zero)

    np.save(os.path.join(OUT_PATH, 'mi.{}'.format(i)),
            np.swapaxes(mutual_information, 0, 1))


def kde_estimate_p_0(size=500, traj_length=5000):
    signal_init = np.random.normal(
        size=size, loc=mean_s, scale=np.sqrt(mean_s))
    sig = generate_signals_sim(
        size, length=traj_length, initial_values=signal_init)

    components = np.empty((size, 2), dtype=np.int16)
    components[:, 0] = sig['components'][:, 0, 0]
    components[:, 1] = mean_x

    until = sig['timestamps'][:, -1]

    stochastic_sim.simulate_until(
        until, components, reactions, ext_timestamps=sig['timestamps'], ext_components=sig['components'])

    return gaussian_kde(components.T)


def main():
    output_path = pathlib.Path(OUT_PATH)
    output_path.mkdir(exist_ok=False)

    runinfo = configuration.get()
    runinfo['run'] = {
        'started': datetime.now(timezone.utc)
    }

    with (output_path / 'info.toml').open('w') as f:
        toml.dump(runinfo, f)

    kde_estimate = kde_estimate_p_0(size=num_signals)
    print("generating signals...")
    combined_signal = generate_signals_sim(num_signals, length=50000)

    num_responses = CONFIGURATION['num_responses']
    pbar = tqdm(total=num_responses, smoothing=0.9, desc='simulated responses')
    response_batch = multiprocessing.cpu_count()

    for i in range(num_responses // response_batch):
        calculate(i, response_batch, combined_signal, kde_estimate)
        pbar.update(response_batch)
    i = num_responses // response_batch
    remaining_responses = num_responses % response_batch
    calculate(i, remaining_responses, combined_signal, kde_estimate)
    pbar.update(remaining_responses)


if __name__ == '__main__':
    main()
