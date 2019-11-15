#!/usr/bin/env python3

import multiprocessing
from analyzer import analyzer, stochastic_sim, ornstein_uhlenbeck
from analyzer import configuration
import numba
from numba import njit
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

OUT_PATH = pathlib.Path(configuration.get()['output'])
num_signals = configuration.get()['num_signals']

SIGNAL_NETWORK, RESPONSE_NETWORK = configuration.read_reactions()


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
    timestamps = np.zeros((count, length))
    trajectory = np.zeros((count, 1, length), dtype=np.int16)
    reaction_events = np.zeros((count, length - 1), dtype=np.uint8)

    # initial values
    if initial_values is not None:
        trajectory[:, 0, 0] = initial_values
    else:
        raise RuntimeError('Need to specify initial values')

    stochastic_sim.simulate(
        timestamps, trajectory, reaction_events=reaction_events, reactions=SIGNAL_NETWORK)

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
        raise RuntimeError('Need to specify initial values')

    stochastic_sim.simulate(timestamps, trajectory, reaction_events=reaction_events,
                            reactions=RESPONSE_NETWORK, ext_components=signal_comps, ext_timestamps=signal_timestamps)

    return {
        'components': trajectory,
        'timestamps': timestamps,
        'reaction_events': reaction_events,
    }


@njit(parallel=True, fastmath=True, cache=True)
def log_evaluate_kde(points, dataset, inv_cov):
    d, n = dataset.shape
    num_r, _, m = points.shape
    assert len(inv_cov.shape) == 2

    log_norm_factor = -0.5 * np.log(np.linalg.det(inv_cov / (2 * np.pi)))

    result = np.empty((num_r, m))
    for j in numba.prange(num_r):
        tmp = np.zeros((n, m))
        for i in range(n):
            diff = np.empty((d, m))
            for k in range(m):
                for l in range(d):
                    diff[l, k] = dataset[l, i] - points[j, l, k]
            tdiff = np.dot(inv_cov, diff)
            tmp[i] = -np.sum(diff*tdiff, axis=0) / 2.0

        result[j] = analyzer.logsumexp(tmp) - np.log(n) - log_norm_factor

    return result


def calculate(i, num_responses, averaging_signals, kde_estimate, log_p0_signal):
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

    response_len = configuration.get()['response']['length']

    joined_distr = kde_estimate['joined']
    signal_distr = kde_estimate['signal']

    # generate responses from signals

    # first we sample the initial points from the joined distribution
    initial_comps = joined_distr.resample(num_responses)
    sig = generate_signals_sim(
        num_responses, length=response_len, initial_values=initial_comps[0])
    responses = generate_responses(
        num_responses, sig['timestamps'], sig['components'], length=response_len, initial_values=initial_comps[1])

    log_p_x_zero_this = joined_distr.logpdf(
        initial_comps) - signal_distr.logpdf(initial_comps[0])

    points = np.empty((num_responses, 2, num_signals))
    points[:, 0, :] = averaging_signals['components'][np.newaxis, :, 0, 0]
    points[:, 1, :] = responses['components'][:, 0, 0, np.newaxis]

    log_p_x_zero = log_evaluate_kde(
        points, joined_distr.dataset, joined_distr.inv_cov)
    log_p_x_zero -= log_p0_signal

    # TODO: make this configurable
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
                            response_timestamps, reaction_events, RESPONSE_NETWORK, out=mutual_information[1])
    mutual_information[1] += log_p_x_zero_this[:, np.newaxis]

    signal_components = averaging_signals['components']
    signal_timestamps = averaging_signals['timestamps']
    mutual_information[1] -= analyzer.log_averaged_likelihood(
        traj_lengths, signal_components, signal_timestamps, response_components, response_timestamps, reaction_events, RESPONSE_NETWORK, log_p_x_zero.T)

    np.save(OUT_PATH / 'mi.{}'.format(i),
            np.swapaxes(mutual_information, 0, 1))


def kde_estimate_p_0(size, traj_length, signal_init, response_init):
    sig = generate_signals_sim(
        size, length=traj_length, initial_values=signal_init)

    components = np.empty((size, 2), dtype=np.int16)
    components[:, 0] = sig['components'][:, 0, 0]
    components[:, 1] = response_init

    until = sig['timestamps'][:, -1]

    stochastic_sim.simulate_until(
        until, components, RESPONSE_NETWORK, ext_timestamps=sig['timestamps'], ext_components=sig['components'])

    np.save(OUT_PATH / 'equilibrated', components)

    return {
        'joined': gaussian_kde(components.T),
        'signal': gaussian_kde(components[:, 0])
    }


def main():
    OUT_PATH.mkdir(exist_ok=False)
    conf = configuration.get()

    runinfo = configuration.get()
    runinfo['run'] = {
        'started': datetime.now(timezone.utc)
    }

    with (OUT_PATH / 'info.toml').open('x') as f:
        toml.dump(runinfo, f)

    kde_estimate = kde_estimate_p_0(
        size=conf['kde_estimate']['size'],
        traj_length=conf['kde_estimate']['signal']['length'],
        signal_init=conf['kde_estimate']['signal']['initial'],
        response_init=conf['kde_estimate']['response']['initial']
    )
    print("generating signals...")
    initial_values = kde_estimate['signal'].resample(size=num_signals)
    signal_length = configuration.get()['signal']['length']
    combined_signal = generate_signals_sim(
        num_signals, length=signal_length, initial_values=initial_values)

    signal_distr = kde_estimate['signal']
    points = combined_signal['components'][np.newaxis, np.newaxis, :, 0, 0]
    log_p0_signal = log_evaluate_kde(
        points, signal_distr.dataset, signal_distr.inv_cov)

    num_responses = configuration.get()['num_responses']
    pbar = tqdm(total=num_responses, smoothing=0.9, desc='simulated responses')
    response_batch = multiprocessing.cpu_count()

    for i in range(num_responses // response_batch):
        calculate(i, response_batch, combined_signal,
                  kde_estimate, log_p0_signal)
        pbar.update(response_batch)

    i = num_responses // response_batch
    remaining_responses = num_responses % response_batch
    calculate(i, remaining_responses, combined_signal,
              kde_estimate, log_p0_signal)
    pbar.update(remaining_responses)

    runinfo['run']['ended'] = datetime.now(timezone.utc)
    with (OUT_PATH / 'info.toml').open('w') as f:
        toml.dump(runinfo, f)


if __name__ == '__main__':
    main()
