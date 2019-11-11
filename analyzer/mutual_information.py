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
from scipy.stats import gaussian_kde


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
    trajectory = np.zeros((count, 1, length))
    reaction_events = np.zeros((count, length - 1), dtype='i1')

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


def generate_histogram(signal_start, size=80, traj_length=5000):
    signals = generate_signals_sim(
        size, length=traj_length, initial_values=signal_start)

    # reverse the signals
    sig_t = -signals['timestamps'][..., ::-1]
    sig_t -= sig_t[..., [0]]
    sig_comp = signals['components'][..., ::-1]

    components = np.empty((size, 2), dtype=np.int32)
    components[:, 0] = sig_comp[:, 0, 0]
    components[:, 1] = 1000

    until = sig_t[:, -1]

    stochastic_sim.simulate_until(
        until, components, reactions, ext_timestamps=sig_t, ext_components=sig_comp)

    return gaussian_kde(components[:, 1])


def generate_histograms(num_signals, progress=False):
    signal_t0 = np.random.normal(
        size=num_signals, loc=mean_s, scale=np.sqrt(mean_s))

    if progress:
        pbar = tqdm(total=num_signals, smoothing=0.9,
                    desc='generated histograms')

    histograms = []
    for s0 in signal_t0:
        hist = generate_histogram(s0)
        histograms.append(hist)
        if progress:
            pbar.update()

    if progress:
        pbar.close()

    return signal_t0, histograms


def evaluate_hist(x, hist):
    return hist(x)
    # density, bins = hist
    # indices = np.digitize(x, bins)
    # out_of_bounds = np.logical_or(indices == 0, indices == len(bins))
    # conditions = [out_of_bounds, np.logical_not(out_of_bounds)]
    # choices = [0.0, np.take(density, indices, mode='wrap')]
    # return np.select(conditions, choices)


def sample_hist(hist, size=1):
    return hist.resample(size=size)
    # density, bins = hist
    # cdf = np.cumsum(density)
    # cdf = cdf / cdf[-1]
    # rv = np.random.random_sample(size=size)
    # indices = np.digitize(rv, bins)


def calculate(i, num_responses, averaging_signals, histograms):
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
    signal_init, new_hists = generate_histograms(num_responses)
    sig = generate_signals_sim(
        num_responses, length=response_len, initial_values=signal_init)
    responses = generate_responses(
        num_responses, sig['timestamps'], sig['components'], length=response_len, initial_values=np.concatenate([hist.resample(size=1).reshape((1)) for hist in new_hists]))

    log_p_x_zero_this = np.empty(num_responses)
    for r, hist in enumerate(new_hists):
        log_p_x_zero_this[r] = hist.logpdf(responses['components'][r, 0, 0])

    log_p_x_zero = np.empty((num_signals, num_responses))
    for s, hist in enumerate(histograms):
        log_p_x_zero[s] = hist.logpdf(responses['components'][:, 0, 0])

    result_size = response_len - 1
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
    mutual_information[1] += log_p_x_zero_this[:, np.newaxis]

    signal_components = averaging_signals['components']
    signal_timestamps = averaging_signals['timestamps']
    mutual_information[1] -= analyzer.log_averaged_likelihood(
        signal_components, signal_timestamps, response_components, response_timestamps, reaction_events, reactions, log_p_x_zero)

    np.save(os.path.join(OUT_PATH, 'mi.{}'.format(i)),
            np.swapaxes(mutual_information, 0, 1))


def kde_estimate_p_0(size=500, traj_length=5000):
    signal_init = np.random.normal(size=size, loc=mean_s, scale=np.sqrt(mean_s))
    sig = generate_signals_sim(
        size, length=traj_length, initial_values=signal_init)
    res = generate_responses(
        size, sig['timestamps'], sig['components'], length=traj_length, initial_values=mean_x)

    points = np.array([sig['components'][:,0,-1], res['components'][:,0,-1]])
    return gaussian_kde(points)


def main():
    pathlib.Path(OUT_PATH).mkdir(exist_ok=True)

    init_vals, histograms = generate_histograms(num_signals, progress=True)
    print("generating signals...")
    combined_signal = generate_signals_sim(
        num_signals, length=50000, initial_values=init_vals)

    num_responses = int(CONFIGURATION['num_responses'])
    pbar = tqdm(total=num_responses, smoothing=0.9, desc='simulated responses')
    response_batch = multiprocessing.cpu_count()

    for i in range(num_responses // response_batch):
        calculate(i, response_batch, combined_signal, histograms)
        pbar.update(response_batch)
    i = num_responses // response_batch
    remaining_responses = num_responses % response_batch
    calculate(i, remaining_responses, combined_signal, histograms)
    pbar.update(remaining_responses)


if __name__ == '__main__':
    main()
