#!/usr/bin/env python3
import settings
from analyzer import analyzer
import pathlib
import numpy as np
import os
import multiprocessing
import math


CONFIGURATION = settings.configuration['mutual_information']
mean = CONFIGURATION['signal_mean']
corr_time = CONFIGURATION['signal_correlation_time']
diffusion = CONFIGURATION['signal_diffusion']


def write_signal(i_signal):
    signal_name = CONFIGURATION['signals'].format(sig=i_signal)
    analyzer.save_signal(CONFIGURATION['signal_component_name'], signal_name,
                         x0=mean, mean=mean, correlation_time=corr_time, diffusion=diffusion)


def write_response(i_signal, i_response):
    signal_name = CONFIGURATION['signals'].format(sig=i_signal)
    analyzer.simulate_trajectory('response', CONFIGURATION['responses'].format(sig=i_signal, res=i_response), [signal_name],
                                 seed=i_response)


def delete_response(i_signal, i_response):
    os.remove(CONFIGURATION['responses'].format(sig=i_signal, res=i_response))


def load_signal(i):
    return analyzer.load_trajectory(pathlib.Path(CONFIGURATION['signals'].format(sig=i)))


def load_response(sig, res):
    return analyzer.load_trajectory(pathlib.Path(CONFIGURATION['responses'].format(sig=sig, res=res)))


# this is shared across processes
combined_signal = None


def calculate(s):
    num_signals = CONFIGURATION['num_signals']
    num_responses = CONFIGURATION['num_responses']

    combined_response = None
    for r in range(num_responses):
        write_response(s, r)
        response = load_response(s, r)
        delete_response(s, r)

        if combined_response is None:
            length = len(response['timestamps'])
            combined_response = analyzer.empty_trajectory(
                response['components'].keys(), num_responses, length, events=True)
            combined_response['reactions'] = response['reactions']

        for name, val in response['components'].items():
            combined_response['components'][name][r] = val
        combined_response['timestamps'][r] = response['timestamps']
        combined_response['reaction_events'][r] = response['reaction_events']

    this_signal = {'components': {},
                   'timestamps': combined_signal['timestamps'][s]}
    for name, val in combined_signal['components'].items():
        this_signal['components'][name] = val[s]

    r_given_s = analyzer.likelihoods_given_signal(
        combined_response, this_signal)

    signals_without_this = {'components': {}, 'timestamps': None}
    for name, comp in combined_signal['components'].items():
        signals_without_this['components'][name] = np.delete(comp, s, axis=0)
    signals_without_this['timestamps'] = np.delete(
        combined_signal['timestamps'], s, axis=0)

    mean_product = []
    for r in range(num_responses):
        response = {'components': {},
                    'timestamps': None, 'reaction_events': None}
        for name, comp in combined_response['components'].items():
            response['components'][name] = np.expand_dims(comp[r], axis=0)
        response['timestamps'] = np.expand_dims(
            combined_response['timestamps'][r], axis=0)
        response['reaction_events'] = np.expand_dims(
            combined_response['reaction_events'][r], axis=0)
        response['reactions'] = combined_response['reactions']

        r_given_s_prime = analyzer.likelihoods_given_signal(
            response, signals_without_this)

        cumulative_product = np.cumprod(
            r_given_s_prime / r_given_s[r], axis=-1)

        mean_product.append(np.mean(cumulative_product, axis=-2))

    return -np.log(np.array(mean_product))


def main():
    global combined_signal
    num_signals = CONFIGURATION['num_signals']

    for s in range(num_signals):
        write_signal(s)

    for s in range(num_signals):
        signal = load_signal(s)
        if combined_signal is None:
            length = len(signal['timestamps'])
            names = signal['components'].keys()
            combined_signal = analyzer.empty_trajectory(
                names, num_signals, length)

        for name, values in signal['components'].items():
            combined_signal['components'][name][s] = values
        combined_signal['timestamps'][s] = signal['timestamps']

    pool = multiprocessing.Pool()
    mutual_information = pool.map(calculate, range(num_signals))

    np.savez_compressed(CONFIGURATION['output'], mutual_information)
    print("DONE")


if __name__ == '__main__':
    main()
