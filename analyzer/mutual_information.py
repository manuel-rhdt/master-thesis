#!/usr/bin/env python3
import settings
from analyzer import analyzer
import pathlib
import numpy as np
import os
import multiprocessing
import queue
import math
from tqdm import tqdm


CONFIGURATION = settings.configuration['mutual_information']
mean = CONFIGURATION['signal_mean']
corr_time = CONFIGURATION['signal_correlation_time']
diffusion = CONFIGURATION['signal_diffusion']
resolution = CONFIGURATION['signal_resolution']


def write_signal(i_signal):
    signal_name = CONFIGURATION['signals'].format(sig=i_signal)
    analyzer.save_signal(CONFIGURATION['signal_component_name'], signal_name,
                         x0=mean, mean=mean, step_size=1/resolution, correlation_time=corr_time, diffusion=diffusion, duration=CONFIGURATION['signal_duration'])
    # analyzer.simulate_trajectory(
    #    'signal.inp', signal_name, cwd = '../runs/signal', seed = i_signal)


def write_response(i_signal, i_response):
    signal_name = CONFIGURATION['signals'].format(sig=i_signal)
    analyzer.simulate_trajectory(CONFIGURATION['input_name'], CONFIGURATION['responses'].format(sig=i_signal, res=i_response), [signal_name],
                                 seed=i_response)


def delete_response(i_signal, i_response):
    os.remove(CONFIGURATION['responses'].format(sig=i_signal, res=i_response))


def load_signal(i):
    return analyzer.load_trajectory(pathlib.Path(CONFIGURATION['signals'].format(sig=i)))


def load_response(sig, res):
    return analyzer.load_trajectory(pathlib.Path(CONFIGURATION['responses'].format(sig=sig, res=res)))


def calculate(s,):
    num_signals = CONFIGURATION['num_signals']
    num_responses = CONFIGURATION['num_responses']

    combined_signal = calculate.combined_signal
    combined_response = None
    response_len = None
    for r in range(num_responses):
        write_response(s, r)
        response = load_response(s, r)
        delete_response(s, r)
        calculate.queue.put(3)
        response_len = response['timestamps'].shape[0]

        if combined_response is None:
            combined_response = analyzer.empty_trajectory(
                response['components'].keys(), num_responses, response_len, events=True)
            combined_response['reactions'] = response['reactions']

        for name, val in response['components'].items():
            combined_response['components'][name][r] = val
        combined_response['timestamps'][r] = response['timestamps']
        combined_response['reaction_events'][r] = response['reaction_events']

    this_signal = {'components': {},
                   'timestamps': combined_signal['timestamps'][[s]]}
    for name, val in combined_signal['components'].items():
        this_signal['components'][name] = val[[s]]

    r_given_s = analyzer.log_likelihoods_given_signal(
        combined_response, this_signal)
    calculate.queue.put(1)

    signals_without_this = {'components': {}, 'timestamps': None}
    for name, comp in combined_signal['components'].items():
        signals_without_this['components'][name] = np.delete(comp, s, axis=0)
    signals_without_this['timestamps'] = np.delete(
        combined_signal['timestamps'], s, axis=0)

    # ignore underflow in exponentials small numbers
    np.seterr(all='warn', under='ignore')
    mutual_information = np.zeros((num_responses, 2, response_len))
    for r in range(num_responses):
        response = {'components': {},
                    'timestamps': None, 'reaction_events': None}
        for name, comp in combined_response['components'].items():
            response['components'][name] = comp[[r]]
        response['timestamps'] = combined_response['timestamps'][[r]]
        response['reaction_events'] = np.expand_dims(
            combined_response['reaction_events'][r], axis=0)
        response['reactions'] = combined_response['reactions']

        r_given_s_prime = analyzer.log_likelihoods_given_signal(
            response, signals_without_this)

        cumulative_sum = np.cumsum(r_given_s_prime - r_given_s[0, r], axis=-1)
        information_gain = np.mean(np.exp(cumulative_sum), axis=0)
        mutual_information[r, 0] = combined_response['timestamps'][r]
        mutual_information[r, 1] = -np.log(information_gain)

        calculate.queue.put(num_signals - 1)

    name = CONFIGURATION['output']
    np.savez('{}.{}'.format(name, s), mutual_information)


def generate_signals():
    num_signals = CONFIGURATION['num_signals']

    for s in range(num_signals):
        write_signal(s)

    combined_signal = None
    for s in range(num_signals):
        signal = load_signal(s)
        if combined_signal is None:
            length = len(signal['timestamps'])
            names = ['S']
            combined_signal = analyzer.empty_trajectory(
                names, num_signals, length)

        combined_signal['components']['S'][s] = signal['components']['S']
        combined_signal['timestamps'][s] = signal['timestamps']

    return combined_signal


def test_run():
    calculate.combined_signal = generate_signals()
    calculate.queue = multiprocessing.Queue()
    calculate(0)


def process_init(q, s):
    calculate.queue = q
    calculate.combined_signal = s


def main():
    num_signals = CONFIGURATION['num_signals']
    num_responses = CONFIGURATION['num_responses']

    pathlib.Path(CONFIGURATION['output']).parent.mkdir(exist_ok=True)

    combined_signal = generate_signals()

    q = multiprocessing.Queue()

    pool = multiprocessing.Pool(None, process_init, (q, combined_signal))
    result = pool.map_async(calculate, range(num_signals))

    num_tasks = num_signals * num_responses * (num_signals + 2)
    pbar = tqdm(total=num_tasks)
    done_tasks = 0
    while done_tasks < num_tasks:
        update_val = q.get(True)
        pbar.update(update_val)
        done_tasks += update_val

    mutual_information = result.get()
    print("DONE")


if __name__ == '__main__':
    # test_run()
    main()
