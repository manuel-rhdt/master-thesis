#!/usr/bin/env python3

import multiprocessing
from analyzer import analyzer, stochastic_sim
from analyzer import configuration
from numba import njit
import pathlib
import numpy as np
from tqdm import tqdm
from scipy.stats import gaussian_kde
from datetime import datetime, timezone
import platform
import sys
import toml
import concurrent.futures

OUT_PATH = pathlib.Path(configuration.get()["output"])
num_signals = configuration.get()["num_signals"]
try:
    NUM_PROCESSES = configuration.get()["num_processes"]
except KeyError:
    NUM_PROCESSES = multiprocessing.cpu_count()

SIGNAL_NETWORK, RESPONSE_NETWORK = configuration.read_reactions()


def generate_signals_sim(count, length=100000, initial_values=None):
    timestamps = np.zeros((count, length), dtype=np.single)
    trajectory = np.zeros((count, 1, length), dtype=np.int16)
    reaction_events = np.zeros((count, length - 1), dtype=np.uint8)

    # initial values
    if initial_values is not None:
        trajectory[:, 0, 0] = initial_values
    else:
        raise RuntimeError("Need to specify initial values")

    if count > 1:
        with concurrent.futures.ThreadPoolExecutor(
            max_workers=NUM_PROCESSES
        ) as executor:
            for i in range(count):
                executor.submit(
                    stochastic_sim.simulate_one,
                    timestamps[i],
                    trajectory[i],
                    reaction_events=reaction_events[i],
                    reactions=SIGNAL_NETWORK,
                )
    else:
        stochastic_sim.simulate_one(
            timestamps[0],
            trajectory[0],
            reaction_events=reaction_events[0],
            reactions=SIGNAL_NETWORK,
        )

    return {
        "timestamps": timestamps,
        "components": trajectory,
        "reaction_events": reaction_events,
    }


def generate_responses(
    count, signal_timestamps, signal_comps, length=100000, initial_values=None
):
    timestamps = np.zeros((count, length))
    trajectory = np.zeros((count, 1, length), dtype=np.int16)
    reaction_events = np.zeros((count, length - 1), dtype=np.uint8)

    # initial values
    if initial_values is not None:
        trajectory[:, 0, 0] = initial_values
    else:
        raise RuntimeError("Need to specify initial values")

    stochastic_sim.simulate(
        timestamps,
        trajectory,
        reaction_events=reaction_events,
        reactions=RESPONSE_NETWORK,
        ext_components=signal_comps,
        ext_timestamps=signal_timestamps,
    )

    return {
        "components": trajectory,
        "timestamps": timestamps,
        "reaction_events": reaction_events,
    }


@njit(fastmath=True, cache=True)
def log_evaluate_kde(points, dataset, inv_cov):
    d, n = dataset.shape
    num_r, _, m = points.shape
    assert len(inv_cov.shape) == 2

    log_norm_factor = -0.5 * np.log(np.linalg.det(inv_cov / (2 * np.pi)))

    result = np.empty((num_r, m))
    for j in range(num_r):
        tmp = np.zeros((n, m))
        for i in range(n):
            diff = np.empty((d, m))
            for k in range(m):
                for l in range(d):
                    diff[l, k] = dataset[l, i] - points[j, l, k]
            tdiff = np.dot(inv_cov, diff)
            tmp[i] = -np.sum(diff * tdiff, axis=0) / 2.0

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

    response_len = configuration.get()["response"]["length"]

    joined_distr = kde_estimate["joined"]
    signal_distr = kde_estimate["signal"]

    # generate responses from signals

    # first we sample the initial points from the joined distribution
    initial_comps = joined_distr.resample(num_responses)
    sig = generate_signals_sim(
        num_responses, length=response_len, initial_values=initial_comps[0]
    )
    responses = generate_responses(
        num_responses,
        sig["timestamps"],
        sig["components"],
        length=response_len,
        initial_values=initial_comps[1],
    )

    log_p_x_zero_this = joined_distr.logpdf(initial_comps) - signal_distr.logpdf(
        initial_comps[0]
    )

    points = np.empty((num_responses, 2, num_signals))
    points[:, 0, :] = averaging_signals["components"][np.newaxis, :, 0, 0]
    points[:, 1, :] = responses["components"][:, 0, 0, np.newaxis]

    log_p_x_zero = log_evaluate_kde(points, joined_distr.dataset, joined_distr.inv_cov)
    log_p_x_zero -= log_p0_signal

    # TODO: make this configurable
    result_size = 5000
    traj_lengths = np.geomspace(0.01, 1000, num=result_size, dtype=np.single)
    # we create an array with the following dimensions to hold the results of our
    # calculations
    # dimension1 = 2: mutual_information[0] holds the trajectory length times
    #                 mutual_information[1] holds the mutual information values
    # dimension2: arrays of responses
    # dimension3: arrays of trajectories
    mutual_information = np.empty((num_responses, result_size), dtype=np.single)

    response_components = responses["components"]
    response_timestamps = responses["timestamps"]
    reaction_events = responses["reaction_events"]

    analyzer.log_likelihood(
        traj_lengths,
        sig["components"],
        sig["timestamps"],
        response_components,
        response_timestamps,
        reaction_events,
        RESPONSE_NETWORK,
        out=mutual_information,
    )
    mutual_information += log_p_x_zero_this[:, np.newaxis]

    signal_components = averaging_signals["components"]
    signal_timestamps = averaging_signals["timestamps"]
    mutual_information -= analyzer.log_averaged_likelihood(
        traj_lengths,
        signal_components,
        signal_timestamps,
        response_components,
        response_timestamps,
        reaction_events,
        RESPONSE_NETWORK,
        log_p_x_zero.T,
    )

    return {"trajectory_length": traj_lengths, "mutual_information": mutual_information}
    # np.save(OUT_PATH / "mi.{}".format(i), np.swapaxes(mutual_information, 0, 1))


def kde_estimate_p_0(size, traj_length, signal_init, response_init):
    sig = generate_signals_sim(size, length=traj_length, initial_values=signal_init)

    components = np.empty((size, 2), dtype=np.int16)
    components[:, 0] = sig["components"][:, 0, 0]
    components[:, 1] = response_init

    until = sig["timestamps"][:, -1]

    with concurrent.futures.ThreadPoolExecutor(max_workers=NUM_PROCESSES) as executor:
        for i in range(size):
            executor.submit(
                stochastic_sim.simulate_until_one,
                until[i],
                components[i],
                reactions=RESPONSE_NETWORK,
                ext_timestamps=sig["timestamps"][i],
                ext_components=sig["components"][i],
            )

    np.save(OUT_PATH / "equilibrated", components)

    return {
        "joined": gaussian_kde(components.T),
        "signal": gaussian_kde(components[:, 0]),
    }


def worker(tasks, done_queue, signal_path, kde_estimate):
    signal_distr = kde_estimate["signal"]
    pregenerated_signals = np.load(signal_path, mmap_mode="r")
    points = pregenerated_signals["components"][np.newaxis, np.newaxis, :, 0, 0]
    log_p0_signal = log_evaluate_kde(points, signal_distr.dataset, signal_distr.inv_cov)

    for i in iter(tasks.get, "STOP"):
        result = calculate(
            i,
            num_responses=1,
            averaging_signals=pregenerated_signals,
            kde_estimate=kde_estimate,
            log_p0_signal=log_p0_signal,
        )
        done_queue.put(result)


def main():
    OUT_PATH.mkdir(exist_ok=False)
    conf = configuration.get()

    runinfo = {}
    runinfo["run"] = {"started": datetime.now(timezone.utc)}
    if platform.node():
        runinfo["run"]["node"] = platform.node()

    with (OUT_PATH / "info.toml").open("x") as f:
        f.write(configuration.CONF_STRING)
        f.write("\n\n")
        toml.dump(runinfo, f)

    print("Estimate initial distribution...", file=sys.stderr)
    kde_estimate = kde_estimate_p_0(
        size=conf["kde_estimate"]["size"],
        traj_length=conf["kde_estimate"]["signal"]["length"],
        signal_init=conf["kde_estimate"]["signal"]["initial"],
        response_init=conf["kde_estimate"]["response"]["initial"],
    )
    print("generating signals...", file=sys.stderr)
    initial_values = kde_estimate["signal"].resample(size=num_signals)
    signal_length = configuration.get()["signal"]["length"]
    combined_signal = generate_signals_sim(
        num_signals, length=signal_length, initial_values=initial_values
    )
    np.savez(OUT_PATH / "signals.npz", **combined_signal)
    del combined_signal
    print("Done!", file=sys.stderr)

    num_responses = configuration.get()["num_responses"]
    task_queue = multiprocessing.Queue()
    for task in range(num_responses):
        task_queue.put(task)
    done_queue = multiprocessing.Queue()

    workers = [
        multiprocessing.Process(
            target=worker,
            args=(task_queue, done_queue, OUT_PATH / "signals.npz", kde_estimate),
        )
        for _ in range(NUM_PROCESSES)
    ]

    pbar = tqdm(total=num_responses, smoothing=0.1, desc="simulated responses")
    for process in workers:
        process.start()
        task_queue.put("STOP")

    results = []
    for _ in range(num_responses):
        results.append(done_queue.get())
        pbar.update(1)

    for process in workers:
        process.join()

    mutual_information = np.concatenate([res["mutual_information"] for res in results])
    np.savez(
        OUT_PATH / "mutual_information",
        trajectory_length=results[0]["trajectory_length"],
        mutual_information=mutual_information,
    )

    runinfo["run"]["ended"] = datetime.now(timezone.utc)
    runinfo["run"]["duration"] = str(
        runinfo["run"]["ended"] - runinfo["run"]["started"]
    )
    with (OUT_PATH / "info.toml").open("w") as f:
        f.write(configuration.CONF_STRING)
        f.write("\n\n")
        toml.dump(runinfo, f)


if __name__ == "__main__":
    main()
