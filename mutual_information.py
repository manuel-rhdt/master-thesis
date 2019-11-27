#!/usr/bin/env python3

import argparse
import concurrent.futures
import logging
import multiprocessing
import os
import pathlib
import platform
import sys
from datetime import datetime, timezone

import numpy as np
import toml
from numba import njit
from scipy.stats import gaussian_kde

from gillespie import configuration, likelihood, stochastic_sim

OUT_PATH = pathlib.Path(configuration.get()["output"])
NUM_SIGNALS = configuration.get()["num_signals"]
try:
    NUM_PROCESSES = configuration.get()["num_processes"]
except KeyError:
    NUM_PROCESSES = os.cpu_count()

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

    # use multithreading to make use of multiple cpu cores. This works because the
    # stochastic simulation releases the Python GIL.
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

        result[j] = likelihood.logsumexp(tmp) - np.log(n) - log_norm_factor

    return result


def calculate(i, averaging_signals, kde_estimate, log_p0_signal):
    """ Calculates and stores the mutual information for `num_responses` respones.

    This function does the following:
    1. generate signals
    2. generate responses following the signals
    3. calculate the log likelihoods of the responses
    4. calculate the log marginal probabilities of the responses
    5. use both quantities to estimate the mutual information
    """
    num_signals = len(averaging_signals["timestamps"])

    response_len = configuration.get()["response"]["length"]

    joined_distr = kde_estimate["joined"]
    signal_distr = kde_estimate["signal"]

    # TODO: make this configurable
    result_size = 5000
    traj_lengths = np.geomspace(0.01, 1000, num=result_size, dtype=np.single)

    # generate responses from signals

    # first we sample the initial points from the joined distribution
    initial_comps = joined_distr.resample(1)
    sig = generate_signals_sim(1, length=response_len, initial_values=initial_comps[0])
    responses = generate_responses(
        1,
        sig["timestamps"],
        sig["components"],
        length=response_len,
        initial_values=initial_comps[1],
    )

    log_p_x_zero_this = joined_distr.logpdf(initial_comps) - signal_distr.logpdf(
        initial_comps[0]
    )

    conditional_entropy = -(
        likelihood.log_likelihood(
            traj_lengths,
            sig["components"],
            sig["timestamps"],
            responses["components"],
            responses["timestamps"],
            responses["reaction_events"],
            reactions=RESPONSE_NETWORK,
            dtype=np.single,
        )
        + log_p_x_zero_this[:, np.newaxis]
    )

    if num_signals > 0:
        points = np.empty((2, num_signals))
        points[0, :] = averaging_signals["components"][:, 0, 0]
        points[1, :] = responses["components"][0, 0, 0, np.newaxis]

        # calculate conditional distribution
        log_p_x_zero = joined_distr.logpdf(points) - log_p0_signal

        response_entropy = -likelihood.log_averaged_likelihood(
            traj_lengths,
            averaging_signals["components"],
            averaging_signals["timestamps"],
            responses["components"],
            responses["timestamps"],
            responses["reaction_events"],
            reactions=RESPONSE_NETWORK,
            p_zero=np.expand_dims(log_p_x_zero, 0),
            dtype=np.single,
        )

    if num_signals > 0:
        mutual_information = response_entropy - conditional_entropy
        return {
            "trajectory_length": np.expand_dims(traj_lengths, axis=0),
            "mutual_information": mutual_information,
            "conditional_entropy": conditional_entropy,
            "response_entropy": response_entropy,
        }
    else:
        return {
            "trajectory_length": np.expand_dims(traj_lengths, axis=0),
            "conditional_entropy": conditional_entropy,
        }


def generate_p0_distributed_values(size, traj_length, signal_init, response_init):
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

    return components


WORKER_VARS = {}


def worker_init(signals, kde_estimate):
    """ The task carried out by each individual process.
    """
    # prevent more threads from being created!
    try:
        import mkl
    except ImportError:
        logging.warn("Not using MKL!")
    else:
        mkl.set_num_threads(1)

    signal_distr = kde_estimate["signal"]

    pregenerated_signals = {}
    for key, path in signals.items():
        pregenerated_signals[key] = np.load(path, mmap_mode="r")

    points = pregenerated_signals["components"][:, 0, 0]
    log_p0_signal = signal_distr.logpdf(points)

    global WORKER_VARS
    WORKER_VARS["signal"] = pregenerated_signals
    WORKER_VARS["kde_estimate"] = kde_estimate
    WORKER_VARS["log_p0_signal"] = log_p0_signal


def worker_work(i):
    return calculate(
        i,
        averaging_signals=WORKER_VARS["signal"],
        kde_estimate=WORKER_VARS["kde_estimate"],
        log_p0_signal=WORKER_VARS["log_p0_signal"],
    )


def get_or_generate_distribution():
    conf = configuration.get()
    distribution_path = pathlib.Path(
        conf.get("distribution_path", OUT_PATH / "distribution.npy")
    )

    if distribution_path.exists():
        logging.info(f"Using distribution from {distribution_path}")
        components = np.load(distribution_path)
    else:
        logging.info("Estimate initial distribution...")
        components = generate_p0_distributed_values(
            size=conf["kde_estimate"]["size"],
            traj_length=conf["kde_estimate"]["signal"]["length"],
            signal_init=conf["kde_estimate"]["signal"]["initial"],
            response_init=conf["kde_estimate"]["response"]["initial"],
        )
        np.save(distribution_path, components)
        logging.info("...Done")

    return {
        "joined": gaussian_kde(components.T),
        "signal": gaussian_kde(components[:, 0]),
    }


def get_or_generate_signals(kde_estimate):
    signal_path = pathlib.Path(
        configuration.get().get("signal_path", OUT_PATH / "signals.npz")
    )
    if signal_path.exists():
        logging.info(f"Using signals from {signal_path}")
        return signal_path
    logging.info("Generate signals...")
    initial_values = kde_estimate["signal"].resample(size=NUM_SIGNALS)
    signal_length = configuration.get()["signal"]["length"]
    combined_signal = generate_signals_sim(
        NUM_SIGNALS, length=signal_length, initial_values=initial_values
    )
    np.savez_compressed(signal_path, **combined_signal)
    del combined_signal
    logging.info("...Done")
    return signal_path


def signal_share_mem(signal_path):
    shared_mem_sig = {}
    with np.load(signal_path, allow_pickle=False) as signal:
        for key, value in signal.items():
            path = f"/dev/shm/signal_{key}.npy"
            np.save(path, value)
            shared_mem_sig[key] = path
    return shared_mem_sig


def endrun(runinfo):
    conf = configuration.get()
    runinfo["run"]["ended"] = datetime.now(timezone.utc)
    runinfo["run"]["duration"] = str(
        runinfo["run"]["ended"] - runinfo["run"]["started"]
    )
    with (OUT_PATH / "info.toml").open("w") as f:
        f.write(conf["original"])
        f.write("\n\n")
        toml.dump(runinfo, f)


def main():
    logging.basicConfig(
        level=logging.DEBUG, format="%(levelname)s:%(asctime)s %(message)s"
    )

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--no-responses",
        help="don't calculate the mutual information but do generate signals",
        action="store_true",
    )
    arguments = parser.parse_args()

    OUT_PATH.mkdir(exist_ok=True)
    conf = configuration.get()

    runinfo = {}
    runinfo["run"] = {"started": datetime.now(timezone.utc)}
    runinfo["run"]["invocation"] = sys.argv
    if platform.node():
        runinfo["run"]["node"] = platform.node()

    with (OUT_PATH / "info.toml").open("x") as f:
        f.write(conf["original"])
        f.write("\n\n")
        toml.dump(runinfo, f)

    kde_estimate = get_or_generate_distribution()
    signal_path = get_or_generate_signals(kde_estimate)
    signals_shared = signal_share_mem(signal_path)

    if arguments.no_responses:
        endrun(runinfo)
        return

    num_responses = conf["num_responses"]
    results = []
    try:
        with concurrent.futures.ProcessPoolExecutor(
            max_workers=NUM_PROCESSES,
            mp_context=multiprocessing.get_context("spawn"),
            initializer=worker_init,
            initargs=(signals_shared, kde_estimate),
        ) as executor:
            for res in executor.map(worker_work, range(num_responses)):
                results.append(res)
                if len(results) % NUM_PROCESSES == 0:
                    i = len(results)
                    logging.info(f"{i}/{num_responses} = {i/num_responses*100} % done")
                    timediff = datetime.now(timezone.utc) - runinfo["run"]["started"]
                    with (OUT_PATH / "progress.txt").open(mode="w") as progress_file:
                        print(
                            f"{i} / {num_responses} responses done\n"
                            f"{i/num_responses * 100} % "
                            f"in {str(timediff)}",
                            file=progress_file,
                        )
    except BaseException as error:
        logging.error("Aborting calculation due to exception.")
        runinfo["run"]["error"] = repr(error)
        raise
    finally:
        runinfo["run"]["completed_responses"] = len(results)
        if results:
            output = {}
            for key in results[0].keys():
                output[key] = np.concatenate([res[key] for res in results])

            np.savez(OUT_PATH / "mutual_information", **output)
        endrun(runinfo)


if __name__ == "__main__":
    main()
