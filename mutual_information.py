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

from gillespie import configuration, likelihood, stochastic_sim
from gillespie.kernel_density_estimate import estimate_log_density

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


def calculate(i, averaging_signals, distribution, log_p0_signal):
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

    # TODO: make this configurable
    result_size = 5000
    traj_lengths = np.geomspace(0.01, 1000, num=result_size, dtype=np.single)

    # generate responses from signals

    # first we sample the initial points from the joined distribution
    initial_comps = distribution.take(i, axis=1, mode="wrap")
    assert initial_comps.shape[0] == 2
    sig = generate_signals_sim(1, length=response_len, initial_values=initial_comps[0])
    responses = generate_responses(
        1,
        sig["timestamps"],
        sig["components"],
        length=response_len,
        initial_values=initial_comps[1],
    )

    log_p_x_zero_this = estimate_log_density(
        np.expand_dims(initial_comps, 1), distribution
    ) - estimate_log_density(initial_comps[0], distribution[0])

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
        log_p_x_zero = estimate_log_density(points, distribution) - log_p0_signal

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


def worker_init(signals, distribution):
    """ The task carried out by each individual process.
    """
    # prevent more threads from being created!
    try:
        import mkl
    except ImportError:
        logging.warn("Not using MKL!")
    else:
        mkl.set_num_threads(1)

    pregenerated_signals = {}
    for key, path in signals.items():
        pregenerated_signals[key] = np.load(path, mmap_mode="r")

    points = pregenerated_signals["components"][:, 0, 0]

    log_p0_signal = estimate_log_density(points, distribution[0])

    global WORKER_VARS
    WORKER_VARS["signal"] = pregenerated_signals
    WORKER_VARS["distribution"] = distribution
    WORKER_VARS["log_p0_signal"] = log_p0_signal


def worker_work(i):
    return calculate(
        i,
        averaging_signals=WORKER_VARS["signal"],
        distribution=WORKER_VARS["distribution"],
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
        with distribution_path.open("wb") as dist_file:
            np.save(dist_file, components)
        logging.info("...Done")

    return components.T


def get_or_generate_signals(distribution):
    signal_path = pathlib.Path(
        configuration.get().get("signal_path", OUT_PATH / "signals.npz")
    )
    if signal_path.exists():
        logging.info(f"Using signals from {signal_path}")
        return signal_path
    logging.info("Generate signals...")
    initial_values = distribution[0].take(np.arange(NUM_SIGNALS), mode="wrap")
    signal_length = configuration.get()["signal"]["length"]
    combined_signal = generate_signals_sim(
        NUM_SIGNALS, length=signal_length, initial_values=initial_values
    )
    with signal_path.open("wb") as signal_file:
        np.savez_compressed(signal_file, **combined_signal)
    del combined_signal
    logging.info("...Done")
    return signal_path


def signal_share_mem(signal_path):
    shared_mem_sig = {}
    with np.load(signal_path, allow_pickle=False) as signal:
        for key, value in signal.items():
            path = f"/dev/shm/signal_{key}.npy"
            with open(path, "xb") as file:
                np.save(file, value)
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
        level=logging.DEBUG, format="%(levelname)s: %(asctime)s %(message)s"
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

    distribution = get_or_generate_distribution()
    signal_path = get_or_generate_signals(distribution)
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
            initargs=(signals_shared, distribution),
        ) as executor:
            last_time = datetime.now()
            for res in executor.map(worker_work, range(num_responses)):
                results.append(res)
                if len(results) % NUM_PROCESSES == 0:
                    i = len(results)

                    time_per_iteration = (datetime.now() - last_time) / NUM_PROCESSES
                    last_time = datetime.now()

                    logging.info(
                        "response "
                        f"{i}/{num_responses} = {i/num_responses*100:.1f} % done. "
                        f"{time_per_iteration.total_seconds():.2f} s/it"
                    )
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
        for path in signals_shared.values():
            os.remove(path)

        runinfo["run"]["completed_responses"] = len(results)
        if results:
            output = {}
            for key in results[0].keys():
                output[key] = np.concatenate([res[key] for res in results])

            np.savez_compressed(OUT_PATH / "mutual_information", **output)
        endrun(runinfo)


if __name__ == "__main__":
    main()
