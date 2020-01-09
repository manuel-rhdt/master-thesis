#!/usr/bin/env python3

import argparse
import base64
import concurrent.futures
import hashlib
import logging
import multiprocessing
import os
import pathlib
import platform
import sys
from contextlib import contextmanager
from datetime import datetime, timezone

import numpy as np
import toml

from gillespie import configuration, likelihood, stochastic_sim
from gillespie.kernel_density_estimate import estimate_log_density

OUT_PATH = pathlib.Path(configuration.get()["output"])
try:
    NUM_PROCESSES = configuration.get()["num_processes"]
except KeyError:
    NUM_PROCESSES = os.cpu_count()

SIGNAL_NETWORK, RESPONSE_NETWORK = configuration.read_reactions()


def generate_signals_sim(count, length=100000, initial_values=None, threads=True):
    timestamps = np.zeros((count, length), dtype=np.double)
    trajectory = np.zeros((count, 1, length), dtype=np.int16)
    reaction_events = np.zeros((count, length - 1), dtype=np.uint8)

    # initial values
    if initial_values is not None:
        trajectory[:, 0, 0] = initial_values
    else:
        raise RuntimeError("Need to specify initial values")

    # use multithreading to make use of multiple cpu cores. This works because the
    # stochastic simulation releases the Python GIL.
    if count > 1 and threads:
        with concurrent.futures.ThreadPoolExecutor(
            max_workers=min(NUM_PROCESSES, count)
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
        stochastic_sim.simulate(
            timestamps,
            trajectory,
            reaction_events=reaction_events,
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

    stochastic_sim.simulate_ext(
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


def generate_response_distribution_from_past_signals(
    num_signals, traj_length, signal_init, response_init
):
    """ Returns sampled response values for signal trajectories ending at `signal_init`.
    """
    # for now we exploit that the signal is time-symmetric
    sig = generate_signals_sim(
        num_signals, length=traj_length, initial_values=signal_init, threads=False
    )

    # reverse the trajectory
    sig["timestamps"] = (
        -sig["timestamps"][..., ::-1] + sig["timestamps"][..., -1, np.newaxis]
    )
    sig["components"] = sig["components"][..., ::-1]

    _, num_signal_components, _ = sig["components"].shape

    # TODO: generate multiple responses per signal
    components = np.empty((num_signals, 2), dtype=np.int16)
    components[:, 0] = sig["components"][:, 0, 0]
    components[:, 1] = response_init

    until = sig["timestamps"][:, -1]

    stochastic_sim.simulate_until(
        until,
        components,
        reactions=RESPONSE_NETWORK,
        ext_timestamps=sig["timestamps"],
        ext_components=sig["components"],
    )
    del sig

    return components[:, num_signal_components:]


def calculate(i, averaging_signals, signal_stationary_distr):
    """ Calculates and stores the mutual information for `num_responses` respones.

    This function does the following:
    1. generate signals
    2. generate responses following the signals
    3. calculate the log likelihoods of the responses
    4. calculate the log marginal probabilities of the responses
    5. use both quantities to estimate the mutual information
    """
    conf = configuration.get()

    # TODO: make this configurable
    result_size = 5000
    traj_lengths = np.geomspace(0.01, 1000, num=result_size, dtype=np.double)

    batch_size = conf["response"].get("batch_size", 1)
    response_len = configuration.get()["response"]["length"]

    # first we sample the initial points from the joined distribution
    sig0 = signal_stationary_distr.take(i, axis=0, mode="wrap")[0]
    sig = generate_signals_sim(1, length=response_len, initial_values=sig0)

    response_distribution = generate_response_distribution_from_past_signals(
        max(1000, batch_size * 2),
        1000,
        sig["components"][0, 0, 0],
        conf["kde_estimate"]["response"]["initial"],
    )

    responses = generate_responses(
        batch_size,
        sig["timestamps"],
        sig["components"],
        length=response_len,
        initial_values=response_distribution[:batch_size, 0],
    )

    log_p_x_zero_this = estimate_log_density(
        response_distribution[:batch_size, 0], response_distribution[:, 0]
    )

    conditional_entropies = -(
        likelihood.log_likelihood(
            traj_lengths,
            sig["components"],
            sig["timestamps"],
            responses["components"],
            responses["timestamps"],
            responses["reaction_events"],
            reactions=RESPONSE_NETWORK,
            dtype=np.double,
        )
        + log_p_x_zero_this[:, np.newaxis]
    )
    conditional_entropy = np.mean(conditional_entropies, axis=0, keepdims=True)
    del conditional_entropies  # free memory asap

    num_signals = conf["num_signals"]
    if num_signals > 0:
        points = responses["components"][:, 0, 0]
        log_p_x_zero = np.zeros((batch_size, num_signals))
        conditional_distribution = averaging_signals["conditional_distribution"]

        for s in range(num_signals):
            # calculate conditional distribution
            log_p_x_zero[:, s] = estimate_log_density(
                points, conditional_distribution[s, :, 0]
            )

        marginal_entropies = -likelihood.log_averaged_likelihood(
            traj_lengths,
            averaging_signals["components"],
            averaging_signals["timestamps"],
            responses["components"],
            responses["timestamps"],
            responses["reaction_events"],
            reactions=RESPONSE_NETWORK,
            p_zero=log_p_x_zero,
            dtype=np.double,
        )
        marginal_entropy = np.mean(marginal_entropies, axis=0, keepdims=True)
        del marginal_entropies

    if num_signals > 0:
        return {
            "trajectory_length": np.expand_dims(traj_lengths, axis=0),
            "conditional_entropy": conditional_entropy,
            "response_entropy": marginal_entropy,
        }
    else:
        return {
            "trajectory_length": np.expand_dims(traj_lengths, axis=0),
            "conditional_entropy": conditional_entropy,
        }


def generate_signal_stationary_distribution(num_signals, traj_length, signal_init):
    sig = generate_signals_sim(
        num_signals, length=traj_length, initial_values=signal_init
    )

    return sig["components"][..., -1]


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

    global WORKER_VARS
    WORKER_VARS["signal"] = pregenerated_signals
    WORKER_VARS["distribution"] = distribution


def worker_work(i):
    return calculate(
        i,
        averaging_signals=WORKER_VARS["signal"],
        signal_stationary_distr=WORKER_VARS["distribution"],
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
        components = generate_signal_stationary_distribution(
            num_signals=conf["kde_estimate"]["size"],
            traj_length=conf["kde_estimate"]["signal"]["length"],
            signal_init=conf["kde_estimate"]["signal"]["initial"],
        )
        with distribution_path.open("wb") as dist_file:
            np.save(dist_file, components)
            logging.info(f"Saving distribution in {distribution_path}")
        logging.info("...Done")

    return components.T


def get_or_generate_signals(distribution):
    conf = configuration.get()
    num_signals = conf["num_signals"]
    signal_path = pathlib.Path(conf.get("signal_path", OUT_PATH / "signals.npz"))
    if signal_path.exists():
        logging.info(f"Using signals from {signal_path}")
        return signal_path
    logging.info("Generate signals...")
    initial_values = distribution[0].take(np.arange(num_signals), mode="wrap")
    signal_length = configuration.get()["signal"]["length"]
    signals = generate_signals_sim(
        num_signals, length=signal_length, initial_values=initial_values
    )
    _, num_signal_components, _ = signals["components"].shape
    logging.info(f"Generated {num_signals} signals.")

    num_past_trajectories = 250
    conditional_distribution = np.zeros(
        (num_signals, num_past_trajectories, num_signal_components)
    )

    with concurrent.futures.ThreadPoolExecutor(max_workers=NUM_PROCESSES) as executor:

        def gen(args):
            i, s0 = args
            result = generate_response_distribution_from_past_signals(
                num_past_trajectories,
                1000,
                s0,
                conf["kde_estimate"]["response"]["initial"],
            )
            conditional_distribution[i] = result

        executor.map(gen, enumerate(signals["components"][:, 0, 0]))

    signals["conditional_distribution"] = conditional_distribution
    logging.info(f"Generated conditional distribution of response.")

    with signal_path.open("wb") as signal_file:
        np.savez_compressed(signal_file, **signals)
        logging.info(f"Written signals to {signal_path}")
    del signals
    return signal_path


@contextmanager
def signal_share_mem(signal_path):
    shared_mem_sig = {}
    conf = configuration.get()
    hasher = hashlib.sha256()
    hasher.update(str(conf).encode("utf-8"))
    digest = hasher.digest()
    digest = base64.urlsafe_b64encode(digest)[:8]

    with np.load(signal_path, allow_pickle=False) as signal:
        assert (
            signal["timestamps"].shape[0] >= conf["num_signals"]
        ), f"Not enough signals are stored in file {signal_path}"
        for key, value in signal.items():
            path = f"/dev/shm/signal_{key}_{digest.decode('ascii')}.npy"
            with open(path, "xb") as file:
                np.save(file, value)
                logging.info(f"Created file in shared memory: {path}")
            shared_mem_sig[key] = path

    try:
        yield shared_mem_sig
    finally:
        # clean up the created files after they're no longer needed
        for path in shared_mem_sig.values():
            try:
                os.remove(path)
                logging.info(f"Removed {path}")
            except OSError as e:
                logging.error(f"Could not remove {path}: {e}")


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


def perform_simulation(sim, conf):
    with concurrent.futures.ProcessPoolExecutor(
        max_workers=NUM_PROCESSES,
        mp_context=multiprocessing.get_context("spawn"),
        initializer=worker_init,
        initargs=(signals_shared, distribution),
    ) as executor:
        last_time = datetime.now()
        futures = []
        for r in range(num_responses):
            futures.append(executor.submit(worker_work, r))

        for fut in futures:
            results.append(fut.result())
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

    if arguments.no_responses:
        endrun(runinfo)
        return

    for sim in conf.get("simulation", []):
        perform_simulation(sim)

    with signal_share_mem(signal_path) as signals_shared:
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
                futures = []
                for r in range(num_responses):
                    futures.append(executor.submit(worker_work, r))

                for fut in futures:
                    results.append(fut.result())
                    if len(results) % NUM_PROCESSES == 0:
                        i = len(results)

                        time_per_iteration = (
                            datetime.now() - last_time
                        ) / NUM_PROCESSES
                        last_time = datetime.now()

                        logging.info(
                            "response "
                            f"{i}/{num_responses} = {i/num_responses*100:.1f} % done. "
                            f"{time_per_iteration.total_seconds():.2f} s/it"
                        )
                        timediff = (
                            datetime.now(timezone.utc) - runinfo["run"]["started"]
                        )
                        with (OUT_PATH / "progress.txt").open(
                            mode="w"
                        ) as progress_file:
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

                np.savez_compressed(OUT_PATH / "mutual_information", **output)
                logging.info(f"Saved results to {OUT_PATH / 'mutual_information'}")
            endrun(runinfo)


if __name__ == "__main__":
    main()
