#!/usr/bin/env python3

import pathlib
import sys
from datetime import datetime, timezone

import numpy as np
import toml
from dask import array as da
from dask import delayed
from dask.distributed import Client, LocalCluster
from scipy.special import logsumexp

import gillespie.configuration
from gillespie import trajectory_array as ta
from gillespie.likelihood import log_p, log_p_multi
from gillespie.simulate import simulate


class Run:
    def __init__(self, config):
        self.conf = config
        self.out_path = pathlib.Path(self.conf["output"])
        self.out_path.mkdir(exist_ok=True)

        self.runinfo = {}
        self.runinfo["run"] = {"started": datetime.now(timezone.utc)}
        self.runinfo["run"]["invocation"] = sys.argv

        with (self.out_path / "info.toml").open("x") as f:
            f.write(self.conf["original"])
            f.write("\n\n")
            toml.dump(self.runinfo, f)

    def __del__(self):
        self.runinfo["run"]["ended"] = datetime.now(timezone.utc)
        self.runinfo["run"]["duration"] = str(
            self.runinfo["run"]["ended"] - self.runinfo["run"]["started"]
        )
        with (self.out_path / "info.toml").open("w") as f:
            f.write(self.conf["original"])
            f.write("\n\n")
            toml.dump(self.runinfo, f)

    def start(self):
        try:
            computation_start(self.conf, self.out_path / "result.npz")
        except BaseException as error:
            self.runinfo["run"]["error"] = repr(error)
            raise


def conditional_likelihood(traj_lengths, signal, response, network):
    # result = []
    sig = signal.to_delayed()
    # assert len(sig) == 1
    # for res, chunk in zip(response.to_delayed(), response.chunks[0]):
    #     ll = delayed(log_p)(traj_lengths, sig[0], res, network)
    #     result.append(
    #         da.from_delayed(ll, shape=(chunk,) + traj_lengths.shape, dtype=np.double)
    #     )

    # return da.concatenate(result, axis=0)
    ll = delayed(log_p)(traj_lengths, sig[0], response.to_delayed()[0], network)
    return da.from_delayed(ll, shape=(response.chunks[0][0], response.chunks[1][0]), dtype=np.double)


def marginal_likelihood(traj_lengths, signal, response, network):
    level1 = []

    for res, res_chunk in zip(response.to_delayed(), response.chunks[0]):
        level2 = []
        for sig, sig_chunk in zip(signal.to_delayed(), signal.chunks[0]):
            ll = delayed(log_p_multi)(traj_lengths, sig, res, network)
            level2.append(
                da.from_delayed(
                    ll,
                    shape=(res_chunk, sig_chunk, response.chunks[1][0]),
                    dtype=np.double,
                )
            )
        level1.append(da.concatenate(level2, axis=1))

    return da.concatenate(level1, axis=0)


def simulate_batched(count, batch, length, network, ext_trajectory=None):
    trajectories = []
    for i in range(count // batch):
        traj = delayed(simulate)(
            batch,
            length,
            network,
            ext_trajectory=ext_trajectory,
            initial_value=50,
            pure=False,
        )
        traj = ta.from_delayed(traj, batch, length)
        trajectories.append(traj)
    return ta.concatenate(trajectories)


def simulate_outer(count, length, network, ext_trajectory):
    r = np.empty(len(ext_trajectory), dtype=object)
    for i, s in enumerate(ext_trajectory):
        response = delayed(simulate)(
            count, length, network, ext_trajectory=s, initial_value=50
        )
        r[i] = ta.from_delayed(response, count, length)
    return r


def signals_and_responses(count, batch, length, sig_network, res_network):
    responses = []
    for i in range(count // batch):
        sig = delayed(simulate)(
            batch, length, sig_network, initial_value=50, pure=False
        )
        res = delayed(simulate)(
            batch, length, res_network, ext_trajectory=sig, initial_value=50, pure=False
        )
        res = ta.from_delayed(res, batch, length)
        responses.append(res)

    responses = ta.concatenate(responses)
    return responses


def computation_start(conf, path):
    cluster = conf.get("scheduler_address", None)
    if cluster is None:
        cluster = LocalCluster(processes=False)
        print("Connecting to", cluster.dashboard_link)
    client = Client(cluster)
    print(client)
    # import dask
    # dask.config.set(scheduler='threads')
    sig_network, res_network = gillespie.configuration.read_reactions(conf)

    length = conf["length"]
    ce_num_signals = conf["conditional_entropy"]["num_signals"]
    ce_res_per_signal = conf["conditional_entropy"]["responses_per_signal"]
    me_num_signals = conf["marginal_entropy"]["num_signals"]
    me_num_responses = conf["marginal_entropy"]["num_responses"]
    batch = conf["batch_size"]

    ce_signals = simulate_batched(ce_num_signals, batch, length, sig_network)
    ce_responses = simulate_outer(ce_res_per_signal, length, res_network, ce_signals)

    me_responses = signals_and_responses(
        me_num_responses, batch, length, sig_network, res_network
    )
    me_signals = simulate_batched(me_num_signals, batch, length, sig_network)

    min_length = da.min(me_responses.timestamps[:, -1]).compute()
    print(min_length)
    traj_lengths = np.linspace(0, min_length, conf["num_points"])
    ce = []
    for s, r in zip(ce_signals, ce_responses):
        ce.append(-conditional_likelihood(traj_lengths, s, r, res_network))
    ce = da.stack(ce).mean(axis=1).compute()

    me = marginal_likelihood(traj_lengths, me_signals, me_responses, res_network)
    me = -da.reduction(
        me, chunk=logsumexp, aggregate=logsumexp, axis=1, dtype=np.double
    ) + np.log(me.shape[1])
    me = me.compute()

    mi = ce.mean(axis=0) - me.mean(axis=0)
    np.savez(
        path,
        trajectory_length=traj_lengths,
        marginal_entropy=me,
        conditional_entropy=ce,
        mutual_information=mi,
    )


def main():
    run = Run(config=gillespie.configuration.get("configuration.toml"))
    run.start()


if __name__ == "__main__":
    main()
