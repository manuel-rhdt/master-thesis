#!/usr/bin/env python3

import dask
from dask import delayed
from dask.distributed import Client

import numpy as np

import gillespie.configuration
from gillespie import trajectory_array as ta
from gillespie.simulate import simulate
from gillespie.likelihood import log_p
from dask import array as da

from scipy.special import logsumexp


def conditional_likelihood(traj_lengths, signal, response, network):
    result = []
    for sig, res, chunk in zip(
        signal.to_delayed(), response.to_delayed(), signal.chunks[0]
    ):
        ll = delayed(log_p)(traj_lengths, sig, res, network)
        result.append(
            da.from_delayed(ll, shape=(chunk,) + traj_lengths.shape, dtype=np.double)
        )

    return da.concatenate(result)


def log_p_multi(traj_lengths, signal, response, network):
    result = np.empty((len(response), len(signal), len(traj_lengths)))
    for i, r in enumerate(response):
        result[i] = log_p(traj_lengths, signal, r, network)
    return result


def marginal_likelihood(traj_lengths, signal, response, network):
    level1 = []

    for res, res_chunk in zip(response.to_delayed(), response.chunks[0]):
        level2 = []
        for sig, sig_chunk in zip(signal.to_delayed(), signal.chunks[0]):
            ll = delayed(log_p_multi)(traj_lengths, sig, res, network)
            level2.append(
                da.from_delayed(
                    ll,
                    shape=(res_chunk, sig_chunk) + traj_lengths.shape,
                    dtype=np.double,
                )
            )
        level1.append(da.concatenate(level2, axis=1))

    return da.concatenate(level1)


def signals_and_responses(count, batch, length, sig_network, res_network):
    signals = []
    responses = []
    for i in range(count // batch):
        sig = delayed(simulate)(
            batch, length, sig_network, initial_value=50, pure=False
        )
        res = delayed(simulate)(
            batch, length, res_network, ext_trajectory=sig, initial_value=50, pure=False
        )
        sig = ta.from_delayed(sig, batch, length)
        res = ta.from_delayed(res, batch, length)
        signals.append(sig)
        responses.append(res)

    signals = ta.concatenate(signals)
    responses = ta.concatenate(responses)
    return signals, responses


def main():
    # dask.config.set(scheduler='synchronous')
    _ = Client()
    conf = gillespie.configuration.get("configuration.toml")
    sig_network, res_network = gillespie.configuration.read_reactions(conf)

    length = 1_000
    count_signals = 1_000
    count_res_per_signal = 1000
    count_avrg_signals = 20_000
    batch = 125

    signals = []
    responses = []
    for s in range(count_signals // batch):
        sig = delayed(simulate)(
            batch, length, sig_network, initial_value=50, pure=False
        )
        sig = ta.from_delayed(sig, batch, length)
        response = np.empty(batch, dtype=object)
        for r in range(batch):
            res = delayed(simulate)(
                count_res_per_signal,
                length,
                res_network,
                ext_trajectory=sig[r],
                initial_value=50,
                pure=False,
            )
            response[r] = ta.from_delayed(res, count_res_per_signal, length)
        responses.append(response)
        signals.append(sig)

    signals = ta.concatenate(signals)
    responses = np.concatenate(responses)

    _, responses2 = signals_and_responses(
        30_000, batch, length, sig_network, res_network
    )

    averaging_signals = []
    for i in range(count_avrg_signals // batch):
        sig = delayed(simulate)(
            batch, length, sig_network, initial_value=50, pure=False
        )
        sig = ta.from_delayed(sig, batch, length)
        averaging_signals.append(sig)
    averaging_signals = ta.concatenate(averaging_signals)

    last_ts = da.concatenate(
        [r.timestamps[:, -1] for r in responses] + [responses2.timestamps[:, -1]]
    )
    min_length = da.min(last_ts).compute()

    traj_lengths = np.linspace(0, min_length, 1000)
    ce = []
    for s, r in zip(signals, responses):
        ce.append(-conditional_likelihood(traj_lengths, s, r, res_network))
    ce = da.concatenate(ce)

    me = marginal_likelihood(traj_lengths, averaging_signals, responses2, res_network)
    me = -da.reduction(
        me, chunk=logsumexp, aggregate=logsumexp, axis=1, dtype=np.double
    ) + np.log(me.shape[1])
    mi = ce.mean(axis=0) - me.mean(axis=0)

    ce, me, mi = dask.compute(ce.mean(axis=0), me.mean(axis=0), mi)
    np.savez(
        "result",
        trajectory_length=traj_lengths,
        marginal_entropy=me,
        conditional_entropy=ce,
        mutual_information=mi,
    )


if __name__ == "__main__":
    main()
