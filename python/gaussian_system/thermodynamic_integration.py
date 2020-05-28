import numpy as np
from numpy.random import random_sample, randint

import pandas as pd

from scipy.stats import multivariate_normal
from scipy.special import logsumexp

from .system import System, multivariate_gaussian_logpdf, decompose

from numba import njit, objmode

import time
from tqdm import tqdm


@njit
def potential(conf, e_val_prior, prec_U_prior, e_val_joint, prec_U_joint, theta):
    n_dim = len(conf) // 2
    signal = conf[:n_dim]

    result = 0.0
    if theta > 0.0:
        result += theta * multivariate_gaussian_logpdf(conf, e_val_joint, prec_U_joint)
    if (1.0 - theta) > 0.0:
        result += (1.0 - theta) * multivariate_gaussian_logpdf(
            signal, e_val_prior, prec_U_prior
        )

    return -result


@njit
def propose_conf(previous_conf, scale):
    n_dim = len(previous_conf) // 2
    offset = np.zeros_like(previous_conf)
    # offset[:n_dim] = (random_sample(n_dim) - 0.5) * scale
    offset[randint(0, n_dim)] = (random_sample() - 0.5) * scale
    return previous_conf + offset


@njit
def generate_samples_mcmc(initial_conf, c_z, scale, equilibrate=1000, theta=1.0):
    accepted = 0
    rejected = 0

    n_dim = len(initial_conf) // 2

    e_val_prior, prec_U_prior = decompose(c_z[:n_dim, :n_dim])
    e_val_joint, prec_U_joint = decompose(c_z)

    current_conf = initial_conf
    current_pot = potential(
        current_conf, e_val_prior, prec_U_prior, e_val_joint, prec_U_joint, theta
    )

    while True:
        proposed_conf = propose_conf(current_conf, scale)
        proposed_pot = potential(
            proposed_conf, e_val_prior, prec_U_prior, e_val_joint, prec_U_joint, theta
        )

        if random_sample() < np.exp(current_pot - proposed_pot):
            accepted += 1
            current_conf = proposed_conf
            current_pot = proposed_pot
            if accepted % equilibrate == 0 and rejected > 0:
                acceptance_rate = np.divide(accepted, rejected)
                yield current_conf, acceptance_rate
                accepted = 0
                rejected = 0
        else:
            rejected += 1


def estimate_marginal_density(
    initial_configuration, num_samples, system, t, scale, equilibrate=1000, theta=1.0
):
    current_time = time.perf_counter()

    samples = np.zeros((num_samples, len(initial_configuration)))
    acceptance = np.zeros(num_samples)

    generator = generate_samples_mcmc(
        initial_conf=initial_configuration,
        c_z=system.corr_z(t),
        scale=scale,
        equilibrate=equilibrate,
        theta=theta,
    )

    for sample_slot, rate_slot, (sample, rate) in zip(
        samples, acceptance[:, np.newaxis], generator
    ):
        sample_slot[:] = sample
        rate_slot[:] = rate

    elapsed = time.perf_counter() - current_time

    n_dim = t.shape[0]
    responses = samples[:, n_dim:]
    signals = samples[:, :n_dim]

    ll = system.log_likelihood(responses, signals, t)

    return pd.DataFrame(
        {
            "log_likelihood": ll,
            "acceptance_rates": acceptance,
            "time": elapsed,
            "scale": scale,
            "skip": equilibrate,
            "theta": theta,
        }
    )

    # ll = system.log_likelihood(responses, signals, t)

    # estimate = -logsumexp(-ll, b=1 / len(ll))

    # error = system.log_marginal(responses[0], t)[0] - estimate

    # return {
    #     "num_samples": num_samples,
    #     "log_marginal": estimate,
    #     "error": error,
    #     "acceptance_rate": acceptance.mean(),
    #     "time": elapsed,
    #     "scale": scale,
    #     "skip": equilibrate,
    #     "theta": theta,
    # }
