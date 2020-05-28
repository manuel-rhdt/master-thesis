import numpy as np
import pandas as pd
from scipy.special import logsumexp
from scipy.stats import multivariate_normal

from .system import System, time_matrix


def estimate_log_marginal_at(x, num_s, system: System, t):
    xdim = x.shape[-1]
    x = x.reshape((-1, 1, xdim))
    num_x = x.shape[0]

    marg_s = multivariate_normal(cov=system.corr_ss(t))
    s_samples = marg_s.rvs((num_x, num_s)).reshape((num_x, num_s, -1))
    return system.log_likelihood(x, s_samples, t)


def estimate_log_marginal(num_x: int, num_s: int, system: System, t):
    marg_x = multivariate_normal(cov=system.corr_xx(t))
    x_samples = marg_x.rvs(num_x).reshape((num_x, 1, -1))
    log_likelihood = estimate_log_marginal_at(x_samples, num_s, system, t)
    return x_samples, logsumexp(log_likelihood, b=1/num_s, axis=-1)


def monte_carlo_sim(dim, delta_t, num_x: int, num_s, system: System):
    _, p_x = estimate_log_marginal(
        num_x, num_s, system, time_matrix(dim, delta_t))
    return pd.DataFrame(
        {
            "dim": dim,
            "delta_t": delta_t,
            "num_responses": 1,
            "num_signals": num_s,
            "marginal_entropy": -p_x,
        }
    )


def estimate_marginal_entropy(num_x, num_s, sys: System, t):
    _, estimated_marginal = estimate_log_marginal(num_x, num_s, sys, t)
    return -estimated_marginal.mean(axis=-1)
