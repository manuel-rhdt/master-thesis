import numpy as np
from scipy.special import logsumexp
from scipy.stats import multivariate_normal

from .correlation_funcs import System


def estimate_log_marginal(num_x, num_s, system: System, t):
    marg_x = multivariate_normal(cov=system.corr_xx(t))
    marg_s = multivariate_normal(cov=system.corr_ss(t))

    x_samples = marg_x.rvs(num_x).reshape((num_x, 1, -1))
    s_samples = marg_s.rvs((num_x, num_s)).reshape((num_x, num_s, -1))

    log_likelihood = system.log_likelihood(x_samples, s_samples, t)
    return x_samples, logsumexp(log_likelihood, axis=-1) - np.log(num_s)


def estimate_marginal_entropy(num_x, num_s, sys: System, t):
    _, estimated_marginal = estimate_log_marginal(num_x, num_s, sys, t)
    return -estimated_marginal.mean(axis=-1)
