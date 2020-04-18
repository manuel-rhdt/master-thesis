import numpy
from scipy.special import logsumexp
from scipy.stats import multivariate_normal

from .correlation_funcs import System

# log likelihood (evaluation of the log density of a conditional Gaussian distribution)
# evaluates P(x|s)
def log_likelihood(x, s, sys: System, t):
    c_ss = sys.corr_ss(t)
    c_sx = sys.corr_sx(t)
    c_xs = sys.corr_xs(t)
    c_xx = sys.corr_xx(t)
    regression_coef = c_xs @ numpy.linalg.inv(c_ss)
    mean = numpy.inner(s, regression_coef)

    # generate the covariance matrix and compute the square root of its inverse
    p_x_given_s_cov = c_xx - regression_coef @ c_sx
    e_val, e_vec = numpy.linalg.eigh(p_x_given_s_cov)
    prec_U = numpy.sqrt(numpy.reciprocal(e_val)) * e_vec

    maha = numpy.sum(numpy.square(numpy.inner(x - mean, prec_U.T)), axis=-1)
    return -0.5 * (
        numpy.log(2 * numpy.pi) * len(e_val) + numpy.sum(numpy.log(e_val)) + maha
    )


def log_marginal_entropy_power(sys: System, t, num_x: int, num_s: int, power: int):
    marg_x = multivariate_normal(cov=sys.corr_xx(t))
    marg_s = multivariate_normal(cov=sys.corr_ss(t))

    x_samples = marg_x.rvs(num_x).reshape((num_x, 1, 1, -1))

    s_shape = (num_x, num_s, power)
    s_samples = marg_s.rvs(s_shape).reshape(s_shape + (-1,))

    log_l = log_likelihood(x_samples, s_samples, sys, t)
    return logsumexp(numpy.sum(log_l, axis=-1), axis=(0, 1)) - numpy.log(
        numpy.double(num_x * num_s)
    )


log_marginal_entropy_power = numpy.vectorize(
    log_marginal_entropy_power, excluded=["sys", "t", 0, 1]
)
