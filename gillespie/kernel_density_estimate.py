"""Implement gaussian kernel density estimation, strongly inspired by code from scipy.

The main difference is that the implementation here is much more limited. However it
uses less memory for some computations which turns out to make a crucial difference
for our purposes.
"""

import numpy
from numba import njit

from . import likelihood


@njit
def silverman_factor(n, dimensions):
    """Compute the Silverman factor.
    Returns
    -------
    s : float
        The silverman factor.
    """
    return numpy.power(n * (dimensions + 2) / 4, -1 / (dimensions + 4))


@njit(cache=True)
def estimate_log_density(points, dataset):
    """Calculates a log probability estimate for every point in `points` where the
    probability estimate is based on gaussian kernel density estimation from `dataset`.

    Arguments:
        points {numpy.ndarray} -- The points for which the density estimate is
                                  evaluated.
        dataset {numpy.ndarray} -- The dataset used for the density estimation

    Returns:
        numpy.ndarray -- An array of log probabilities
    """
    dataset = numpy.atleast_2d(dataset)
    points = numpy.atleast_2d(numpy.asarray(points))

    d, n = dataset.shape
    d_points, m = points.shape

    assert d == d_points, "dimension of dataset must match dimension of points"

    cov_factor = silverman_factor(n, d)
    if d == 1:
        covariance = numpy.atleast_2d(numpy.cov(dataset[0]))
    else:
        covariance = numpy.atleast_2d(numpy.cov(dataset))
    inv_cov_scaled = numpy.linalg.inv(covariance) / cov_factor ** 2

    log_norm_factor = 0.5 * numpy.log(numpy.linalg.det(inv_cov_scaled / (2 * numpy.pi)))

    tmp = numpy.zeros((n, m))
    for i in range(n):
        diff = numpy.empty((d, m))
        for k in range(m):
            for l in range(d):
                diff[l, k] = dataset[l, i] - points[l, k]
        tdiff = numpy.dot(inv_cov_scaled, diff)
        tmp[i] = -numpy.sum(diff * tdiff, axis=0) / 2.0

    return likelihood.logsumexp(tmp) - numpy.log(n) + log_norm_factor
