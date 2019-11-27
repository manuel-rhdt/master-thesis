import unittest

import numpy
import scipy.stats
from numpy.random import default_rng

from gillespie import kernel_density_estimate


class TestSim(unittest.TestCase):
    def test_kde(self):
        size = 100
        dataset = numpy.transpose(
            default_rng().multivariate_normal(
                mean=[100.0, 200.0], cov=[[10, 5], [5, 100]], size=size
            )
        )

        points = numpy.array([numpy.linspace(0, 400, 100), numpy.linspace(0, 400, 100)])

        scipy_estimate = scipy.stats.gaussian_kde(dataset, bw_method="silverman")
        our_estimate = kernel_density_estimate.estimate_log_density(points, dataset)

        for i, (a, b) in enumerate(zip(scipy_estimate.logpdf(points), our_estimate)):
            self.assertAlmostEqual(a, b, msg=f"value {i} did not match!")


if __name__ == "__main__":
    unittest.main()
