import unittest
from gaussian_system import System, time_matrix
import numpy as np
from scipy.stats import multivariate_normal


class TestLikelihood(unittest.TestCase):
    def test_trivial(self):
        self.assertEqual("a", "a")

    def test_log_likelihood(self):
        system = System(0.1, 0.2, 0.3, 0.4)
        x = np.random.random_sample((100, 100))
        s = np.random.random_sample((100, 100))
        t = time_matrix(100, 0.1)
        result = system.log_likelihood(x, s, t)

        c_ss = system.corr_ss(t)
        c_sx = system.corr_sx(t)
        c_xs = system.corr_xs(t)
        c_xx = system.corr_xx(t)
        regression_coef = c_sx @ np.linalg.inv(c_ss)
        p_x_given_s_cov = c_xx - regression_coef @ c_xs

        test = []
        for x, s in zip(x, s):
            likelihood_distr = multivariate_normal(
                cov=p_x_given_s_cov, mean=regression_coef @ s
            )
            likelihood = likelihood_distr.logpdf(x)
            test.append(likelihood)

        for val1, val2 in zip(result, test):
            self.assertAlmostEqual(val1, val2)

    def test_mutual_information(self):
        system = System(0.1, 0.2, 0.3, 0.4)
        t = time_matrix(100, 0.1)
        self.assertAlmostEqual(
            system.mutual_information(t),
            system.marginal_entropy(t) - system.conditional_entropy(t),
        )

    def test_distributions(self):
        system = System(0.1, 0.2, 0.3, 0.4)
        t = time_matrix(100, 0.1)

        sample = np.random.random_sample((100, 200))
        s = sample[:, :100]
        x = sample[:, 100:]

        # three mathematically identical ways to compute the joint logpdf
        # log(P(s,x)) = log(P(x|s)) + log(P(s)) = log(P(s|x)) + log(P(x))
        val1 = system.log_likelihood(x, s, t) + system.log_prior(s, t)
        val2 = system.log_posterior(s, x, t) + system.log_marginal(x, t)
        val3 = system.log_joint(s, x, t)
        for s1, s2, s3 in zip(val1, val2, val3):
            self.assertAlmostEqual(s1, s2)
            self.assertAlmostEqual(s1, s3)


if __name__ == "__main__":
    unittest.main()
