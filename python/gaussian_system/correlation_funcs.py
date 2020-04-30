import numpy as np
import scipy
from numba import guvectorize, float64


def time_matrix(N: int, delta_t: float):
    time_stamps = np.expand_dims(np.linspace(0, (N - 1) * delta_t, N), 0)
    return np.atleast_2d(time_stamps - time_stamps.T)


@guvectorize(
    [(float64[:], float64[:], float64[:, :], float64[:, :], float64[:], float64[:])],
    "(n),(m),(n,m),(n,n),(n)->()",
    nopython=True,
    target="parallel",
    fastmath=True,
)
def log_likelihood_jit(x, s, regression_coef, prec_U, e_val, result=None):
    mean = regression_coef @ s
    maha = np.sum(np.square(prec_U.T @ (x - mean)))
    result[0] = -0.5 * (np.log(2 * np.pi) * len(e_val) +
                        np.sum(np.log(e_val)) + maha)


# from https://blogs.sas.com/content/iml/2012/10/31/compute-the-log-determinant-of-a-matrix.html
def logdet_symmetric(matrix):
    lower_triangular = np.linalg.cholesky(matrix)
    return 2 * np.sum(np.log(np.diag(lower_triangular)))


class System:
    def __init__(self, lamda, kappa, rho, mu):
        self.lamda = lamda
        self.kappa = kappa
        self.rho = rho
        self.mu = mu

    def corr_ss(self, t):
        return self.kappa / self.lamda * np.exp(-np.abs(t) * self.lamda)

    def _corr_xs_pos(self, t):
        return (
            self.rho
            * self.kappa
            / self.lamda
            / (self.lamda + self.mu)
            * np.exp(-self.lamda * t)
        )

    def _corr_sx_pos(self, t):
        a = self.rho * self.kappa / self.lamda / (self.lamda - self.mu)
        b1 = (1 + (self.lamda - self.mu) / (self.lamda + self.mu)) * np.exp(
            -self.mu * t
        )
        b2 = -np.exp(-self.lamda * np.abs(t))
        return a * (b1 + b2)

    def corr_xs(self, t):
        return np.where(t >= 0, self._corr_xs_pos(t), self._corr_sx_pos(-t))

    def corr_sx(self, t):
        return np.where(t >= 0, self._corr_sx_pos(t), self._corr_xs_pos(-t))

    def corr_xx(self, t):
        c1 = np.exp(-self.mu * np.abs(t)) - np.exp(-self.lamda * np.abs(t))
        c2 = np.exp(-self.mu * np.abs(t))
        d1 = self.rho ** 2 / (self.lamda ** 2 - self.mu **
                              2) * self.kappa / self.lamda
        d2 = (
            (1 + self.rho / (self.lamda + self.mu))
            * self.kappa
            / self.lamda
            * self.rho
            / self.mu
        )
        return d1 * c1 + d2 * c2

    def corr_z(self, t):
        c_ss = self.corr_ss(t)
        c_sx = self.corr_sx(t)
        c_xs = self.corr_xs(t)
        c_xx = self.corr_xx(t)
        return np.block([[c_ss, c_xs], [c_sx, c_xx]])

    def log_joint(self, s, x, t):
        c_z = self.corr_z(t)
        joint = scipy.stats.multivariate_normal(cov=c_z)
        s_br, x_br = np.broadcast_arrays(np.atleast_2d(s), np.atleast_2d(x))
        joint_sample = np.concatenate((s_br, x_br), axis=-1)
        return joint.logpdf(joint_sample).reshape(s_br.shape[:-1])

    def log_prior(self, s, t):
        c_ss = self.corr_ss(t)
        s2d = np.atleast_2d(s)
        return scipy.stats.multivariate_normal.logpdf(s2d, cov=c_ss).reshape(
            s2d.shape[:-1]
        )

    def log_marginal(self, x, t):
        c_xx = self.corr_xx(t)
        x2d = np.atleast_2d(x)
        return scipy.stats.multivariate_normal.logpdf(x2d, cov=c_xx).reshape(
            x2d.shape[:-1]
        )

    # log likelihood (evaluation of the log density of a conditional Gaussian distribution)
    # evaluates P(x|s)
    def log_likelihood(self, x, s, t):
        c_ss = self.corr_ss(t)
        c_sx = self.corr_sx(t)
        c_xs = self.corr_xs(t)
        c_xx = self.corr_xx(t)
        regression_coef = c_sx @ np.linalg.inv(c_ss)

        # generate the covariance matrix and compute the square root of its inverse
        p_x_given_s_cov = c_xx - regression_coef @ c_xs
        e_val, e_vec = np.linalg.eigh(p_x_given_s_cov)
        prec_U = np.sqrt(np.reciprocal(e_val)) * e_vec

        return log_likelihood_jit(x, s, regression_coef, prec_U, e_val)

    def log_posterior(self, s, x, t):
        c_ss = self.corr_ss(t)
        c_sx = self.corr_sx(t)
        c_xs = self.corr_xs(t)
        c_xx = self.corr_xx(t)
        regression_coef = c_xs @ np.linalg.inv(c_xx)

        # generate the covariance matrix and compute the square root of its inverse
        p_s_given_x_cov = c_ss - regression_coef @ c_sx
        e_val, e_vec = np.linalg.eigh(p_s_given_x_cov)
        prec_U = np.sqrt(np.reciprocal(e_val)) * e_vec

        return log_likelihood_jit(s, x, regression_coef, prec_U, e_val)

    def mutual_information(self, t):
        c_ss = np.atleast_2d(self.corr_ss(t))
        c_xx = np.atleast_2d(self.corr_xx(t))
        z = np.atleast_2d(self.corr_z(t))

        det_c_ss = logdet_symmetric(c_ss)
        det_c_xx = logdet_symmetric(c_xx)
        det_z = logdet_symmetric(z)

        return 0.5 * (det_c_ss + det_c_xx - det_z)

    def conditional_entropy(self, t):
        c_ss = self.corr_ss(t)
        z = self.corr_z(t)

        det_c_ss = logdet_symmetric(c_ss)
        n = c_ss.shape[0]
        det_z = logdet_symmetric(z)

        return 0.5 * (det_z - det_c_ss + n * np.log(2 * np.pi * np.e))

    def marginal_entropy(self, t):
        c_xx = self.corr_xx(t)
        return scipy.stats.multivariate_normal(cov=c_xx).entropy()
