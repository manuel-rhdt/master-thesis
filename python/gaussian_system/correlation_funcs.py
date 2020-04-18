import numpy as np
import scipy


def time_matrix(N: int, delta_t: float):
    time_stamps = np.expand_dims(np.linspace(0, (N - 1) * delta_t, N), 0)
    return np.atleast_2d(time_stamps - time_stamps.T)


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
        return np.where(t >= 0, self._corr_xs_pos(t), self._corr_sx_pos(-t),)

    def corr_sx(self, t):
        return np.where(t >= 0, self._corr_sx_pos(t), self._corr_xs_pos(-t),)

    def corr_xx(self, t):
        c1 = np.exp(-self.mu * np.abs(t)) - np.exp(-self.lamda * np.abs(t))
        c2 = np.exp(-self.mu * np.abs(t))
        d1 = self.rho ** 2 / (self.lamda ** 2 - self.mu ** 2) * self.kappa / self.lamda
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
        return np.block([[c_ss, c_sx], [c_xs, c_xx]])

    def marginal_entropy(self, t):
        c_xx = self.corr_xx(t)
        return scipy.stats.multivariate_normal(cov=c_xx).entropy()
