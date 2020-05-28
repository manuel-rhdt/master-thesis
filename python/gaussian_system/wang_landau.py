from .system import System, decompose, multivariate_gaussian_logpdf
from scipy.stats import multivariate_normal
from numpy.random import random_sample
import numpy as np
from numba import njit, objmode

from matplotlib import pyplot as plt


@njit
def is_flat(histogram, flatness=0.95):
    return np.min(histogram) / np.max(histogram) >= flatness


@njit
def flatness(histogram):
    return np.min(histogram) / np.max(histogram)


def energy_histogram(response, system, t):
    marg_s = multivariate_normal(cov=system.corr_ss(t))

    num_signals = 500
    signals = marg_s.rvs(num_signals).reshape((num_signals, -1))

    return system.energy(signals, response, t)


@njit
def wang_landau_jit(
    response, initial_signal, corr_z, epsilon, energy_bins, scale, min_flatness
):

    f_param = 1.0
    response = response.reshape((-1,))
    signal = initial_signal.reshape((-1,))

    configuration = np.zeros(len(response) + len(signal))
    configuration[: len(signal)] = signal
    configuration[len(signal) :] = response

    entropy = np.zeros(len(energy_bins) - 1)
    histogram = np.zeros(len(energy_bins) - 1, dtype=np.uint64)

    e_val, prec_U = decompose(corr_z)

    current_energy_val = -multivariate_gaussian_logpdf(configuration, e_val, prec_U)
    current_energy = np.searchsorted(energy_bins, current_energy_val) - 1
    if current_energy < 0 or current_energy >= len(entropy):
        with objmode():
            print(
                "current energy {} in illegal bin {}".format(
                    current_energy_val, current_energy
                )
            )
        assert current_energy < 0 or current_energy >= len(entropy)

    current_conf = np.copy(configuration)
    offset = np.zeros(len(signal))

    accepted = 0
    rejected = 0

    while f_param > epsilon:
        offset[:] = (random_sample(offset.shape) - 0.5) * scale
        configuration[: len(signal)] = current_conf[: len(signal)] + offset
        proposed_energy_val = -multivariate_gaussian_logpdf(
            configuration, e_val, prec_U
        )
        proposed_energy = np.searchsorted(energy_bins, proposed_energy_val) - 1

        if (
            proposed_energy >= 0
            and proposed_energy < len(entropy)
            and random_sample()
            < np.exp(entropy[current_energy] - entropy[proposed_energy])
        ):
            # If accepted, update the energy and the system:
            current_energy = proposed_energy
            current_conf[:] = configuration
            accepted += 1
        else:
            # If rejected
            rejected += 1

        if rejected > 50000:
            yield histogram, entropy
            rejected = 0
            accepted = 0

        histogram[current_energy] += 1
        entropy[current_energy] += f_param

        if is_flat(histogram, min_flatness):
            with objmode():
                print("update f: {} -> {}".format(f_param, f_param * 0.5))
            histogram[:] = 0
            f_param *= 0.5  # Refine the f parameter
    yield histogram, entropy


def wang_landau(
    response, initial_signal, system, t, epsilon, energy_bins, scale, flatity
):
    marg_s = multivariate_normal(cov=system.corr_ss(t))

    for hist, entropy in wang_landau_jit(
        response, initial_signal, system.corr_z(t), epsilon, energy_bins, scale, flatity
    ):
        pass
        # plt.plot((energy_bins[:-1] + energy_bins[1:]) / 2, hist)
        # plt.show()

    return

    f = 1.0

    entropy = np.zeros(len(energy_bins) - 1)
    histogram = np.zeros(len(energy_bins) - 1)

    current_energy = system.energy(initial_signal, response, t)
    current_energy = np.searchsorted(energy_bins, current_energy) - 1

    accepted = 0
    rejected = 0

    current_conf = np.copy(initial_signal.reshape(-1))
    configuration = np.zeros_like(current_conf)
    offset = np.zeros_like(current_conf)
    while f > epsilon:
        offset[:] = (random_sample(configuration.shape) - 0.5) * scale
        configuration[:] = current_conf + offset
        proposed_energy = system.energy(configuration, response, t)
        proposed_energy = np.searchsorted(energy_bins, proposed_energy) - 1

        if (
            proposed_energy >= 0
            and proposed_energy < len(entropy)
            and random_sample()
            < np.exp(entropy[current_energy] - entropy[proposed_energy])
        ):
            # If accepted, update the energy and the system:
            current_energy = proposed_energy
            current_conf[:] = configuration
            accepted += 1
        else:
            # If rejected
            rejected += 1

        if rejected > 1000:
            print("ratio accepted / rejected " + str(accepted / rejected))
            # print(
            #     f"ratio accepted / rejected {accepted / (rejected+accepted)}, flatness={flatness(histogram)}"
            # )
            rejected = 0
            accepted = 0

        histogram[current_energy] += 1
        entropy[current_energy] += f

        if is_flat(histogram):
            print(f"update f: {f} -> {f * 0.5}")
            histogram[:] = 0.0
            f *= 0.5  # Refine the f parameter
    return entropy
