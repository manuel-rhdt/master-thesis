import scipy
import numpy as np
from numba import jit


@jit(nopython=True)
def iterate_process(x, time_steps, noise_steps, gamma, mean):
    for i in range(1, len(x)):
        x[i] = x[i-1] - gamma * (x[i-1] - mean) * \
            time_steps[i-1] + noise_steps[i-1]


def generate(times, x0, correlation_time, diffusion_constant, mean=0.0):
    """ Simulates a sample path for an Ornstein-Uhlenbeck process."""
    rvs = scipy.stats.norm.rvs(size=len(times) - 1)
    time_steps = np.diff(times)
    noise_steps = (2 * diffusion_constant) * rvs * np.sqrt(time_steps)

    result = np.zeros(len(times))
    result[0] = x0
    iterate_process(result, time_steps, noise_steps,
                    1.0/correlation_time, mean)
    return result
