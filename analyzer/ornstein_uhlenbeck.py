import numpy as np
from numba import jit


@jit(nopython=True)
def generate(times, x0, correlation_time, diffusion_constant, mean=0.0):
    """ Simulates a sample path for an Ornstein-Uhlenbeck process."""
    time_steps = np.diff(times)
    gamma = 1.0/correlation_time
    x = np.zeros(times.shape)
    x[0] = x0
    for i in range(1, len(x)):
        dt = time_steps[i-1]
        dW = np.random.normal(loc=0.0, scale=np.sqrt(dt))
        x[i] = x[i-1] - gamma * (x[i-1] - mean) * dt + \
            np.sqrt(2 * diffusion_constant) * dW
    return x
