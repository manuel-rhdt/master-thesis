import pathlib
import subprocess
import sys

import numpy as np
import scipy
import scipy.stats


def ornstein_uhlenbeck_path(x0, t, mean_rev_speed, mean_rev_level, vola):
    """ Simulates a sample path for an Ornstein-Uhlenbeck process."""
    assert len(t) > 1
    x = scipy.stats.norm.rvs(size=len(t))
    x[0] = x0
    dt = np.diff(t)
    scale = std(dt, mean_rev_speed, vola)
    x[1:] = x[1:] * scale
    for i in range(1, len(x)):
        x[i] += mean(x[i - 1], dt[i - 1], mean_rev_speed, mean_rev_level)
    return x


def std(t, mean_rev_speed, vola):
    return np.sqrt(variance(t, mean_rev_speed, vola))


def variance(t, mean_rev_speed, vola):
    assert mean_rev_speed >= 0
    assert vola >= 0
    return vola * vola * (1.0 - np.exp(- 2.0 * mean_rev_speed * t)) / (2 * mean_rev_speed)


def mean(x0, t, mean_rev_speed, mean_rev_level):
    assert mean_rev_speed >= 0
    return x0 * np.exp(-mean_rev_speed * t) + (1.0 - np.exp(- mean_rev_speed * t)) * mean_rev_level


if __name__ == '__main__':
    signal_path = pathlib.Path('signal')
    response_path = pathlib.Path('response')
    signal_path.mkdir(exist_ok=True)
    response_path.mkdir(exist_ok=True)

    executable = pathlib.Path('../cmake-build-debug/Gillespie')

    for trajNum in range(10):
        command = [str(executable), 'signal.inp',
                   '-o', str(signal_path / (str(trajNum) + '.traj.txt')),
                   '-s', str(trajNum)]
        print(' '.join(command), file=sys.stderr)
        subprocess.run(command, stdout=subprocess.PIPE)
