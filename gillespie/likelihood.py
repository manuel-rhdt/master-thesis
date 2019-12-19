import numpy as np
from numba import generated_jit, jit, types
from numba.typed import List as TypedList
from . import stochastic_sim


# inspired by scipy.special.logsumexp but only supports two-dimensional arrays and
# axis=0
@jit(nopython=True)
def logsumexp(x):
    """ Evaluate log-sum-exp on the outer axis of a 2d array.

    If `z = logsumexp(x)` then `z == log(sum(exp(x), axis = 0))` is approximately
    true.
    """
    x = np.atleast_2d(np.asarray(x))
    _, length = x.shape
    xmax = np.empty(length, dtype=x.dtype)
    for i in range(length):
        xmax[i] = max(x[:, i])

    tmp = np.exp(x - xmax)

    s = np.sum(tmp, axis=0)
    out = np.log(s) + xmax

    return out


@jit(nopython=True)
def calculate_sum_of_reaction_propensities(components, reactions):
    result = np.empty_like(components[0], dtype=np.double)

    for n_reaction in range(reactions.size):
        tmp = np.full_like(result, reactions.k[n_reaction])

        for j_reactant in reactions.reactants[n_reaction]:
            if j_reactant >= 0:
                tmp *= components[j_reactant]

        if n_reaction == 0:
            result[:] = tmp
        else:
            result += tmp

    return result


@jit(nopython=True)
def calculate_selected_reaction_propensities(components, reaction_events, reactions):
    """ Accelerated reaction propensity calculation

    Arguments:
        components {numba.typed.List} -- a list of components
        reaction_k {np.ndarray} -- [description]
        reaction_reactants {np.ndarray} -- a two-dimensional array (num_reactions,
        num_reactants)

    Returns:
        [type] -- [description]
    """
    (length,) = components[0].shape
    assert length == reaction_events.shape[-1]

    propensities = np.empty(reactions.k.shape + components[0].shape, dtype=np.double)
    for n_reaction in range(reactions.size):
        propensities[n_reaction] = reactions.k[n_reaction]

        for j_reactant in reactions.reactants[n_reaction]:
            if j_reactant >= 0:
                propensities[n_reaction] *= components[j_reactant]

    result = np.empty_like(propensities[0])
    for i in range(length):
        event = reaction_events[i]
        result[i] = propensities[event, i]

    return result


@jit(nopython=True)
def time_average(
    trajectory,
    old_timestamps,
    new_timestamps,
    dtype=np.double,
    out=None,
    evaluated=None,
):
    """
    Average of `trajectory` with `old_timestamp` in the time-intervals specified by
    `new_timestamps`.

    Note: This function assumes that both `old_timestamps` and `new_timestamps` are
    ordered.

    # Discussion

    This function is used for the calculation of the mean propensity over a timestep in
    the response trajectory. Since the propensity for a reaction is variable, the
    survival probability depends on the time-integrated value of the total propensity.
    Since the signal trajectory is assumed to be piecewise constant we can calculate the
    average signal between every pair of response timestamps. If we use this average
    signal to compute the total propensity, we don't require the calculation of the time
    integral anymore.

    This function computes the average of a trajectory given by the values `trajectory`
    at the timestamps `old_timestamps` where the averaging happens for the intervals
    between pairs of timestamps in `new_timestamps`.

    Returns a list of averages of size `len(new_timestamps) - 1`.

              |                                    |
              |        +---------------------------| <-- trajectory[k + 1]
              |========|===========================| <== average[i]
              |        |                           |
           ---|--------+                           | <-- trajectory[k]
              | new_timestamps[k]                  |
              |                                    |
              +------------------------------------+---> time
        old_timestamps[i]                  old_timestamps[i+1]
    """
    if out is None:
        out = np.empty(len(new_timestamps) - 1, dtype=dtype)

    def iter_trajectory(old_idx):
        old_idx += 1
        if old_idx > len(old_timestamps) - 1:
            return np.inf, trajectory[-1]
        else:
            return old_timestamps[old_idx], trajectory[old_idx - 1]

    old_idx = 0
    next_trajectory_change, trajectory_value = iter_trajectory(old_idx)
    old_idx += 1
    for idx, (low, high) in enumerate(zip(new_timestamps[:-1], new_timestamps[1:])):
        delta_t = high - low
        assert delta_t > 0.0
        acc = 0.0
        while low < high:
            while next_trajectory_change <= low:
                next_trajectory_change, trajectory_value = iter_trajectory(old_idx)
                old_idx += 1

            acc += trajectory_value * (min(high, next_trajectory_change) - low)
            low = next_trajectory_change

        out[idx] = acc / delta_t
        if evaluated is not None:
            evaluated[idx] = trajectory_value
    return out


@jit(nopython=True, cache=True, nogil=True)
def log_likelihood_inner(
    signal_components,
    signal_timestamps,
    response_components,
    response_timestamps,
    reaction_events,
    reactions,
    dtype=np.dtype(np.double),
    out=None,
):
    """[summary]

    Returns:
        An array of length `(len(response_timestamps) - 1)` that contains the individual
        components of the log-likelihood.
    """
    num_signal_comps, _ = signal_components.shape
    num_response_comps, length = response_components.shape

    if out is not None:
        dtype = out.dtype

    # time-averaged signal components
    rsc = np.empty((num_signal_comps, 2, length - 1), dtype=dtype)
    for i in range(num_signal_comps):
        time_average(
            signal_components[i],
            signal_timestamps,
            response_timestamps,
            dtype=dtype,
            out=rsc[i, 0],
            evaluated=rsc[i, 1],
        )

    components = TypedList()
    for i in range(num_signal_comps):
        components.append(rsc[i, 0])

    for i in range(num_response_comps):
        components.append(response_components[i, 1:].astype(dtype))

    averaged_rates = calculate_sum_of_reaction_propensities(components, reactions)

    for i in range(num_signal_comps):
        components[i] = rsc[i, 1]

    instantaneous_rates = calculate_selected_reaction_propensities(
        components, reaction_events, reactions
    )

    # return the logarithm of `np.cumprod(instantaneous_rates * np.exp(-averaged_rates \
    # * dt))`

    dt = response_timestamps[1:] - response_timestamps[:-1]
    likelihoods = np.log(instantaneous_rates)
    likelihoods -= averaged_rates * dt

    result = out if out is not None else np.empty(length - 1, dtype=dtype)
    accumulator = 0.0
    for i in range(length - 1):
        accumulator += likelihoods[i]
        result[i] = accumulator

    return result


@generated_jit(nopython=True)
def expand_3d(array):
    if isinstance(array, types.Array):
        if array.ndim == 1:
            return lambda array: array.reshape((1, 1, -1))
        elif array.ndim == 2:
            return lambda array: np.expand_dims(array, 0)
        else:
            return lambda array: array


@jit(nopython=True, cache=True, nogil=True)
def log_likelihood(
    traj_lengths,
    signal_components,
    signal_timestamps,
    response_components,
    response_timestamps,
    reaction_events,
    reactions,
    dtype=None,
    out=None,
):
    response_components = expand_3d(response_components)
    signal_components = expand_3d(signal_components)
    response_timestamps = np.atleast_2d(response_timestamps)
    signal_timestamps = np.atleast_2d(signal_timestamps)
    reaction_events = np.atleast_2d(reaction_events)

    num_r, _, _ = response_components.shape
    num_s, _, _ = signal_components.shape
    (length,) = traj_lengths.shape

    assert num_r == num_s or num_s == 1 or num_r == 1

    result = out if out is not None else np.zeros((max(num_r, num_s), length), dtype=dtype)

    for i in range(result.shape[0]):
        r_index = i % num_r
        s_index = i % num_s
        rc = response_components[r_index]
        rt = response_timestamps[r_index]
        sc = signal_components[s_index]
        st = signal_timestamps[s_index]

        log_p = log_likelihood_inner(
            sc, st, rc, rt, reaction_events[r_index], reactions, dtype=dtype
        )

        indices = np.digitize(traj_lengths, rt)
        for j, index in enumerate(indices):
            if index >= len(log_p):
                # the index is out of bounds... just extrapolate for now
                result[i, j] = result[i, j - 1] + np.mean(log_p)
            else:
                result[i, j] = log_p[index]

    return result


@jit(nopython=True, cache=True)
def log_averaged_likelihood(
    traj_lengths,
    signal_components,
    signal_timestamps,
    response_components,
    response_timestamps,
    reaction_events,
    reactions,
    p_zero,
    dtype=None,
    out=None,
):
    """
    Calculates the log likelihoods of the responses for various signals and averages
    over the signals.
    """
    response_components = expand_3d(response_components)
    signal_components = expand_3d(signal_components)
    response_timestamps = np.atleast_2d(response_timestamps)
    signal_timestamps = np.atleast_2d(signal_timestamps)
    reaction_events = np.atleast_2d(reaction_events)

    num_r, _, _ = response_components.shape
    num_s, _, _ = signal_components.shape
    (length,) = traj_lengths.shape

    assert num_r == p_zero.shape[0]
    assert num_s == p_zero.shape[1]

    if dtype is None:
        if out is not None:
            dtype = out.dtype

    result = out if out is not None else np.empty((num_r, length), dtype=dtype)
    for r in range(num_r):
        rc = response_components[r]
        rt = response_timestamps[r]
        tmp = np.empty((num_s, length), dtype=dtype)
        for s in range(num_s):
            tmp[s] = p_zero[r, s]

        indices = np.digitize(traj_lengths, rt)
        for s in range(num_s):
            sc = signal_components[s]
            st = signal_timestamps[s]

            log_p = log_likelihood_inner(
                sc, st, rc, rt, reaction_events[r], reactions, dtype=dtype
            )
            for i, index in enumerate(indices):
                if index >= len(log_p):
                    # the index is out of bounds... just extrapolate for now
                    tmp[s, i] += tmp[s, i - 1] + np.mean(log_p)
                else:
                    tmp[s, i] += log_p[index]

        # this line performs the averaging in log space (thus we need logsumexp)
        result[r] = logsumexp(tmp) - np.log(num_s)

    return result


def log_p(traj_lengths, signal, response, reactions):
    network = stochastic_sim.create_reaction_network(**reactions)
    
    sc = signal['components']
    st = signal['timestamps']
    rc = response['components']
    rt = response['timestamps']
    events = response['reaction_events']
    
    return log_likelihood(traj_lengths, sc, st, rc, rt, events, network, dtype=np.double)
    
    
