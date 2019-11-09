import numpy
import numpy.random

import numba
from numba import njit, prange


@njit(fastmath=True)
def calc_propensities(components, reaction_k, reaction_reactants, propensities):
    for n in range(len(reaction_k)):
        propensities[n] = reaction_k[n]
        for reactant in reaction_reactants[n]:
            if reactant >= 0:
                propensities[n] *= components[reactant]

    return propensities


@njit(fastmath=True)
def try_propagate_time(random_variate, timestamp, next_ext_timestamp, total_propensity):
    """ Returns `True` if a reaction should be executed before continuing. 

    The second return value is the amount of time that was propagated.

    The random variate should be the negative log of a uniform value sampled from the interval [0.0, 1.0).
    """
    max_time_step = next_ext_timestamp - timestamp

    if max_time_step <= 0.0:
        return (False, 0.0)

    assert total_propensity >= 0.0

    if max_time_step * total_propensity < random_variate:
        # no progress was made
        return (False, max_time_step)
    else:
        time_step = random_variate / total_propensity
        return (True, time_step)


@njit(fastmath=True)
def select_reaction(propensities):
    r = numpy.random.random_sample()
    propensities = numpy.asarray(propensities)

    probabilities = propensities / numpy.sum(propensities)

    selected_reaction = 0
    acc = probabilities[0]
    while r > acc and selected_reaction + 1 < len(probabilities):
        selected_reaction += 1
        acc += probabilities[selected_reaction]

    return selected_reaction


@njit(fastmath=True)
def update_components(selected_reaction, components, reaction_reactants, reaction_products):
    for reactant in reaction_reactants[selected_reaction]:
        if reactant >= 0:
            components[reactant] -= 1
    for product in reaction_products[selected_reaction]:
        if product >= 0:
            components[product] += 1


@njit(fastmath=True)
def simulate_one(timestamps, trajectory, ext_components, ext_timestamps, reaction_k, reaction_reactants, reaction_products, reaction_events):
    """Simulate one trajectory

    Arguments:
        timestamps {numpy.ndarray} -- [description]
        initial_components {[type]} -- [description]
        signal_components {[type]} -- [description]
        signal_timestamps {[type]} -- [description]
        response_components {[type]} -- [description]
        reaction_k {[type]} -- [description]
        reaction_reactants {[type]} -- [description]
        reaction_products {[type]} -- [description]
        reaction_events {[type]} -- [description]

    Returns:
        [type] -- [description]
    """
    length = len(timestamps)
    if length < 2:
        return
    random_variates = numpy.array(
        [numpy.random.random_sample() for _ in range(length)])
    random_variates = -numpy.log(random_variates)

    num_comps = len(trajectory)
    num_ext_comps = len(ext_components) if ext_components is not None else 0
    ext_length = len(ext_timestamps) if ext_timestamps is not None else 0

    # set the trajectory's initial values
    timestamps[0] = 0.0
    trajectory[:, progress] = components[num_ext_comps:]

    # define values used during iteration
    progress = 0
    ext_progress = 0
    components = numpy.zeros(num_comps + num_ext_comps)
    for i in range(num_ext_comps):
        if ext_components is not None:
            components[i] = ext_components[i][ext_progress]
    for comp in range(num_comps):
        components[comp + num_ext_comps] = trajectory[comp][progress]
    propensities = numpy.zeros(reaction_k.shape)
    accumulated_time = 0.0
    current_time = 0.0

    while progress + 1 < length:
        if ext_progress >= ext_length:
            next_ext_timestamp = numpy.Inf
        elif ext_timestamps is not None:
            next_ext_timestamp = ext_timestamps[ext_progress]

        calc_propensities(components, reaction_k,
                          reaction_reactants, propensities)
        total_propensity = numpy.sum(propensities)

        perform_reaction, timestep = try_propagate_time(
            random_variates[progress], current_time, next_ext_timestamp, total_propensity)
        
        current_time += timestep

        if perform_reaction:
            # Update components
            selected_reaction = select_reaction(propensities)
            update_components(selected_reaction, components,
                              reaction_reactants, reaction_products)

            # update trajectory
            progress += 1
            timestamps[progress] = current_time
            reaction_events[progress - 1] = selected_reaction
            trajectory[:, progress] = components[num_ext_comps:]
        else:
            random_variates[progress] -= timestep * total_propensity
            # update the external trajectory
            ext_progress += 1
            for i in range(num_ext_comps):
                if ext_components is not None:
                    components[i] = ext_components[i][min(
                        ext_progress, len(ext_components[i]) - 1)]


@njit(parallel=True, fastmath=True, cache=True)
def simulate(timestamps, trajectory, ext_components, ext_timestamps, reaction_k, reaction_reactants, reaction_products, reaction_events):
    assert len(timestamps.shape) == 2
    assert len(trajectory.shape) == 3
    assert len(reaction_events.shape) == 2

    assert timestamps.shape[0] == trajectory.shape[0] == reaction_events.shape[0]

    for r in prange(timestamps.shape[0]):
        if ext_components is not None:
            simulate_one(timestamps[r], trajectory[r], ext_components[r], ext_timestamps[r],
                         reaction_k, reaction_reactants, reaction_products, reaction_events[r])
        else:
            simulate_one(timestamps[r], trajectory[r], None, None,
                         reaction_k, reaction_reactants, reaction_products, reaction_events[r])
