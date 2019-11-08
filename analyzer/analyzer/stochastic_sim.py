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
def try_propagate_time(progress, random_variates, timestamps, next_ext_timestamp, components, reaction_k, reaction_reactants, propensities):
    """ Returns `True` if time was successfully propagated. 

    If successfull the new timestamp will be saved in `timestamps`.
    """
    max_time_step = next_ext_timestamp - timestamps[progress]

    if max_time_step <= 0.0:
        return False

    calc_propensities(components, reaction_k, reaction_reactants, propensities)
    total_propensity = numpy.sum(propensities)

    assert total_propensity >= 0.0

    if max_time_step * total_propensity < random_variates[progress]:
        random_variates[progress] -= max_time_step * total_propensity
        timestamps[progress + 1] += max_time_step
        # no progress was made
        return False
    else:
        time_step = random_variates[progress] / total_propensity
        timestamps[progress + 1] += time_step
        return True


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
            components[reactant] -= 1.0
    for product in reaction_products[selected_reaction]:
        if product >= 0:
            components[product] += 1.0


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

    progress = 0
    ext_progress = 0
    components = numpy.zeros(num_comps + num_ext_comps)

    # use the initial values
    timestamps[:] = 0.0
    for i in range(num_ext_comps):
        if ext_components is not None:
            components[i] = ext_components[i][ext_progress]
    for comp in range(num_comps):
        components[comp + num_ext_comps] = trajectory[comp][progress]

    trajectory[:, progress] = components[num_ext_comps:]

    propensities = numpy.zeros(reaction_k.shape)

    while progress + 1 < length:
        if ext_progress >= ext_length:
            next_ext_timestamp = numpy.Inf
        elif ext_timestamps is not None:
            next_ext_timestamp = ext_timestamps[ext_progress]

        if try_propagate_time(progress, random_variates, timestamps, next_ext_timestamp,
                              components, reaction_k, reaction_reactants, propensities):
            progress += 1

            # Update components
            selected_reaction = select_reaction(propensities)
            reaction_events[progress - 1] = selected_reaction
            update_components(selected_reaction, components,
                              reaction_reactants, reaction_products)

            trajectory[:, progress] = components[num_ext_comps:]

            if progress + 1 < length:
                timestamps[progress + 1] = timestamps[progress]
        else:
            # update the external trajectory
            ext_progress += 1
            for i in range(num_ext_comps):
                if ext_components is not None:
                    components[i] = ext_components[i][min(
                        ext_progress, len(ext_components[i]) - 1)]


@njit(parallel=True, fastmath=True)
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
