import numpy
import numpy.random

import numba
from numba import njit, prange, float64, int32, jitclass

spec = [
    ('k', float64[:]),
    ('reactants', int32[:, :]),
    ('products', int32[:, :])
]


@jitclass(spec)
class ReactionNetwork(object):
    def __init__(self, num_reactions):
        self.k = numpy.zeros(num_reactions)
        self.reactants = numpy.zeros((num_reactions, 1), dtype=numpy.int32)
        self.products = numpy.zeros((num_reactions, 1), dtype=numpy.int32)

    @property
    def size(self):
        return self.k.size


@njit(fastmath=True)
def calc_propensities(components, propensities, reactions):
    for n in range(reactions.size):
        propensities[n] = reactions.k[n]
        for reactant in reactions.reactants[n]:
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
def update_components(selected_reaction, components, reactions):
    for reactant in reactions.reactants[selected_reaction]:
        if reactant >= 0:
            components[reactant] -= 1
    for product in reactions.products[selected_reaction]:
        if product >= 0:
            components[product] += 1


@njit(fastmath=True)
def timestep_generate(components, ext_timestamps, ext_components, reactions):
    num_ext_comps = len(ext_components) if ext_components is not None else 0
    ext_length = len(ext_timestamps) if ext_timestamps is not None else 0
    current_time = 0.0
    ext_progress = 0
    propensities = numpy.zeros(reactions.size)
    random_variate = -numpy.log(numpy.random.random_sample())
    while True:
        if ext_progress >= ext_length:
            next_ext_timestamp = numpy.Inf
        elif ext_timestamps is not None:
            next_ext_timestamp = ext_timestamps[ext_progress]

        calc_propensities(components, propensities, reactions)
        total_propensity = numpy.sum(propensities)

        perform_reaction, timestep = try_propagate_time(
            random_variate, current_time, next_ext_timestamp, total_propensity)
        current_time += timestep

        if perform_reaction:
            # Update components
            selected_reaction = select_reaction(propensities)
            update_components(selected_reaction, components, reactions)

            yield current_time, selected_reaction
            random_variate = -numpy.log(numpy.random.random_sample())
        else:
            random_variate -= timestep * total_propensity
            # update the external trajectory
            ext_progress += 1
            for i in range(num_ext_comps):
                if ext_components is not None:
                    components[i] = ext_components[i][min(
                        ext_progress, len(ext_components[i]) - 1)]


@njit(fastmath=True)
def simulate_one(timestamps, trajectory, reaction_events, reactions, ext_timestamps=None, ext_components=None):
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

    num_comps = len(trajectory)
    num_ext_comps = len(ext_components) if ext_components is not None else 0

    # set the trajectory's initial values
    timestamps[0] = 0.0

    # define values used during iteration
    progress = 0
    components = numpy.zeros(num_comps + num_ext_comps)
    for i in range(num_ext_comps):
        if ext_components is not None:
            components[i] = ext_components[i][0]
    for comp in range(num_comps):
        components[comp + num_ext_comps] = trajectory[comp][0]

    for time, selected_reaction in timestep_generate(components, ext_timestamps, ext_components, reactions):
        # update trajectory
        progress += 1
        timestamps[progress] = time
        reaction_events[progress - 1] = selected_reaction
        trajectory[:, progress] = components[num_ext_comps:]

        if progress >= length - 1:
            break


@njit(fastmath=True)
def simulate_until_one(until, initial_values, reactions, ext_timestamps=None, ext_components=None):
    prev_time = 0.0
    components = initial_values
    prev_components = components
    for time, _ in timestep_generate(components, ext_timestamps, ext_components, reactions):
        if time >= until:
            return prev_time, prev_components

        prev_time = time
        prev_components = components


@njit(parallel=True, fastmath=True, cache=True)
def simulate_until(until, initial_values, reactions, ext_timestamps=None, ext_components=None):
    num_trajectories, = until.shape

    for r in prange(num_trajectories):
        if ext_components is not None:
            simulate_until_one(until[r], initial_values[r], reactions,
                               ext_timestamps[r], ext_components[r])
        else:
            simulate_until_one(until[r], initial_values[r], reactions)


@njit(parallel=True, fastmath=True, cache=True)
def simulate(timestamps, trajectory, reaction_events, reactions, ext_timestamps=None, ext_components=None):
    assert len(timestamps.shape) == 2
    assert len(trajectory.shape) == 3
    assert len(reaction_events.shape) == 2

    assert timestamps.shape[0] == trajectory.shape[0] == reaction_events.shape[0]

    for r in prange(timestamps.shape[0]):
        if ext_components is not None:
            simulate_one(timestamps[r], trajectory[r],
                         reaction_events[r], reactions, ext_timestamps[r], ext_components[r])
        else:
            simulate_one(timestamps[r], trajectory[r],
                         reaction_events[r], reactions)
