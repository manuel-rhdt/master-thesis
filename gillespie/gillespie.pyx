cimport numpy
cimport cython
import numpy.random
from numpy.math cimport INFINITY

cdef expand_dims(array, int dims=2):
    if array is None:
        return None

    array_dim = array.ndim
    if array_dim > dims:
        raise Exception("Array has too many dimensions!")

    shape = array.shape
    return array.reshape(((1,) * (dims - array_dim) + shape))

cdef expand_3d(array):
    return expand_dims(array, 3)

cdef class ReactionNetwork():
    cdef public numpy.double_t[:] k
    cdef public numpy.int32_t[:, :] reactants
    cdef public numpy.int32_t[:, :] products

    def __init__(self, num_reactions):
        self.k = numpy.zeros(num_reactions, dtype=numpy.double)
        self.reactants = numpy.zeros((num_reactions, 1), dtype=numpy.int32)
        self.products = numpy.zeros((num_reactions, 1), dtype=numpy.int32)

    cdef unsigned int size(self) nogil:
        return len(self.k)


def create_reaction_network(k, reactants, products):
    network = ReactionNetwork(len(k))
    network.k = numpy.asarray(k, dtype=numpy.double)
    max_num_reactants = max(len(react) for react in reactants)
    network.reactants = numpy.full((len(k), max_num_reactants), -1, dtype=numpy.int32)
    for i, react in enumerate(reactants):
        for j, r in enumerate(react):
            network.reactants[i, j] = r

    max_num_products = max(len(prod) for prod in products)
    network.products = numpy.full((len(k), max_num_products), -1, dtype=numpy.int32)
    for i, prod in enumerate(products):
        for j, p in enumerate(prod):
            network.products[i, j] = p

    return network


@cython.boundscheck(False)
cdef void calc_propensities(int[:] components, 
        double[:] propensities, 
        ReactionNetwork reactions) nogil:
    for n in range(reactions.size()):
        propensities[n] = reactions.k[n]
        for r in range(len(reactions.reactants[n])):
            if reactions.reactants[n, r] >= 0:
                propensities[n] *= <double>components[reactions.reactants[n][r]]

cdef (bint, double) try_propagate_time(double random_variate, 
        double timestamp, 
        double next_ext_timestamp, 
        double total_propensity):
    """
    Returns `True` if a reaction should be executed before continuing.

    The second return value is the amount of time that was propagated.

    The random variate should be the negative log of a uniform value sampled from the
    interval [0.0, 1.0).
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


cdef unsigned int select_reaction(double[:] propensities):
    r = numpy.random.random_sample()
    propensities = numpy.asarray(propensities)

    probabilities = propensities / numpy.sum(propensities)

    selected_reaction = 0
    acc = probabilities[0]
    while r > acc and selected_reaction + 1 < len(probabilities):
        selected_reaction += 1
        acc += probabilities[selected_reaction]

    return selected_reaction


cdef class StochasticSim():
    cdef public double current_time
    cdef public unsigned int ext_progress
    cdef public double[:] propensities
    cdef public int[:] components
    cdef public double[:] ext_timestamps
    cdef public int[:,:] ext_components
    cdef public ReactionNetwork reactions

    def __init__(
        self, components, propensities, ext_timestamps, ext_components, reactions
    ):
        self.current_time = 0.0
        self.ext_progress = 0
        self.propensities = propensities
        self.components = components
        self.ext_timestamps = ext_timestamps
        self.ext_components = ext_components
        self.reactions = reactions

    cdef int ext_len(self) nogil:
        return len(self.ext_timestamps) if self.ext_timestamps is not None else 0

    cdef int num_ext_comps(self) nogil:
        return len(self.ext_components) if self.ext_components is not None else 0

    @cython.boundscheck(False)
    cdef double next_ext_timestamp(self) nogil:
        if self.ext_progress >= self.ext_len():
            return INFINITY
        elif self.ext_timestamps is not None:
            return self.ext_timestamps[self.ext_progress]
        else:
            raise RuntimeError("no external trajectory")

    @cython.boundscheck(False)
    cdef int next_ext_comp(self, int i) nogil:
        if self.ext_components is not None:
            index = min(self.ext_progress, len(self.ext_components[i]) - 1)
            return self.ext_components[i][index]
        else:
            raise RuntimeError("no external trajectory")

    @cython.boundscheck(False)
    cdef void update_components(self, int selected_reaction) nogil:
        cdef int i = 0
        cdef numpy.int32_t[:] reactants = self.reactions.reactants[selected_reaction]
        for i in range(len(reactants)):
            if reactants[i] >= 0:
                self.components[reactants[i]] -= 1

        cdef numpy.int32_t[:] products = self.reactions.products[selected_reaction]
        for i in range(len(products)):
            if products[i] >= 0:
                self.components[products[i]] += 1

    cdef (double, int) propagate_time(self) nogil:
        cdef double random_variate = -numpy.log(numpy.random.random_sample())

        while True:
            calc_propensities(self.components, self.propensities, self.reactions)
            total_propensity = numpy.sum(self.propensities)

            perform_reaction, timestep = try_propagate_time(
                random_variate,
                self.current_time,
                self.next_ext_timestamp(),
                total_propensity,
            )
            self.current_time += timestep

            if perform_reaction:
                selected_reaction = select_reaction(self.propensities)
                self.update_components(selected_reaction)
                return self.current_time, selected_reaction
            else:
                random_variate -= timestep * total_propensity
                # update the external trajectory
                self.ext_progress += 1
                for i in range(self.num_ext_comps()):
                    if self.ext_components is not None:
                        self.components[i] = self.next_ext_comp(i)


def simulate_one(
    timestamps,
    trajectory,
    reaction_events,
    reactions,
    ext_timestamps=None,
    ext_components=None,
    components=None,
    propensities=None,
):
    """Simulate one trajectory
    """
    length = len(timestamps)
    if length < 2:
        return

    num_comps, _ = trajectory.shape
    num_ext_comps, _ = ext_components.shape if ext_components is not None else (0, 0)

    # set the trajectory's initial values
    timestamps[0] = 0.0

    # define values used during iteration
    progress = 0
    propensities = numpy.zeros(reactions.size) if propensities is None else propensities
    components = (
        numpy.zeros(num_comps + num_ext_comps, dtype=numpy.uint32)
        if components is None
        else components
    )
    for i in range(num_ext_comps):
        if ext_components is not None:
            components[i] = ext_components[i][0]
    for comp in range(num_comps):
        components[comp + num_ext_comps] = trajectory[comp][0]

    sim = StochasticSim(
        components, propensities, ext_timestamps, ext_components, reactions
    )
    while (progress + 1) < length:
        progress += 1
        time, selected_reaction = sim.propagate_time()
        timestamps[progress] = time
        reaction_events[progress - 1] = selected_reaction
        trajectory[:, progress] = sim.components[num_ext_comps:]


def simulate_until_one(
    until, initial_values, reactions, ext_timestamps=None, ext_components=None
):
    prev_time = 0.0
    components = initial_values
    prev_components = components

    sim = StochasticSim(components, ext_timestamps, ext_components, reactions)
    while True:
        time, _ = sim.propagate_time()
        if time >= until:
            return prev_time, prev_components

        prev_time = time
        prev_components = sim.components


def simulate_until(
    until, initial_values, reactions, ext_timestamps=None, ext_components=None
):
    (num_trajectories,) = until.shape

    for r in range(num_trajectories):
        if ext_components is not None:
            simulate_until_one(
                until[r],
                initial_values[r],
                reactions,
                ext_timestamps[r],
                ext_components[r],
            )
        else:
            simulate_until_one(until[r], initial_values[r], reactions)


def simulate(timestamps, trajectory, reaction_events, reactions):
    timestamps = numpy.atleast_2d(timestamps)
    reaction_events = numpy.atleast_2d(reaction_events)
    trajectory = expand_3d(trajectory)

    assert timestamps.shape[0] == trajectory.shape[0] == reaction_events.shape[0]

    num_r, _ = timestamps.shape

    for r in range(num_r):
        t = timestamps[r]
        c = trajectory[r]
        events = reaction_events[r]
        simulate_one(t, c, events, reactions)


def simulate_ext(
    timestamps, trajectory, reaction_events, reactions, ext_timestamps, ext_components
):
    timestamps = numpy.atleast_2d(timestamps)
    reaction_events = numpy.atleast_2d(reaction_events)
    trajectory = expand_3d(trajectory)

    ext_timestamps = numpy.atleast_2d(ext_timestamps)
    ext_components = expand_3d(ext_components)

    assert timestamps.shape[0] == trajectory.shape[0] == reaction_events.shape[0]
    assert ext_timestamps.shape[0] == ext_components.shape[0]

    num_r, num_comps, _ = trajectory.shape
    num_s, num_ext_comps, _ = ext_components.shape

    components = numpy.zeros(num_comps + num_ext_comps, dtype=numpy.uint32)
    propensities = numpy.zeros(reactions.size)

    # numpy broadcasting rules
    assert num_s == num_r or num_s == 1

    for r in range(timestamps.shape[0]):
        t = timestamps[r]
        c = trajectory[r]
        events = reaction_events[r]
        ext_t = ext_timestamps[r % num_s]
        ext_c = ext_components[r % num_s]
        simulate_one(t, c, events, reactions, ext_t, ext_c, components, propensities)
