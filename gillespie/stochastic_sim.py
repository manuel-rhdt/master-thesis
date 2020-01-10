import numba
import numpy
import numpy.random
from numba import float64, generated_jit, int32, jitclass, njit, types

spec = [("k", float64[:]), ("reactants", int32[:, :]), ("products", int32[:, :])]


@generated_jit(nopython=True)
def expand_3d(array):
    if isinstance(array, types.Array):
        if array.ndim == 1:
            return lambda array: array.reshape((1, 1, -1))
        elif array.ndim == 2:
            return lambda array: numpy.expand_dims(array, 0)
        else:
            return lambda array: array


@jitclass(spec)
class ReactionNetwork(object):
    def __init__(self, num_reactions):
        self.k = numpy.zeros(num_reactions, dtype=numpy.double)
        self.reactants = numpy.zeros((num_reactions, 1), dtype=numpy.int32)
        self.products = numpy.zeros((num_reactions, 1), dtype=numpy.int32)

    @property
    def size(self):
        return self.k.size


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


@njit
def calc_propensities(components, propensities, reactions):
    for n in range(reactions.size):
        propensities[n] = reactions.k[n]
        for reactant in reactions.reactants[n]:
            if reactant >= 0:
                propensities[n] *= components[reactant]

    return propensities


@njit
def try_propagate_time(random_variate, timestamp, next_ext_timestamp, total_propensity):
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


@njit
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


dummy = ReactionNetwork(0)
spec_ss = [
    ("current_time", types.float64),
    ("ext_progress", types.int32),
    ("propensities", types.float64[:]),
    ("components", types.uint32[:]),
    ("ext_timestamps", types.Optional(float64[:])),
    ("ext_components", types.Optional(types.uint16[:, :])),
    ("reactions", numba.typeof(dummy)),
]


@jitclass(spec_ss)
class StochasticSim(object):
    def __init__(self, components, ext_timestamps, ext_components, reactions):
        self.current_time = 0.0
        self.ext_progress = 0
        self.propensities = numpy.zeros(reactions.size)
        self.components = components
        self.ext_timestamps = ext_timestamps
        self.ext_components = ext_components
        self.reactions = reactions

    @property
    def ext_len(self):
        return len(self.ext_timestamps) if self.ext_timestamps is not None else 0

    @property
    def num_ext_comps(self):
        return len(self.ext_components) if self.ext_components is not None else 0

    def next_ext_timestamp(self):
        if self.ext_progress >= self.ext_len:
            return numpy.Inf
        elif self.ext_timestamps is not None:
            return self.ext_timestamps[self.ext_progress]
        else:
            raise RuntimeError("no external trajectory")

    def next_ext_comp(self, i):
        if self.ext_components is not None:
            index = min(self.ext_progress, len(self.ext_components[i]) - 1)
            return self.ext_components[i][index]
        else:
            raise RuntimeError("no external trajectory")

    def update_components(self, selected_reaction):
        for reactant in self.reactions.reactants[selected_reaction]:
            if reactant >= 0:
                self.components[reactant] -= 1
        for product in self.reactions.products[selected_reaction]:
            if product >= 0:
                self.components[product] += 1

    def propagate_time(self):
        random_variate = -numpy.log(numpy.random.random_sample())

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
                for i in range(self.num_ext_comps):
                    if self.ext_components is not None:
                        self.components[i] = self.next_ext_comp(i)


@njit(cache=True, nogil=True)
def simulate_one(
    timestamps,
    trajectory,
    reaction_events,
    reactions,
    ext_timestamps=None,
    ext_components=None,
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
    components = numpy.zeros(num_comps + num_ext_comps, dtype=numpy.uint32)
    for i in range(num_ext_comps):
        if ext_components is not None:
            components[i] = ext_components[i][0]
    for comp in range(num_comps):
        components[comp + num_ext_comps] = trajectory[comp][0]

    sim = StochasticSim(components, ext_timestamps, ext_components, reactions)
    while (progress + 1) < length:
        progress += 1
        time, selected_reaction = sim.propagate_time()
        timestamps[progress] = time
        reaction_events[progress - 1] = selected_reaction
        trajectory[:, progress] = sim.components[num_ext_comps:]


@njit(cache=True, nogil=True)
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


@njit(cache=True, nogil=True)
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


@njit(cache=True, nogil=True)
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


@njit(cache=True, nogil=True)
def simulate_ext(
    timestamps, trajectory, reaction_events, reactions, ext_timestamps, ext_components
):
    timestamps = numpy.atleast_2d(timestamps)
    reaction_events = numpy.atleast_2d(reaction_events)
    trajectory = expand_3d(trajectory)

    ext_timestamps = numpy.atleast_2d(ext_timestamps)
    ext_components = expand_3d(ext_components)

    assert timestamps.shape[0] == trajectory.shape[0] == reaction_events.shape[0]

    num_r, _ = timestamps.shape
    num_s, _ = ext_timestamps.shape

    # numpy broadcasting rules
    assert num_s == num_r or num_s == 1

    for r in range(timestamps.shape[0]):
        t = timestamps[r]
        c = trajectory[r]
        events = reaction_events[r]
        ext_t = ext_timestamps[r % num_s]
        ext_c = ext_components[r % num_s]
        simulate_one(t, c, events, reactions, ext_t, ext_c)
