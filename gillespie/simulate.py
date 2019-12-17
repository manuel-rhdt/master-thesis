import numpy

from . import stochastic_sim


def simulate(length, reactions, ext_trajectory=None):
    timestamps = numpy.zeros(length, dtype=numpy.double)
    trajectory = numpy.zeros((1, length), dtype=numpy.uint16)
    reaction_events = numpy.zeros(length - 1, dtype=numpy.uint8)
    network = stochastic_sim.create_reaction_network(**reactions)

    if ext_trajectory is not None:
        ext_timestamps = ext_trajectory['timestamps']
        ext_components = ext_trajectory['components']
    else:
        ext_timestamps = None
        ext_components = None

    stochastic_sim.simulate_one(
        timestamps, trajectory, reaction_events, network, ext_timestamps, ext_components
    )
    return {
        "timestamps": timestamps,
        "components": trajectory,
        "reaction_events": reaction_events,
    }
