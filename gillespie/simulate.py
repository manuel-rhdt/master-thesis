import numpy

import accelerate
from . import stochastic_sim
from .trajectory_array import TrajectoryArray


def simulate(count, length, reactions, ext_trajectory=None, initial_value=0):
    timestamps = numpy.zeros((count, length), dtype=numpy.double)
    trajectory = numpy.zeros((count, 1, length), dtype=numpy.uint32)
    reaction_events = numpy.zeros((count, length - 1), dtype=numpy.uint32)

    trajectory[..., 0] = initial_value

    if ext_trajectory is not None:
        ext_timestamps = ext_trajectory.timestamps
        ext_components = ext_trajectory.components
    else:
        ext_timestamps = None
        ext_components = None
    
    accelerate.simulate_trajectories(
        timestamps,
        trajectory,
        reaction_events,
        reactions,
        ext_timestamps,
        ext_components,
    )

    return TrajectoryArray(timestamps, trajectory, reaction_events)


def simulate_until(until, reactions, ext_trajectory=None, initial_value=0):
    network = stochastic_sim.create_reaction_network(**reactions)

    components = numpy.zeros(2, dtype=numpy.int16)
    components[0] = ext_trajectory.components[0, 0]
    components[1] = initial_value

    if ext_trajectory is not None:
        ext_timestamps = ext_trajectory.timestamps
        ext_components = ext_trajectory.components
    else:
        ext_timestamps = None
        ext_components = None

    stochastic_sim.simulate_until_one(
        until, components, network, ext_timestamps, ext_components
    )

    return components[1:]


def reverse_trajectory(trajectory):
    trajectory.timestamps = (
        -trajectory.timestamps[..., ::-1]
        + trajectory.timestamps[..., -1, numpy.newaxis]
    )
    trajectory.components = trajectory.components[..., ::-1]

    if trajectory.reaction_events is not None:
        trajectory.reaction_events = trajectory.reaction_events[..., ::-1]

    return trajectory
