import numpy

from . import stochastic_sim


def simulate(count, length, reactions, ext_trajectory=None, initial_value=0):
    timestamps = numpy.zeros((count, length), dtype=numpy.double)
    trajectory = numpy.zeros((count, 1, length), dtype=numpy.uint16)
    reaction_events = numpy.zeros((count, length - 1), dtype=numpy.uint8)
        
    trajectory[...,0] = initial_value
    network = stochastic_sim.create_reaction_network(**reactions)

    if ext_trajectory is not None:
        ext_timestamps = ext_trajectory['timestamps']
        ext_components = ext_trajectory['components']
        stochastic_sim.simulate_ext(timestamps, trajectory, reaction_events, network, ext_timestamps, ext_components)
    else:
        stochastic_sim.simulate(timestamps, trajectory, reaction_events, network)

    return {
        "timestamps": timestamps,
        "components": trajectory,
        "reaction_events": reaction_events,
    }


def simulate_until(until, reactions, ext_trajectory=None, initial_value=0):
    network = stochastic_sim.create_reaction_network(**reactions)
    
    components = np.zeros(2, dtype=np.int16)
    components[0] = ext_trajectory["components"][0, 0]
    components[1] = initial_value

    if ext_trajectory is not None:
        ext_timestamps = ext_trajectory['timestamps']
        ext_components = ext_trajectory['components']
    else:
        ext_timestamps = None
        ext_components = None

    stochastic_sim.simulate_until_one(
        until, components, network, ext_timestamps, ext_components
    )
    
    return components[1:]



def reverse_trajectory(trajectory):
    trajectory["timestamps"] = (
        -trajectory["timestamps"][..., ::-1] + trajectory["timestamps"][..., -1, numpy.newaxis]
    )
    trajectory["components"] = trajectory["components"][..., ::-1]
    
    if "reaction_events" in trajectory:
        trajectory["reaction_events"] = trajectory["reaction_events"][..., ::-1]
    
    return trajectory
