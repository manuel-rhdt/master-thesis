import numpy
import numpy.random


def calc_propensities(components, reaction_k, reaction_reactants):
    propensities = numpy.zeros(reaction_k.shape)
    for n in range(len(reaction_k)):
        propensities[n] = reaction_k[n]
        for reactant in reaction_reactants[n]:
            propensities[n] *= components[reactant]

    return propensities


def try_propagate_time(progress, random_variates, timestamps, next_ext_timestamp, components, reaction_k, reaction_reactants):
    """ Returns `True` if time was successfully propagated
    """
    max_time_step = next_ext_timestamp - timestamps[progress]

    if max_time_step <= 0.0:
        return False

    propensities = calc_propensities(
        components, reaction_k, reaction_reactants)
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


def simulate(length, initial_components, signal_components, signal_timestamps, response_components, reaction_k, reaction_reactants, reaction_events):
    random_variates = -numpy.log(numpy.random.random_sample(size=length))
    timestamps = numpy.zeros(length)

    components = initial_components

    progress = 0
    ext_progress = 0

    result = numpy.empty((len(response_components), length))
    timestamps = numpy.zeros(length)
    while progress < length - 1:
        next_ext_timestamp = signal_timestamps[ext_progress]

        if try_propagate_time(progress, random_variates, timestamps, next_ext_timestamp,
                              components, reaction_k, reaction_reactants):
            progress += 1
            # TODO: Update components
            result[:, progress] = components[len(signal_components):]
        else:
            # update the external trajectory
            ext_progress += 1
            for i, comp in enumerate(signal_components):
                components[i] = comp[ext_progress]
