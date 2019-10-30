

def propagate_time(reactions, components):

    for reaction in reactions:
        k = reaction['k']

        # propensity
