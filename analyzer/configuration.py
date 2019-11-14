""" Loads the configuration file
"""

import toml

from string import Template
import os

from . import stochastic_sim

CONF = None


def get(path='configuration.toml'):
    global CONF
    if CONF is not None:
        return CONF
    CONF = load(path)
    return CONF


def expand(dictionary):
    """Perform a recursive evaluation of configuration functions.

    Arguments:
        dictionary {dict} -- The dictionary to evaluate
    """
    for key, value in dictionary.items():
        if type(value) == dict:
            expand(value)
            if 'envvar' in value:
                dictionary[key] = int(os.environ[value['envvar']])
            if 'expandvars' in value:
                dictionary[key] = Template(
                    value['expandvars']).substitute(os.environ)


def load(path):
    conf = toml.load(path)
    expand(conf)
    return conf


def read_reactions(conf=get()):
    name_index_table = {}

    for comp in conf['signal']['components']:
        name_index_table[comp] = len(name_index_table)

    k = []
    reactants = []
    products = []
    for reaction in conf['signal']['reactions']:
        k.append(reaction['k'])
        reactants.append([])
        products.append([])
        for reactant in reaction['reactants']:
            reactants[-1].append(name_index_table[reactant])
        for product in reaction['products']:
            products[-1].append(name_index_table[product])

    signal_network = stochastic_sim.create_reaction_network(
        k, reactants, products)

    for comp in conf['response']['components']:
        name_index_table[comp] = len(name_index_table)

    k = []
    reactants = []
    products = []
    for reaction in conf['response']['reactions']:
        k.append(reaction['k'])
        reactants.append([])
        products.append([])
        for reactant in reaction['reactants']:
            reactants[-1].append(name_index_table[reactant])
        for product in reaction['products']:
            products[-1].append(name_index_table[product])

    response_network = stochastic_sim.create_reaction_network(
        k, reactants, products)

    return signal_network, response_network
