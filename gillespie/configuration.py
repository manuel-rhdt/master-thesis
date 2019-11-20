""" Loads the configuration file
"""

import os
import string

import toml
from jinja2 import Environment, FileSystemLoader

from . import stochastic_sim

CONF = None


def get(path="configuration.toml"):
    global CONF
    if CONF is not None:
        return CONF
    CONF = load(path)
    return CONF


def evaluate_envvar(func):
    val = os.environ[func["envvar"]]
    try:
        return int(val)
    except ValueError:
        try:
            return float(val)
        except ValueError:
            return val


def evaluate_expandvars(func):
    return string.Template(func["expandvars"]).substitute(os.environ)


def evaluate_expression(expr):
    functions = [
        {"key": "envvar", "func": evaluate_envvar},
        {"key": "expandvars", "func": evaluate_expandvars},
    ]

    for fn in functions:
        if fn["key"] in expr:
            return fn["func"](expr)

    # we did not find any matching function so we just return the
    # expression unchanged
    return expr


def expand(dictionary):
    """Perform a recursive evaluation of configuration expression.

    Arguments:
        dictionary {dict} -- The dictionary whose entries to expand
    """
    for key, value in dictionary.items():
        if type(value) == dict:
            expand(value)
            dictionary[key] = evaluate_expression(value)


def process_template(path):
    directory, filename = os.path.split(path)
    env = Environment(loader=FileSystemLoader(directory), trim_blocks=True)
    template = env.get_template(filename)
    return template.render(env=os.environ)


def load(path):
    loaded = process_template(path)
    conf_string = loaded
    conf = toml.loads(loaded)
    conf["original"] = conf_string
    expand(conf)
    return conf


def read_reactions(conf=None):
    if conf is None:
        conf = get()

    name_index_table = {}

    for comp in conf["signal"]["components"]:
        name_index_table[comp] = len(name_index_table)

    k = []
    reactants = []
    products = []
    for reaction in conf["signal"]["reactions"]:
        k.append(reaction["k"])
        reactants.append([])
        products.append([])
        for reactant in reaction["reactants"]:
            reactants[-1].append(name_index_table[reactant])
        for product in reaction["products"]:
            products[-1].append(name_index_table[product])

    signal_network = stochastic_sim.create_reaction_network(k, reactants, products)

    for comp in conf["response"]["components"]:
        name_index_table[comp] = len(name_index_table)

    k = []
    reactants = []
    products = []
    for reaction in conf["response"]["reactions"]:
        k.append(reaction["k"])
        reactants.append([])
        products.append([])
        for reactant in reaction["reactants"]:
            reactants[-1].append(name_index_table[reactant])
        for product in reaction["products"]:
            products[-1].append(name_index_table[product])

    response_network = stochastic_sim.create_reaction_network(k, reactants, products)

    return signal_network, response_network