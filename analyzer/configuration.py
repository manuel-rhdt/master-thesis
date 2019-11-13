""" Loads the configuration file
"""

import toml

from string import Template
import os

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
