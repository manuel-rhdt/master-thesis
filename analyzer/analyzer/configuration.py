""" Loads the configuration file
"""

import toml

CONF = None


def get(path='configuration.toml'):
    global CONF
    if CONF is not None:
        return CONF
    CONF = load(path)
    return CONF


def load(path):
    return toml.load(path)
