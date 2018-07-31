"""Module Containing helper functions for settings."""
from aiida_vasp.utils.extended_dicts import DictWithAttributes


def create_new_settings(settings_dict, default_settings=None):
    """
    Create a new settings object for the VaspParser.

    :param settings_dict: the 'settings' dict from the VaspCalculation.
    :param default_settings: dictionary with default settings.
    :return: A DictWithAttributes object.
    """

    if settings_dict is None:
        settings_dict = {}

    settings = DictWithAttributes(settings_dict)

    if default_settings:
        update_settings_with(settings, default_settings)

    set_nodes(settings)

    return settings


def set_nodes(settings):
    """Set the 'nodes' card of a settings object."""
    # Find all the nodes, that should be added.
    nodes = []
    for key, value in settings.items():
        if not key.startswith('add_'):
            # only keys starting with 'add_' are relevant as nodes.
            continue
        if not value:
            # The quantity should not be added.
            continue
        nodes.append(key[4:])

    settings['nodes'] = nodes


def update_settings_with(settings, update_dict):
    """Selectively update keys from one Dictionary to another."""
    for key, value in update_dict.items():
        if key not in settings:
            settings[key] = value
