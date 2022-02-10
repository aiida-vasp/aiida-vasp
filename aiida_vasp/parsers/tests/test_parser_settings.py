# pylint: disable=unused-import,too-few-public-methods,missing-docstring,no-self-use,no-member
"""Test the Settings utils Module."""

from aiida_vasp.parsers.settings import ParserSettings
from aiida_vasp.parsers.vasp import DEFAULT_SETTINGS

SETTINGS = {'add_wavecar': True, 'critical_notifications': {'add_brmix': False}}


def test_settings_utils():

    settings = ParserSettings(SETTINGS)
    assert len(settings.output_nodes_dict) == 1 and 'wavecar' in settings.output_nodes_dict
    settings = ParserSettings(SETTINGS, DEFAULT_SETTINGS)
    assert 'wavecar' in settings.output_nodes_dict
    assert 'kpoints' not in settings.output_nodes_dict
    assert 'misc' in settings.output_nodes_dict

    assert 'psmaxn' not in settings.critical_notifications_to_check
    assert 'brmix' not in settings.critical_notifications_to_check
    assert 'fock_acc' in settings.critical_notifications_to_check
    assert 'rhosyg' in settings.critical_notifications_to_check
