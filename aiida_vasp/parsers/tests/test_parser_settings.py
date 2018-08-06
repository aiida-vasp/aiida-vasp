# pylint: disable=unused-import,too-few-public-methods,missing-docstring,no-self-use,no-member
"""Test the Settings utils Module."""

from aiida_vasp.parsers.parser_settings import ParserSettings
from aiida_vasp.parsers.vasp import DEFAULT_OPTIONS

SETTINGS = {
    'add_wavecar': True,
}


def test_settings_utils():

    settings = ParserSettings(SETTINGS)
    assert len(settings.nodes) == 1 and 'wavecar' in settings.nodes
    settings = ParserSettings(SETTINGS, DEFAULT_OPTIONS)
    assert 'wavecar' in settings.nodes
