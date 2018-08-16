"""Module defining sets of FileParsers to be used by the VaspParser"""

from aiida_vasp.io.doscar import DosParser
from aiida_vasp.io.eigenval import EigParser
from aiida_vasp.io.kpoints import KpParser
from aiida_vasp.io.outcar import OutcarParser
from aiida_vasp.io.vasprun import VasprunParser
from aiida_vasp.io.chgcar import ChgcarParser
from aiida_vasp.io.wavecar import WavecarParser
from aiida_vasp.io.poscar import PoscarParser

FILE_PARSER_SETS = {
    'default': {
        'DOSCAR': {
            'parser_class': DosParser,
            'is_critical': False,
            'status': 'Unknown'
        },
        'EIGENVAL': {
            'parser_class': EigParser,
            'is_critical': False,
            'status': 'Unknown'
        },
        'IBZKPT': {
            'parser_class': KpParser,
            'is_critical': False,
            'status': 'Unknown'
        },
        'OUTCAR': {
            'parser_class': OutcarParser,
            'is_critical': False,
            'status': 'Unknown'
        },
        'vasprun.xml': {
            'parser_class': VasprunParser,
            'is_critical': True,
            'status': 'Unknown'
        },
        'CHGCAR': {
            'parser_class': ChgcarParser,
            'is_critical': False,
            'status': 'Unknown'
        },
        'WAVECAR': {
            'parser_class': WavecarParser,
            'is_critical': False,
            'status': 'Unknown'
        },
        'CONTCAR': {
            'parser_class': PoscarParser,
            'is_critical': False,
            'status': 'Unknown'
        },
    },
}


class ParserSettings(object):
    """
    Settings object for the VaspParser.

    :param settings: Dict with the 'parser_settings'.
    :param default_settings: Dict with default settings.

    This provides the following properties to other components of the VaspParser:

        * nodes: A list with all requested output nodes.

        * parser_definitions: A Dict with the FileParser definitions.
    """

    def __init__(self, settings, default_settings=None):

        if settings is None:
            settings = {}
        self._settings = settings
        if default_settings:
            self.update_with(default_settings)

        self.nodes = []
        self.set_nodes()

        self.parser_definitions = {}
        self.set_parser_definitions(self._settings.get('file_parser_set'))

    def set_nodes(self):
        """Set the 'nodes' card of a settings object."""
        # Find all the nodes, that should be added.
        nodes = []
        for key, value in self._settings.items():
            if not key.startswith('add_'):
                # only keys starting with 'add_' are relevant as nodes.
                continue
            if not value:
                # The quantity should not be added.
                continue
            nodes.append(key[4:])

        self.nodes = nodes

    def update_with(self, update_dict):
        """Selectively update keys from one Dictionary to another."""
        for key, value in update_dict.items():
            if key not in self._settings:
                self._settings[key] = value

    def get(self, item, default=None):
        return self._settings.get(item, default)

    def set_parser_definitions(self, file_parser_set='default'):
        """Load the parser definitions."""
        from copy import deepcopy

        if file_parser_set not in FILE_PARSER_SETS:
            return
        for file_name, parser_dict in FILE_PARSER_SETS.get(file_parser_set).items():
            self.parser_definitions[file_name] = deepcopy(parser_dict)
