"""Module defining sets of FileParsers to be used by the VaspParser"""

from aiida_vasp.io.doscar import DosParser
from aiida_vasp.io.eigenval import EigParser
from aiida_vasp.io.kpoints import KpParser
from aiida_vasp.io.outcar import OutcarParser
from aiida_vasp.io.vasprun import VasprunParser
from aiida_vasp.io.chgcar import ChgcarParser
from aiida_vasp.io.wavecar import WavecarParser
from aiida_vasp.io.poscar import PoscarParser
from aiida_vasp.utils.extended_dicts import DictWithAttributes

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
            'is_critical': False,
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

NODES = {
    'parameters': {
        'link_name': 'output_parameters',
        'type': 'parameter',
        'quantities': ['total_energies', 'symmetries', 'maximum_stress', 'maximum_force', 'fermi_level'],
    },
    'kpoints': {
        'link_name': 'output_kpoints',
        'type': 'array.kpoints',
        'quantities': ['kpoints'],
    },
    'structure': {
        'link_name': 'output_structure',
        'type': 'structure',
        'quantities': ['structure'],
    },
    'poscar-structure': {
        'link_name': 'output_structure',
        'type': 'structure',
        'quantities': ['poscar-structure'],
    },
    'trajectory': {
        'link_name': 'output_trajectory',
        'type': 'array.trajectory',
        'quantities': ['trajectory'],
    },
    'bands': {
        'link_name': 'output_bands',
        'type': 'array.bands',
        'quantities': ['eigenvalues', 'kpoints', 'occupations'],
    },
    'dos': {
        'link_name': 'output_dos',
        'type': 'array',
        'quantities': ['dos'],
    },
    'energies': {
        'link_name': 'output_energies',
        'type': 'array',
        'quantities': ['energies'],
    },
    'projectors': {
        'link_name': 'output_projectors',
        'type': 'array',
        'quantities': ['projectors'],
    },
    'born_charges': {
        'link_name': 'output_born_charges',
        'type': 'array',
        'quantities': ['born_charges'],
    },
    'dielectrics': {
        'link_name': 'output_dielectrics',
        'type': 'array',
        'quantities': ['dielectrics'],
    },
    'hessian': {
        'link_name': 'output_hessian',
        'type': 'array',
        'quantities': ['hessian'],
    },
    'dynmat': {
        'link_name': 'output_dynmat',
        'type': 'array',
        'quantities': ['dynmat'],
    },
    'chgcar': {
        'link_name': 'output_chgcar',
        'type': 'vasp.chargedensity',
        'quantities': ['chgcar'],
    },
    'wavecar': {
        'link_name': 'output_wavecar',
        'type': 'vasp.wavefun',
        'quantities': ['wavecar'],
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

        self.alternatives = {}
        if settings is None:
            settings = {}
        self._settings = settings
        if default_settings is not None:
            self.update_with(default_settings)

        self.nodes = {}
        self.set_nodes()

        self.parser_definitions = {}
        self.set_parser_definitions(self._settings.get('file_parser_set'))

    def set_nodes(self):
        """
        Set the 'nodes' card of a settings object.

        Nodes can be added by setting:

            'add_node_name' : value

        in 'parser_settings'. The type of value determines the mode for adding the node:

         -  bool: if True and node_name in NODES a default node will be add with default quantities.

                'add_structure': True

         -  list: if node_name in NODES a default node with updated quantities will be added.

                'add_parameters': ['efermi', 'energies', ...]

         -  dict: a custom node will be add. The dict must provide 'type' and 'quantities'. 'link_name' is optional

                'add_custom_node': {'type': 'parameter', 'quantities': ['efermi', 'forces'], 'link_name': 'my_custom_node'}
        """
        from copy import deepcopy

        self.nodes = {}

        # First, find all the nodes, that should be added.
        for key, value in self._settings.items():
            if not key.startswith('add_'):
                # only keys starting with 'add_' are relevant as nodes.
                continue
            if not value:
                # The quantity should not be added (This also excludes nodes with empty list or dict).
                continue

            node_name = key[4:]
            node_dict = deepcopy(NODES.get(node_name, {}))

            if isinstance(value, list):
                node_dict['quantities'] = value

            if isinstance(value, dict):
                node_dict.update(value)

            if 'link_name' not in node_dict:
                node_dict['link_name'] = node_name

            self.add_node(node_name, node_dict)

    def add_node(self, node_name, node_dict=None):
        """Add a definition of node to the nodes dictionary."""
        from copy import deepcopy

        if node_dict is None:
            # Try to get a node_dict from NODES.
            node_dict = deepcopy(NODES.get(node_name, {}))

        # Check, whether the node_dict contains required keys 'type' and 'quantities'
        for key in ['type', 'quantities']:
            if node_dict.get(key) is None:
                return

        self.nodes[node_name] = DictWithAttributes(node_dict)

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

    @property
    def quantities_to_parse(self):
        """Return the combined list of all the quantities, required for the current nodes."""
        quantities = []
        for value in self.nodes.values():
            for quantity in value['quantities']:
                if quantity in quantities:
                    continue
                quantities.append(quantity)
        return quantities
