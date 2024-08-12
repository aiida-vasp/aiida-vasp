"""
Parser settings.

----------------
Defines the object associated with each object parser and the default node
compositions of quantities.
"""

# pylint: disable=import-outside-toplevel
from copy import deepcopy
from dataclasses import dataclass
from typing import List

from pydantic import BaseModel

from aiida_vasp.parsers.content_parsers import *
from aiida_vasp.utils.extended_dicts import update_nested_dict


@dataclass
class ContentParserConfig:
    """Content parser set configuration"""

    parser_class: BaseFileParser
    is_critical: bool
    mode: str
    status: str = 'Unknown'


CONTENT_PARSER_SETS = {
    'default': {
        'DOSCAR': ContentParserConfig(DoscarParser, False, 'r'),
        'EIGENVAL': ContentParserConfig(EigenvalParser, False, 'r'),
        'IBZKPT': ContentParserConfig(KpointsParser, False, 'r'),
        'OUTCAR': ContentParserConfig(OutcarParser, True, 'r'),
        'vasprun.xml': ContentParserConfig(VasprunParser, True, 'rb'),
        'CHGCAR': ContentParserConfig(ChgcarParser, False, 'r'),
        'CONTCAR': ContentParserConfig(PoscarParser, False, 'r'),
        'vasp_output': ContentParserConfig(StreamParser, False, 'r'),
    },
    'neb': {
        'DOSCAR': ContentParserConfig(DoscarParser, False, 'r'),
        'EIGENVAL': ContentParserConfig(EigenvalParser, False, 'r'),
        'IBZKPT': ContentParserConfig(KpointsParser, False, 'r'),
        'OUTCAR': ContentParserConfig(VtstNebOutcarParser, False, 'r'),
        'vasprun.xml': ContentParserConfig(VasprunParser, False, 'rb'),
        'CHGCAR': ContentParserConfig(ChgcarParser, False, 'r'),
        'CONTCAR': ContentParserConfig(PoscarParser, False, 'r'),
        # The STDOUT is rename as 'stdout' for NEB calculations, this is because VASP itself
        # divert STDOUT for each image to <ID>/stdout
        'stdout': ContentParserConfig(StreamParser, False, 'r'),
    },
}


class NodeConfig(BaseModel):
    """Output node configuration"""

    link_name: str
    tyle: str
    quantities: List[str]


class NodeConfig(BaseModel):
    """Configurations for the output nodes"""

    misc: NodeConfig = NodeConfig(
        link_name='misc',
        tyle='core.dict',
        quantities=[
            'total_energies',
            'maximum_stress',
            'maximum_force',
            'notifications',
            'run_status',
            'run_stats',
            'version',
        ],
    )
    kpoints: NodeConfig = NodeConfig(link_name='kpoints', tyle='core.array.kpoints', quantities=['kpoints'])
    structure: NodeConfig = NodeConfig(link_name='structure', tyle='core.structure', quantities=['structure'])
    poscar_structure: NodeConfig = NodeConfig(
        link_name='structure', tyle='core.structure', quantities=['poscar-structure']
    )
    trajectory: NodeConfig = NodeConfig(link_name='trajectory', tyle='core.array.trajectory', quantities=['trajectory'])
    forces: NodeConfig = NodeConfig(link_name='forces', tyle='core.array', quantities=['forces'])
    stress: NodeConfig = NodeConfig(link_name='stress', tyle='core.array', quantities=['stress'])
    bands: NodeConfig = NodeConfig(
        link_name='bands', tyle='core.array.bands', quantities=['eigenvalues', 'kpoints', 'occupancies']
    )
    dos: NodeConfig = NodeConfig(link_name='dos', tyle='core.array', quantities=['dos'])
    energies: NodeConfig = NodeConfig(link_name='energies', tyle='core.array', quantities=['energies'])
    projectors: NodeConfig = NodeConfig(link_name='projectors', tyle='core.array', quantities=['projectors'])
    born_charges: NodeConfig = NodeConfig(link_name='born_charges', tyle='core.array', quantities=['born_charges'])
    dielectrics: NodeConfig = NodeConfig(link_name='dielectrics', tyle='core.array', quantities=['dielectrics'])
    hessian: NodeConfig = NodeConfig(link_name='hessian', tyle='core.array', quantities=['hessian'])
    dynmat: NodeConfig = NodeConfig(link_name='dynmat', tyle='core.array', quantities=['dynmat'])
    charge_density: NodeConfig = NodeConfig(
        link_name='charge_density', tyle='core.array', quantities=['charge_density']
    )
    wavecar: NodeConfig = NodeConfig(link_name='wavecar', tyle='vasp.wavefun', quantities=['wavecar'])
    site_magnetization: NodeConfig = NodeConfig(
        link_name='site_magnetization', tyle='core.dict', quantities=['site_magnetization']
    )


""" NODES """

NODES = {
    'misc': {
        'link_name': 'misc',
        'type': 'core.dict',
        'quantities': [
            'total_energies',
            'maximum_stress',
            'maximum_force',
            'notifications',
            'run_status',
            'run_stats',
            'version',
        ],
    },
    'kpoints': {
        'link_name': 'kpoints',
        'type': 'core.array.kpoints',
        'quantities': ['kpoints'],
    },
    'structure': {
        'link_name': 'structure',
        'type': 'core.structure',
        'quantities': ['structure'],
    },
    'poscar-structure': {
        'link_name': 'structure',
        'type': 'core.structure',
        'quantities': ['poscar-structure'],
    },
    'trajectory': {
        'link_name': 'trajectory',
        'type': 'core.array.trajectory',
        'quantities': ['trajectory'],
    },
    'forces': {
        'link_name': 'forces',
        'type': 'core.array',
        'quantities': ['forces'],
    },
    'stress': {
        'link_name': 'stress',
        'type': 'core.array',
        'quantities': ['stress'],
    },
    'bands': {
        'link_name': 'bands',
        'type': 'core.array.bands',
        'quantities': ['eigenvalues', 'kpoints', 'occupancies'],
    },
    'dos': {
        'link_name': 'dos',
        'type': 'core.array',
        'quantities': ['dos'],
    },
    'energies': {
        'link_name': 'energies',
        'type': 'core.array',
        'quantities': ['energies'],
    },
    'projectors': {
        'link_name': 'projectors',
        'type': 'core.array',
        'quantities': ['projectors'],
    },
    'born_charges': {
        'link_name': 'born_charges',
        'type': 'core.array',
        'quantities': ['born_charges'],
    },
    'dielectrics': {
        'link_name': 'dielectrics',
        'type': 'core.array',
        'quantities': ['dielectrics'],
    },
    'hessian': {
        'link_name': 'hessian',
        'type': 'core.array',
        'quantities': ['hessian'],
    },
    'dynmat': {
        'link_name': 'dynmat',
        'type': 'core.array',
        'quantities': ['dynmat'],
    },
    'charge_density': {
        'link_name': 'charge_density',
        'type': 'core.array',
        'quantities': ['charge_density'],
    },
    'wavecar': {
        'link_name': 'wavecar',
        'type': 'vasp.wavefun',
        'quantities': ['wavecar'],
    },
    'site_magnetization': {
        'link_name': 'site_magnetization',
        'type': 'core.dict',
        'quantities': ['site_magnetization'],
    },
}


class ParserDefinitions:  # pylint: disable=useless-object-inheritance
    """Container of parser definitions"""

    def __init__(self, content_parser_set='default'):
        self._parser_definitions = {}
        if content_parser_set not in CONTENT_PARSER_SETS:
            return
        for name, parser_config in CONTENT_PARSER_SETS.get(content_parser_set).items():
            self._parser_definitions[name] = deepcopy(parser_config)

    @property
    def parser_definitions(self):
        return self._parser_definitions


class ParserSettings:  # pylint: disable=useless-object-inheritance
    """
    Settings object for the VaspParser.

    :param settings: Dict with the 'parser_settings'.
    :param default_settings: Dict with default settings.

    This provides the following properties to other components of the VaspParser:

        * nodes_dict: A list with all requested output nodes.
        * parser_definitions: A Dict with the definitions for each specific content parser.
        * quantities_to_parse: Collection of quantities in nodes_dict.

    """

    NODES = NODES

    def __init__(self, settings, default_settings=None, vasp_parser_logger=None):  # pylint: disable=missing-function-docstring
        self._vasp_parser_logger = vasp_parser_logger
        if settings is None:
            self._settings = {}
        else:
            self._settings = settings

        # If the default is supplied use it as the base and update with the explicity settings
        if default_settings is not None:
            new_settings = deepcopy(default_settings)
            update_nested_dict(new_settings, self._settings)
            self._settings = new_settings

        self._output_nodes_dict = {}
        self._critical_error_list = []
        self._init_output_nodes_dict()
        self._init_critical_error_list()

    @property
    def output_nodes_dict(self):
        return self._output_nodes_dict

    @property
    def quantity_names_to_parse(self):
        """Return the combined list of all the quantities, required for the current nodes."""
        quantity_names_to_parse = []
        for node_dict in self.output_nodes_dict.values():
            for quantity_key in node_dict['quantities']:
                if quantity_key in quantity_names_to_parse:
                    continue
                quantity_names_to_parse.append(quantity_key)
        return quantity_names_to_parse

    @property
    def critical_notifications_to_check(self):
        """Return the list of critical notification names to be checked"""
        return self._critical_error_list

    def add_output_node(self, node_name, node_dict, is_custom_node=False):
        """Add a definition of node to the nodes dictionary."""
        _node_dict = deepcopy(node_dict)

        if is_custom_node:
            if 'link_name' not in _node_dict:
                _node_dict['link_name'] = node_name
                self._vasp_parser_logger.info(f"'{node_name}' was set as 'link_name' in node_dict.")
            if 'quantities' not in node_dict:
                _node_dict['quantities'] = [node_name]
                self._vasp_parser_logger.info(f"'{node_name}' was set as 'quantities' in node_dict.")

        # Check, whether the node_dict contains required keys.
        exist_missing_key = False
        for key in ['type', 'quantities', 'link_name']:
            if key not in _node_dict:
                self._vasp_parser_logger.warning(f"'{key}' was not found in node_dict.")
                exist_missing_key = True
        if not exist_missing_key:
            self._output_nodes_dict[node_name] = _node_dict

    @property
    def settings(self):
        """Return the settings dictionary."""
        return self._settings

    def get(self, item, default=None):
        return self._settings.get(item, default)

    def update_quantities_to_parse(self, new_quantities):
        """Update the quantities to be parsed."""
        try:
            # Update the quantities to be parsed, any extra keys already sitting in settings are preserved
            self._settings['quantities_to_parse'].append(new_quantities)
            self._settings['quantities_to_parse'] = list(dict.fromkeys(self._settings['quantities_to_parse']))
        except KeyError:
            self._settings['quantities_to_parse'] = new_quantities

    def _init_output_nodes_dict(self):
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

                'add_custom_node': {'type': 'parameter', 'quantities': ['efermi', 'forces'], 'link_name':
                                    'my_custom_node'}
        """
        self._output_nodes_dict = {}

        # First, find all the nodes, that should be added.
        for key, value in self._settings.items():
            if not key.startswith('add_'):
                # only keys starting with 'add_' are relevant as nodes.
                continue
            if not value:
                # The quantity should not be added (This also excludes nodes with empty list or dict).
                continue

            node_name = key[4:]
            node_dict = deepcopy(self.NODES.get(node_name, {}))
            # Considered as custom node if node_dict == {}.
            is_custom_node = not bool(node_dict)

            if isinstance(value, list):
                node_dict['quantities'] = value

            if isinstance(value, dict):
                node_dict.update(value)

            self.add_output_node(node_name, node_dict, is_custom_node=is_custom_node)

    def _init_critical_error_list(self):
        """
        Set the critical error list to be checked from a settings object.

        Name of critical notifications can be added by setting:

            'add_name' : True

        in ``parser_settings``'s ``critical_notifications`` field.
        The the notifications can be removed by setting:

            'add_name' : False

        from the default set defined by CRITICAL_NOTIFICATIONS.
        """
        self._critical_error_list = []

        # First, find all the nodes, that should be added.
        for key, value in self._settings.get('critical_notifications', {}).items():
            if key.startswith('add_'):
                # only keys starting with 'add_' are relevant as nodes.
                if value:
                    self._critical_error_list.append(key[4:])
