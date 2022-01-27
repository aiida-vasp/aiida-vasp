"""
Node composer.

--------------
A composer that composes different quantities onto AiiDA data nodes.
"""
# pylint: disable=logging-fstring-interpolation
import traceback

from warnings import warn
import math
import numbers
import numpy as np

from aiida_vasp.utils.aiida_utils import get_data_class

NODES_TYPES = {
    'dict': [
        'total_energies', 'maximum_force', 'maximum_stress', 'symmetries', 'magnetization', 'site_magnetization', 'notifications',
        'band_properties', 'run_status', 'run_stats', 'version'
    ],
    'array.kpoints': ['kpoints'],
    'structure': ['structure'],
    'array.trajectory': ['trajectory'],
    'array.bands': ['eigenvalues', 'kpoints', 'occupancies'],
    'vasp.wavefun': ['wavecar'],
    'array': [],
}


class NodeComposer:
    """
    Prototype for a generic NodeComposer, that will compose output nodes based on parsed quantities.

    Provides methods to compose output_nodes from quantities. Currently supported node types are defined in NODES_TYPES.

    Parameters
    ----------
    """

    def __init__(self, nodes, equivalent_quantity_keys, quantities, logger=None):
        """Initialize."""
        self._equivalent_quantity_keys = equivalent_quantity_keys
        self._quantities = quantities
        self._nodes = nodes
        self._failed_to_create = []
        self._created = {}

        # Set logger
        if logger is not None:
            self._logger = logger
        else:
            import logging  # pylint: disable=import-outside-toplevel
            logging.basicConfig(level=logging.DEBUG)
            self._logger = logging.getLogger('NodeComposer')

        # Compose the nodes
        self.compose_nodes()

    def compose_nodes(self):
        """Compose the nodes according to the specifications."""
        for node_name, node_dict in self._nodes.items():
            inputs = {}
            for key in self._quantities.keys():
                # Make sure we strip prefixes as the quantities can contain
                # multiple equivalent keys, relevant only up to now.
                new_key = key
                if '-' in key:
                    new_key = key.split('-')[1]
                if key in node_dict['quantities'] or new_key in node_dict['quantities']:
                    # The found quantity should sit on this particular node.
                    # The delegated compose methods for each structure does not recognize prefixed keys,
                    # so use the stripped one.
                    inputs[new_key] = self._quantities[key]

            # If the input is empty, we skip creating the node as it is bound to fail
            if not inputs:
                self._failed_to_create.append(node_name)
                self._logger.warning(f'Creating node {node_dict["link_name"]} of type {node_dict["type"]} failed. '
                                     'No parsed data available.')
                continue

            exception = None  # pylint: disable=unused-variable
            # Guard the parsing in case of errors
            try:
                node = self.compose_node(node_dict['type'], inputs)
            except Exception:  # pylint: disable=broad-except
                node = None
                exception = traceback.format_exc()

            if node is not None:
                self._created[node_dict['link_name']] = node
            else:
                self._logger.warning(f'Creating node {node_dict["link_name"]} of type {node_dict["type"]} failed, '
                                     'exception: {exception}')
                self._failed_to_create.append(node_dict['link_name'])

    def compose_node(self, node_type, inputs):
        """
        A wrapper for compose_node with a node definition taken from NODES.

        :param node_type: str holding the type of the node. Must be one of the keys of NODES_TYPES.
        :param quantities: A list of strings with quantities to be used for composing this node.

        :return: An AiidaData object of a type corresponding to node_type.
        """

        if node_type in ('float', 'int', 'str'):
            return self._compose_basic_type(node_type, inputs)

        # Call the correct specialised method for assembling.
        method_name = 'compose_' + node_type.replace('.', '_')
        return getattr(self, method_name)(node_type, inputs)

    def compose_array_bands(self, node_type, inputs):
        """Compose a bands node."""
        typ = 'eigenvalues'
        if typ not in inputs:
            raise ValueError(f'The {typ} are not present after parsing.')
        node = get_data_class(node_type)()
        kpoints = self.compose_array_kpoints('array.kpoints', {'kpoints': inputs['kpoints']})
        node.set_kpointsdata(kpoints)
        if 'total' in inputs['eigenvalues']:
            eigenvalues = np.array([inputs['eigenvalues']['total']])
            occupancies = np.array([inputs['occupancies']['total']])
        else:
            eigenvalues = np.array([inputs['eigenvalues']['up'], inputs['eigenvalues']['down']])
            occupancies = np.array([inputs['occupancies']['up'], inputs['occupancies']['down']])
        node.set_bands(eigenvalues, occupations=occupancies)
        return node

    @staticmethod
    def compose_dict(node_type, inputs):
        """Compose the dictionary node."""
        node = get_data_class(node_type)()
        inputs = clean_nan_values(inputs)
        node.update_dict(inputs)
        return node

    @staticmethod
    def compose_structure(node_type, inputs):
        """
        Compose a structure node from consumable inputs.

        :param node_type: str Contains which type of AiiDA data structure you want for the node, i.e. ``structure``.
        :param inputs: Containing a key ``structure`` with keys ``sites`` and ``unitcell``.
        Here `sites` should contains an entry for each site with keys ``position``, ``symbol`` and ``kind_name``,
        which contain the position, the atomic symbol (e.g. ``Fe``) and the name to separate e.g.
        different symbols from each other, e.g. ``Fe1``.

        """
        typ = 'structure'
        if typ not in inputs:
            raise ValueError(f'The {typ} is not present after parsing.')
        node = get_data_class(node_type)()
        node.set_cell(inputs['structure']['unitcell'])
        for site in inputs['structure']['sites']:
            node.append_atom(position=site['position'], symbols=site['symbol'], name=site['kind_name'])
        return node

    @staticmethod
    def compose_array(node_type, inputs):
        """Compose an array node."""
        node = get_data_class(node_type)()
        for item in inputs:
            for key, value in inputs[item].items():
                node.set_array(key, value)
        return node

    @staticmethod
    def _compose_basic_type(node_type, inputs):
        """Compose a basic type node (int, float, str)."""
        node = None
        for key in inputs:
            # Technically this dictionary has only one key. to
            # avoid problems with python 2/3 it is done with the loop.
            node = get_data_class(node_type)(inputs[key])
        return node

    @staticmethod
    def _compose_vasp_wavefun(node_type, inputs):
        """Compose a wave function node."""
        node = None
        for key in inputs:
            # Technically this dictionary has only one key. to
            # avoid problems with python 2/3 it is done with the loop.
            node = get_data_class(node_type)(inputs[key])
        return node

    @staticmethod
    def compose_array_kpoints(node_type, inputs):
        """Compose an array.kpoints node based on inputs."""
        typ = 'kpoints'
        if typ not in inputs:
            raise ValueError(f'The {typ} are not present after parsing.')
        node = get_data_class(node_type)()
        for key in inputs:
            mode = inputs[key]['mode']
            if mode == 'explicit':
                kpoints = inputs[key].get('points')
                weights = inputs[key].get('weights')
                cartesian = inputs[key].get('cartesian')
                if weights[0] is None:
                    weights = None
                node.set_kpoints(kpoints, weights=weights, cartesian=cartesian)
            if mode == 'automatic':
                mesh = inputs[key].get('divisions')
                shifts = inputs[key].get('shifts')
                node.set_kpoints_mesh(mesh, offset=shifts)
        return node

    @staticmethod
    def compose_array_trajectory(node_type, inputs):
        """
        Compose a trajectory node.

        Parameters
        ----------
        node_type : str
            'array.trajectory'
        inputs : dict
            trajectory data is stored at VasprunParser. The keys are
            'cells', 'positions', 'symbols', 'forces', 'stress', 'steps'.

        Returns
        -------
        node : TrajectoryData
            To store the data, due to the definition of TrajectoryData in
            aiida-core v1.0.0, data, using the same keys are those from inputs,
            for 'symbols', the value is stored by set_attribute and
            for the others, the values are stored by by set_array.

        """
        node = get_data_class(node_type)()
        for item in inputs:
            for key, value in inputs[item].items():
                if key == 'symbols':
                    node.set_attribute(key, value)
                else:
                    node.set_array(key, value)
        return node

    @property
    def failed(self):
        return self._failed_to_create

    @property
    def successful(self):
        return self._created


def clean_nan_values(inputs: dict) -> dict:
    """
    Recursively replace NaN, Inf values (np.float) into None in place.

    This is because AiiDA does not support serializing these values
    as node attributes.
    """
    for key, value in inputs.items():
        if isinstance(value, dict):
            clean_nan_values(value)
        if isinstance(value, numbers.Real) and (math.isnan(value) or math.isinf(value)):
            warn('Key <{}> has value <{}> replaced by <{}>'.format(key, value, str(value)))
            inputs[key] = str(value)
    return inputs
