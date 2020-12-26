"""
Node composer.

--------------
A composer that composes different quantities onto AiiDA data nodes.
"""
# pylint: disable=no-member, useless-object-inheritance, import-outside-toplevel
# Reason: pylint erroneously complains about non existing member 'get_quantity', which will be set in __init__.

from copy import deepcopy
from aiida_vasp.utils.aiida_utils import get_data_class
from aiida_vasp.parsers.quantity import ParsableQuantities
"""NODE_TYPES"""  # pylint: disable=pointless-string-statement

NODES_TYPES = {
    'dict': ['total_energies', 'maximum_force', 'maximum_stress', 'symmetries', 'magnetization', 'site_magnetization', 'notifications'],
    'array.kpoints': ['kpoints'],
    'structure': ['structure'],
    'array.trajectory': ['trajectory'],
    'array.bands': ['eigenvalues', 'kpoints', 'occupancies'],
    'vasp.chargedensity': ['chgcar'],
    'vasp.wavefun': ['wavecar'],
    'array': [],
}


class NodeComposer(object):
    """
    Prototype for a generic NodeComposer, that will compose output nodes based on parsed quantities.

    Provides methods to compose output_nodes from quantities. Currently supported node types are defined in NODES_TYPES.
    """

    def __init__(self, parsable_quantities=None, parsed_quantities=None, file_parsers=None):
        self._parsable_quantities = parsable_quantities
        self._parsed_quantities = parsed_quantities
        if file_parsers is not None:
            self._init_with_file_parsers(file_parsers)

    def _init_with_file_parsers(self, file_parsers):
        """Init with a list of file parsers."""
        self._parsable_quantities = ParsableQuantities()
        self._parsed_quantities = {}
        for parser in file_parsers:
            for key, value in parser.parsable_items.items():
                self._parsable_quantities.add_parsable_quantity(key, deepcopy(value))
                self._parsed_quantities[key] = parser.get_quantity(key)

    def compose(self, node_type, quantity_names=None):
        """
        A wrapper for compose_node with a node definition taken from NODES.

        :param node_type: str holding the type of the node. Must be one of the keys of NODES_TYPES.
        :param quantities: A list of strings with quantities to be used for composing this node.

        :return: An AiidaData object of a type corresponding to node_type.
        """

        if quantity_names is None:
            _quantity_names = NODES_TYPES.get(node_type)
        else:
            _quantity_names = quantity_names

        inputs = {}
        for quantity_name in _quantity_names:
            quantity = self._parsable_quantities.get_by_name(quantity_name)
            output_node = self._parsed_quantities.get(quantity_name)
            if output_node is None:
                for item in self._parsable_quantities.get_equivalent_quantities(quantity_name):
                    if item.original_name in self._parsed_quantities:
                        inputs[quantity.name] = self._parsed_quantities.get(item.original_name)
            else:
                inputs[quantity.name] = output_node

        # Call the correct specialised method for assembling.
        return getattr(self, '_compose_' + node_type.replace('.', '_'))(node_type, inputs)

    @staticmethod
    def _compose_dict(node_type, inputs):
        """Compose the dictionary node."""
        node = get_data_class(node_type)()
        node.update_dict(inputs)
        return node

    @staticmethod
    def _compose_structure(node_type, inputs):
        """Compose a structure node."""
        node = get_data_class(node_type)()
        for key in inputs:
            node.set_cell(inputs[key]['unitcell'])
            for site in inputs[key]['sites']:
                node.append_atom(position=site['position'], symbols=site['symbol'], name=site['kind_name'])
        return node

    @staticmethod
    def _compose_array(node_type, inputs):
        """Compose an array node."""
        node = get_data_class(node_type)()
        for item in inputs:
            for key, value in inputs[item].items():
                node.set_array(key, value)
        return node

    @staticmethod
    def _compose_vasp_wavefun(node_type, inputs):
        """Compose a wave function node."""
        node = None
        for key in inputs:
            # Technically this dictionary has only one key. to
            # avoid problems with python 2/3 it is done with the loop.
            node = get_data_class(node_type)(file=inputs[key])
        return node

    @staticmethod
    def _compose_vasp_chargedensity(node_type, inputs):
        """Compose a charge density node."""
        node = None
        for key in inputs:
            # Technically this dictionary has only one key. to
            # avoid problems with python 2/3 it is done with the loop.
            node = get_data_class(node_type)(file=inputs[key])
        return node

    def _compose_array_bands(self, node_type, inputs):
        """Compose a bands node."""
        node = get_data_class(node_type)()
        kpoints = self._compose_array_kpoints('array.kpoints', {'kpoints': inputs['kpoints']})
        node.set_kpointsdata(kpoints)
        node.set_bands(inputs['eigenvalues'], occupations=inputs['occupancies'])
        return node

    @staticmethod
    def _compose_array_kpoints(node_type, inputs):
        """Compose an array.kpoints node based on inputs."""
        node = get_data_class(node_type)()
        for key in inputs:
            mode = inputs[key]['mode']
            if mode == 'explicit':
                kpoints = inputs[key].get('points')
                cartesian = not kpoints[0].get_direct()
                kpoint_list = []
                weights = []
                for kpoint in kpoints:
                    kpoint_list.append(kpoint.get_point().tolist())
                    weights.append(kpoint.get_weight())

                if weights[0] is None:
                    weights = None

                node.set_kpoints(kpoint_list, weights=weights, cartesian=cartesian)

            if mode == 'automatic':
                mesh = inputs[key].get('divisions')
                shifts = inputs[key].get('shifts')
                node.set_kpoints_mesh(mesh, offset=shifts)
        return node

    @staticmethod
    def _compose_array_trajectory(node_type, inputs):
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
