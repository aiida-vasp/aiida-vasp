"""A Node Composer for Aiida data nodes."""
# pylint: disable=no-member
# Reason: pylint erroneously complains about non existing member 'get_quantity', which will be set in __init__.

from aiida_vasp.utils.aiida_utils import get_data_class
from aiida_vasp.utils.delegates import delegate_method_kwargs, Delegate
from aiida_vasp.parsers.quantity import ParsableQuantities

NODES_TYPES = {
    'parameter': ['total_energies', 'maximum_forces', 'maximum_stress', 'fermi_level'],
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

    Provides methods to compose output_nodes from quantities. Currently supported node types are defined in NODE_TYPES.
    """

    def __init__(self, **kwargs):
        # create the delegate for getting quantities.
        setattr(self, 'get_quantity', Delegate())
        self.quantites = None
        self.init_with_kwargs(**kwargs)

    @delegate_method_kwargs(prefix='_init_with_')
    def init_with_kwargs(self, **kwargs):
        """Delegate initialization to _init_with - methods."""

    def _init_with_file_parsers(self, file_parsers):
        """Init with a list of file parsers."""
        from copy import deepcopy

        if not file_parsers:
            return

        self.quantites = ParsableQuantities()

        # Add all the FileParsers get_quantity methods to our get_quantity delegate.
        for parser in file_parsers:
            for key, value in parser.parsable_items.items():
                self.quantites.add_parsable_quantity(key, deepcopy(value))

            self.get_quantity.append(parser.get_quantity)

    def _init_with_vasp_parser(self, vasp_parser):
        """Init with a VaspParser object."""
        self.get_quantity.append(vasp_parser.get_inputs)
        self.quantites = vasp_parser.quantities

    def compose(self, node_type, quantities=None):
        """
        A wrapper for compose_node with a node definition taken from NODES.

        :param node_type: str holding the type of the node. Must be one of the keys of NODES_TYPES.
        :param quantities: A list of strings with quantities to be used for composing this node.

        :return: An AiidaData object of a type corresponding to node_type.
        """

        if quantities is None:
            quantities = NODES_TYPES.get(node_type)

        current_node = get_data_class(node_type)()

        inputs = {}
        for quantity_name in quantities:
            quantity = self.quantites.get_by_name(quantity_name)
            inputs[quantity.name] = self.get_quantity(quantity_name)[quantity_name]

        # Call the correct specialised method for assembling.
        getattr(self, "_compose_with_" + node_type.replace(".", "_"))(current_node, inputs)

        return current_node

    @staticmethod
    def _compose_with_parameter(node, inputs):
        node.update_dict(inputs)

    @staticmethod
    def _compose_with_structure(node, inputs):
        for key in inputs:
            node.set_cell(inputs[key]['unitcell'])
            for site in inputs[key]['sites']:
                node.append_atom(position=site['position'], symbols=site['symbol'], name=site['kind_name'])

    @staticmethod
    def _compose_with_array(node, inputs):
        for item in inputs:
            for key, value in inputs[item].items():
                node.set_array(key, value)

    @staticmethod
    def _compose_with_vasp_wavefun(node, inputs):
        for key in inputs:
            node.set_file(inputs[key])

    @staticmethod
    def _compose_with_vasp_chargedensity(node, inputs):
        for key in inputs:
            node.set_file(inputs[key])

    def _compose_with_array_bands(self, node, inputs):

        kpoints = get_data_class('array.kpoints')()
        self._compose_with_array_kpoints(kpoints, {'kpoints': inputs['kpoints']})
        node.set_kpointsdata(kpoints)
        node.set_bands(inputs['eigenvalues'], occupations=inputs['occupancies'])

    @staticmethod
    def _compose_with_array_kpoints(node, inputs):
        """Compose an array.kpoints node based on inputs."""
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

    @staticmethod
    def _compose_with_array_trajectory(node, inputs):
        for item in inputs:
            for key, value in inputs[item].items():
                node.set_array(key, value)
