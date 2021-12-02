"""
Node composer.

--------------
A composer that composes different quantities onto AiiDA data nodes.
"""

from warnings import warn
import math
import numbers

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
    'vasp.chargedensity': ['chgcar'],
    'vasp.wavefun': ['wavecar'],
    'array': [],
}


def get_node_composer_inputs(equivalent_quantity_keys, parsed_quantities, quantity_names_in_node_dict):
    """
    Collect parsed quantities for the NodeCompoer input.

    When multiple equivalent quantities are found, the first one found in the
    equivalent_quantity_keys is chosen.

    """
    inputs = {}
    for quantity_name in quantity_names_in_node_dict:
        if quantity_name in equivalent_quantity_keys:
            for quantity_key in equivalent_quantity_keys[quantity_name]:
                if quantity_key in parsed_quantities:
                    inputs[quantity_name] = parsed_quantities[quantity_key]
                    break
    return inputs


def get_node_composer_inputs_from_file_parser(file_parser, quantity_keys=None):  # pylint: disable=invalid-name
    """Assemble necessary data from file_parser"""
    inputs = {}
    for key, value in file_parser.parsable_items.items():
        if quantity_keys is not None:
            if key not in quantity_keys:
                continue
        inputs[value['name']] = file_parser.get_quantity(key)
    return inputs


class NodeComposer:
    """
    Prototype for a generic NodeComposer, that will compose output nodes based on parsed quantities.

    Provides methods to compose output_nodes from quantities. Currently supported node types are defined in NODES_TYPES.
    """

    @classmethod
    def compose(cls, node_type, inputs):
        """
        A wrapper for compose_node with a node definition taken from NODES.

        :param node_type: str holding the type of the node. Must be one of the keys of NODES_TYPES.
        :param quantities: A list of strings with quantities to be used for composing this node.

        :return: An AiidaData object of a type corresponding to node_type.
        """

        if node_type in ('float', 'int', 'str'):
            return cls._compose_basic_type(node_type, inputs)

        # Call the correct specialised method for assembling.
        method_name = '_compose_' + node_type.replace('.', '_')
        return getattr(cls, method_name)(node_type, inputs)

    @staticmethod
    def _compose_dict(node_type, inputs):
        """Compose the dictionary node."""
        node = get_data_class(node_type)()
        inputs = clean_nan_values(inputs)
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

    @classmethod
    def _compose_array_bands(cls, node_type, inputs):
        """Compose a bands node."""
        node = get_data_class(node_type)()
        kpoints = cls._compose_array_kpoints('array.kpoints', {'kpoints': inputs['kpoints']})
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
