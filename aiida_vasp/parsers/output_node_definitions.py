"""Definitions for output_nodes and a NodeAssembler."""
from aiida_vasp.utils.aiida_utils import get_data_class

NODES = {
    'parameters': {
        'link_name': 'output_parameters',
        'type': 'parameter',
        'quantities': ['outcar-volume', 'outcar-energies', 'outcar-fermi_level', 'symmetries', 'parameters'],
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
    'trajectory': {
        'link_name': 'output_trajectory',
        'type': 'array.trajectory',
        'quantities': ['trajectory'],
    },
    'bands': {
        'link_name': 'output_bands',
        'type': 'array.bands',
        'quantities': ['bands'],
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
        'link_name': 'chgcar',
        'type': 'vasp.chargedensity',
        'quantities': ['chgcar'],
    },
    'wavecar': {
        'link_name': 'wavecar',
        'type': 'vasp.wavefun',
        'quantities': ['wavecar'],
    },
}


class NodeAssembler(object):
    """
    Prototype for a generic NodeAssembler, that will assemble output nodes based on parsed quantities.

    Provides methods to assemble output_nodes from quantities. Currently supported node types are:

        * StructureData,
        * vasp.wavefun,
        * vasp.chargedensity
        * array,
        * parameter,

    missing types are:

        * array.kpoints
        * array.bands
    """

    def __init__(self):
        self._current_node = None

    def assemble(self, node, quantities):
        """
        Assemble an output_node from a set of quantities.

        :param node: AttributeDict with the nodes definition.
        :param quantities: List of quantities required for assembling this node.
        :return: An AiidaData object.
        """

        self._current_node = get_data_class(node.type)()
        inputs = {}
        for quantity in node.quantities:
            inputs[quantity] = quantities.get(quantity)

        # Call the correct specialised method for assembling.
        getattr(self, "_assemble_with_" + node.type.replace(".", "_"))(inputs)

        return self._current_node

    def _assemble_with_parameter(self, inputs):
        self._current_node.update_dict(inputs)

    def _assemble_with_structure(self, inputs):
        self._current_node.set_cell(inputs['structure']['unitcell'])
        for site in inputs['structure']['sites']:
            self._current_node.append_atom(position=site['position'], symbols=site['symbol'], name=site['kind_name'])

    def _assemble_with_vasp_wavefun(self, inputs):
        self._current_node.set_file(inputs['wavecar'])

    def _assemble_with_vasp_chargedensity(self, inputs):
        self._current_node.set_file(inputs['chgcar'])
