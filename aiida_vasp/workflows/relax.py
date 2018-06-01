"""Structure relaxation workchain for VASP."""
from aiida.common.extendeddicts import AttributeDict
from aiida.orm.data.base import Float, Int, Bool
from aiida.work import WorkChain

from aiida_vasp.utils.aiida_utils import get_data_class, get_data_node
from aiida_vasp.workflows.base import VaspBaseWf


RELAXATION_INCAR_TEMPLATE = AttributeDict({
    "istart": 0,
    "ismear": -1,
    "sigma": 0.2,
    "algo": 'NORMAL',
    "ediff": 1E-6,
    "prec": 'H',
    "ibrion": 2,
    "nsw": 200,
    "melmin": 4,
    "ediffg": -0.01,
    "isif": 4
})


def relax_parameters(ions, cell, shape):
    """Set IBRION and ISIF according to the choice of what should be relaxed."""


class VaspRelaxWf(VaspBaseWf):
    """Structure relaxation workchain for VASP."""

    @classmethod
    def(cls, spec):
        spec.expose_inputs(VaspBaseWf, include=['code', 'structure', 'potcar_family', 'potcar_mapping'])
        spec.input_namespace('kpoints')
        spec.input('mesh', valid_type=get_data_class('array.kpoints'), namespace='kpoints', required=False)
        spec.input('distance', valid_type=Float, namespace='kpoints', required=False)
        spec.input('incar_add', valid_type=get_data_class('parameter'), required=False)
        spec.input_namespace('relax')
        spec.input('ions', valid_type=Bool, required=False)
        spec.input('shape', valid_type=Bool, required=False)
        spec.input('volume', valid_type=Bool, required=False)
        spec.input('options', valid_type=get_data_class('parameter'), required=False)
        spec.expose_inputs(VaspBaseWf, namespace='restart', include=['max_iterations'])

        spec.outline(
            cls.setup,
            cls.validate_inputs,
            cls.run_relax,
            cls.results
        )

        spec.output('output_structure', valid_type=get_data_class('structure'))
        spec.output('output_parameters', valid_type=get_data_class('parameter'))
        spec.output('remote_folder', valid_type=get_data_class('remote'))
        spec.output('retrieved', valid_type=get_data_class('folder'))

    def setup(self):
        self.ctx.inputs = AttributeDict()
        self.ctx.inputs.structure = self.inputs.structure
        self.ctx.inputs.potcar_family = self.inputs.potcar_family
        self.ctx.inputs.potcar_mapping = self.inputs.potcar_mapping

        if 'mesh' in self.inputs.kpoints:
            self.ctx.inputs.kpoints = self.inputs.kpoints.mesh
        elif 'distance' in self.inputs.kpoints:
            self.ctx.inputs.kpoints = get_data_node(
                'array.kpoints', cell_from_structure=self.inputs.structure, kpoints_mesh_from_density=self.inputs.kpionts.distance.value)
        else:
            self._fail_compat(exception=ValueError('either kpoints.mesh or kpoints.distance is required'))

        self.ctx.inputs.incar = RELAXATION_INCAR_TEMPLATE
        if 'incar_add' in self.inputs:
            incar_additions = self.inputs.incar_add.get_dict()
            incar_add_keys = [key.lower() for key in incar_additions.keys()]
            if 'ibrion' in incar_add_keys:
                self._fail_compat(exception=ValueError('overriding IBRION not allowed, use relax.ions to control'))
            incar = self.ctx.inputs.incar_add
