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
    "prec": 'Normal',
    "ibrion": 2,
    "nsw": 200,
    "melmin": 4,
    "ediffg": -0.01,
    "isif": 4
})


def relax_parameters(pos, shape, volume):
    """Set IBRION and ISIF according to the choice of what should be relaxed."""
    if not (pos or volume or shape):
        raise ValueError('relaxation workflow requires at least one of "positions", "cell shape", "cell volume" to be alowed to change.')
    params = AttributeDict({'ibrion': 2, 'isif': 4})
    if not pos:
        params.ibrion = 0
        if shape and volume:
            params.isif = 6
        elif not volume:
            params.isif = 5
        elif not shape:
            params.isif = 7
    else:
        ibrion = 2
        if shape and volume:
            params.isif = 3
        elif not volume:
            params.isif = 4
        else:
            raise ValueError('impossible to relax positions and cell volume without relaxing cell shape')
    return params


def encut_factor(volume):
    if volume:
        return 1.3
    return 1


def clean_incar_overrides(inputs):
    overrides = AttributeDict(inputs.incar_add.get_dict())
    if 'ibrion' in overrides:
        raise ValueError('overriding IBRION not allowed, use relax.xxx inputs to control')
    if 'isif' in overrides:
        raise ValueError('overriding ISIF not allowed, use relax.xxx inputs to control')
    if 'nsw' in overrides:
        if inputs.relax.positions.value and overrides.nsw < 1:
            raise ValueError('NSW (num ionic steps) was set to 0 but relaxing positions was requested')

    return overrides


class VaspRelaxWf(VaspBaseWf):
    """Structure relaxation workchain for VASP."""

    @classmethod
    def define(cls, spec):
        spec.expose_inputs(VaspBaseWf, include=['code', 'structure', 'potcar_family', 'potcar_mapping'])
        spec.input_namespace('kpoints')
        spec.input('mesh', valid_type=get_data_class('array.kpoints'), namespace='kpoints', required=False)
        spec.input('distance', valid_type=Float, namespace='kpoints', required=False)
        spec.input('incar_add', valid_type=get_data_class('parameter'), required=False)
        spec.input_namespace('relax')
        spec.input('positions', valid_type=Bool, required=False)
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

    def _clean_incar(self, base_incar):
        """Update incar parameters based on other inputs."""
        incar = base_incar.copy()
        try:
            incar.update(relax_parameters(
                pos=self.inputs.relax.positions.value, shape=self.inputs.relax.shape.value, volume=self.inputs.relax.shape.value))
        except ValueError as err:
            self._fail_compat(exception=err)

        if self.inputs.relax.shape.value:
            potcars = MultiPotcarIo.from_structure(self.inputs.structure, get_data_class('potcar').get_potcars_from_structure(structure=self.inputs.structure, family_name=self.inputs.potcar_family.value, mapping=self.inputs.potcar_mapping.get_dict()))
            incar.encut = encut_factor(self.inputs.relax.volume) * potcars.max_enmax

        if 'incar_add' in self.inputs:
            try:
                incar_overrides = clean_incar_overrides(self.inputs)
                incar.update(incar_overrides)
            except ValueError as err:
                self._fail_compat(exception=err)
        return incar

    def _clean_kpoints(self):
        kpoints = None
        if 'mesh' in self.inputs.kpoints:
            kpoints = self.inputs.kpoints.mesh
        elif 'distance' in self.inputs.kpoints:
            kpoints = get_data_node(
                'array.kpoints', cell_from_structure=self.inputs.structure, kpoints_mesh_from_density=self.inputs.kpionts.distance.value)
        else:
            self._fail_compat(exception=ValueError('either kpoints.mesh or kpoints.distance is required'))
        return kpoints

    def setup(self):
        self.ctx.inputs = AttributeDict()
        self.ctx.inputs.structure = self.inputs.structure
        self.ctx.inputs.potcar_family = self.inputs.potcar_family
        self.ctx.inputs.potcar_mapping = self.inputs.potcar_mapping
        self.ctx.inputs.kpoints = self._clean_kpoints()
        self.ctx.inputs.incar = self._clean_incar(RELAXATION_INCAR_TEMPLATE)
