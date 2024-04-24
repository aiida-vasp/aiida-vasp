"""
Bands workchain with a more flexible input
TODO:

- Add option to use alternative pathways obtained using sumo-interface
- Improve the hybrid workchain by performing local dryrun to extract the full kpoints
    - If running SOC, the ISYM should be turned to 0 or -1.
"""
from copy import deepcopy
from typing import List

import numpy as np

from aiida.common.extendeddicts import AttributeDict
from aiida.common.links import LinkType
from aiida.engine import WorkChain, append_, calcfunction, if_
import aiida.orm as orm
from aiida.orm.nodes.data.base import to_aiida_type
from aiida.plugins import WorkflowFactory

from aiida_vasp.utils.aiida_utils import get_data_class

from .common import OVERRIDE_NAMESPACE, nested_update, nested_update_dict_node
from .common.transform import magnetic_structure_decorate, magnetic_structure_dedecorate
from .mixins import WithVaspInputSet

try:
    from aiida_vasp.parsers.content_parsers.vasprun import VasprunParser
except ImportError:
    from aiida_vasp.parsers.file_parsers.vasprun import VasprunParser


class VaspBandsWorkChain(WorkChain, WithVaspInputSet):
    """
    Workchain for running bands calculations.

    This workchain does the following:

    1. Relax the structure if requested (eg. inputs passed to the relax namespace).
    2. Do a SCF singlepoint calculation.
    3. Do a non-scf calculation for bands and dos.

    Inputs must be passed for the SCF calculation, others are optional. The dos calculation
    will only run if the kpoints for DOS are passed or a full specification is given under the
    `dos` input namesace.

    The SCF calculation may be skipped by passing a CHGCAR file/remote folder. In which case the SCF inputs
    are carried on for non-scf calculations.

    The band structure calculation will run unless `only_dos` is set to `Bool(True)`.

    For magnetic structures, the workchain will internally create additional species for the symmetry
    analysis and regenerate "undecorated" structures with corresponding initial magnetic moments. This
    works for both FM and AFM species. Care should be taken when the MAGMOM is obtained from site projected
    values in case of unexpected symmetry breaking.

    Input for bands and dos calculations are optional. However, if they are needed, the full list of inputs must
    be passed. For the `parameters` node, one may choose to only specify those fields that need to be updated.
    """

    _base_wk_string = 'vasp.v2.vasp'
    _relax_wk_string = 'vasp.v2.relax'
    DEFAULT_BAND_MODE = 'bradcrack'
    DEFAULT_LINE_DENSITY = 20
    DEFAULT_SYMPREC = 0.01

    @classmethod
    def define(cls, spec):
        """Initialise the WorkChain class"""
        super().define(spec)
        relax_work = WorkflowFactory(cls._relax_wk_string)
        base_work = WorkflowFactory(cls._base_wk_string)

        spec.input('structure', help='The input structure', valid_type=orm.StructureData)
        spec.input(
            'bs_kpoints',
            help='Explicit kpoints for the bands. Will not generate kpoints if supplied.',
            valid_type=orm.KpointsData,
            required=False,
        )
        spec.input(
            'bands_kpoints_distance',
            help='Spacing for band distances for automatic kpoints generation, used by seekpath-aiida mode.',
            valid_type=orm.Float,
            required=False,
        )
        spec.input(
            'line_density',
            help='Density of the point along the path, used by sumo interface.',
            valid_type=orm.Float,
            required=False,
        )
        spec.input(
            'band_mode',
            help=
            'Mode for generating the band path. Choose from: bradcrack, pymatgen, seekpath, seekpath-aiida and latimer-munro.',
            required=False,
            valid_type=orm.Str,
        )
        spec.input(
            'symprec',
            help='Precision of the symmetry determination',
            valid_type=orm.Float,
            required=True,
        )
        spec.input(
            'dos_kpoints_density',
            help='Kpoints for running DOS calculations in A^-1 * 2pi. Will perform non-SCF DOS calculation is supplied.',
            serializer=to_aiida_type,
            required=False,
            valid_type=orm.Float,
        )
        spec.expose_inputs(
            relax_work,
            namespace='relax',
            exclude=('structure',),
            namespace_options={
                'required': False,
                'populate_defaults': False,
                'help': 'Inputs for Relaxation workchain, if needed',
            },
        )
        spec.expose_inputs(
            base_work,
            namespace='scf',
            exclude=('structure',),
            namespace_options={
                'required': True,
                'populate_defaults': True,
                'help': 'Inputs for SCF workchain, mandatory',
            },
        )
        spec.expose_inputs(
            base_work,
            namespace='bands',
            exclude=('structure', 'kpoints'),
            namespace_options={
                'required': False,
                'populate_defaults': False,
                'help': 'Inputs for bands calculation, if needed',
            },
        )
        spec.expose_inputs(
            base_work,
            namespace='dos',
            exclude=('structure',),
            namespace_options={
                'required': False,
                'populate_defaults': False,
                'help': 'Inputs for DOS calculation, if needed',
            },
        )
        spec.input(
            'clean_children_workdir',
            valid_type=orm.Str,
            serializer=to_aiida_type,
            help='What part of the called children to clean',
            required=False,
            default=lambda: orm.Str('none'),
        )
        spec.input(
            'only_dos',
            required=False,
            help='Flag for running only DOS calculations',
            valid_type=orm.Bool,
            serializer=to_aiida_type,
        )
        spec.input(
            'chgcar',
            required=False,
            valid_type=get_data_class('vasp.chargedensity'),
            help='Explicit CHGCAR file used for DOS/Bands calculations',
        )
        spec.input(
            'restart_folder',
            required=False,
            valid_type=get_data_class('remote'),
            help='A remote folder containing the CHGCAR file to be used',
        )
        spec.outline(
            cls.setup,
            if_(cls.should_do_relax)(
                cls.run_relax,
                cls.verify_relax,
            ),
            if_(cls.should_generate_path)(cls.generate_path),
            if_(cls.should_run_scf)(
                cls.run_scf,
                cls.verify_scf,
            ),
            cls.run_bands_dos,
            cls.inspect_bands_dos,
        )

        spec.output(
            'primitive_structure',
            required=False,
            help='Primitive structure used for band structure calculations',
        )
        spec.output('band_structure', required=False, help='Computed band structure with labels')
        spec.output('seekpath_parameters', help='Parameters used by seekpath', required=False)
        spec.output('dos', required=False)
        spec.output('projectors', required=False)

        spec.exit_code(501, 'ERROR_SUB_PROC_RELAX_FAILED', message='Relaxation workchain failed')
        spec.exit_code(502, 'ERROR_SUB_PROC_SCF_FAILED', message='SCF workchain failed')
        spec.exit_code(
            503,
            'ERROR_SUB_PROC_BANDS_FAILED',
            message='Band structure workchain failed',
        )
        spec.exit_code(504, 'ERROR_SUB_PROC_DOS_FAILED', message='DOS workchain failed')
        spec.exit_code(
            601,
            'ERROR_INPUT_STRUCTURE_NOT_PRIMITIVE',
            message='The input structure is not the primitive one!',
        )

    def select_chgcar_from_inputs(self):
        """Setup CHGCAR from inputs"""
        if self.inputs.get('chgcar'):
            self.ctx.chgcar = self.inputs.chgcar
            self.report(f'Using CHGCAR {self.inputs.chgcar} from input')
        else:
            self.ctx.chgcar = None

        if self.inputs.get('restart_folder'):
            self.ctx.restart_folder = self.inputs.restart_folder
            self.report(f'Using remote folder {self.inputs.restart_folder} for restart')
        else:
            self.ctx.restart_folder = None

    def setup(self):
        """Setup the calculation"""
        self.ctx.current_structure = self.inputs.structure
        self.ctx.bs_kpoints = self.inputs.get('bs_kpoints')
        param = self.inputs.scf.parameters.get_dict()
        if 'magmom' in param[OVERRIDE_NAMESPACE] and not self.inputs.get('only_dos'):
            self.report('Magnetic system passed for BS')
            self.ctx.magmom = param[OVERRIDE_NAMESPACE]['magmom']
        else:
            self.ctx.magmom = None

    def should_do_relax(self):
        """Wether we should do relax or not"""
        return 'relax' in self.inputs

    def run_relax(self):
        """Run the relaxation"""
        relax_work = WorkflowFactory(self._relax_wk_string)
        inputs = self.exposed_inputs(relax_work, 'relax', agglomerate=True)
        inputs = AttributeDict(inputs)
        inputs.metadata.call_link_label = 'relax'
        inputs.structure = self.ctx.current_structure

        running = self.submit(relax_work, **inputs)
        return self.to_context(workchain_relax=running)

    def verify_relax(self):
        """Verify the relaxation"""
        relax_workchain = self.ctx.workchain_relax
        if not relax_workchain.is_finished_ok:
            self.report('Relaxation finished with Error')
            return self.exit_codes.ERROR_SUB_PROC_RELAX_FAILED

        # Use the relaxed structure as the current structure
        self.ctx.current_structure = relax_workchain.outputs.relax__structure

    def should_run_scf(self):
        """Wether we should run SCF calculation"""
        # Setup the CHGCAR and remote folder input if necessary
        self.select_chgcar_from_inputs()
        # Only need to run SCF calculation when no explicity CHGCAR or folder set
        return not (self.ctx.chgcar or self.ctx.restart_folder)

    def should_generate_path(self):
        """
        Seekpath should only run if no explicit bands is provided or we are just
        running for DOS, in which case the original structure is used.
        """
        return 'bs_kpoints' not in self.inputs and (not self.inputs.get('only_dos', False))

    def generate_path(self):
        """
        Run seekpath to obtain the primitive structure and bands
        """

        current_structure_backup = self.ctx.current_structure

        mode = self.inputs.get('band_mode')
        if mode is None:
            mode = self.DEFAULT_BAND_MODE
        else:
            mode = mode.value

        if mode == 'seekpath-aiida':
            inputs = {
                'reference_distance': self.inputs.get('bands_kpoints_distance', None),
                'metadata': {
                    'call_link_label': 'seekpath'
                },
            }
            func = seekpath_structure_analysis
        else:
            # Using sumo interface
            from .common.sumo_kpath import kpath_from_sumo
            inputs = {
                'line_density': self.inputs.get('line_density', orm.Float(self.DEFAULT_LINE_DENSITY)),
                'symprec': self.inputs.get('symprec', orm.Float(self.DEFAULT_SYMPREC)),
                'mode': orm.Str(mode),
                'metadata': {
                    'call_link_label': 'sumo_kpath'
                },
            }
            func = kpath_from_sumo

        magmom = self.ctx.get('magmom', None)

        # For magnetic structures, create different kinds for the analysis in case that the
        # symmetry should be lowered. This also makes sure that the magnetic moments are consistent
        if magmom:
            decorate_result = magnetic_structure_decorate(self.ctx.current_structure, orm.List(list=magmom))
            decorated = decorate_result['structure']
            # Run seekpath on the decorated structure
            kpath_results = func(decorated, **inputs)
            decorated_primitive = kpath_results['primitive_structure']
            # Convert back to undecorated structures and add consistent magmom input
            dedecorate_result = magnetic_structure_dedecorate(decorated_primitive, decorate_result['mapping'])
            self.ctx.magmom = dedecorate_result['magmom'].get_list()
            self.ctx.current_structure = dedecorate_result['structure']
        else:
            kpath_results = func(self.ctx.current_structure, **inputs)
            self.ctx.current_structure = kpath_results['primitive_structure']

        if not np.allclose(self.ctx.current_structure.cell, current_structure_backup.cell):
            if self.inputs.scf.get('kpoints'):
                self.report(
                    'The primitive structure is not the same as the input structure but explicit kpoints are supplied - aborting the workchain.'
                )
                return self.exit_codes.ERROR_INPUT_STRUCTURE_NOT_PRIMITIVE  # pylint: disable=no-member
            self.report(
                'The primitive structure is not the same as the input structure - using the former for all calculations from now.'
            )
        self.ctx.bs_kpoints = kpath_results['explicit_kpoints']
        self.out('primitive_structure', self.ctx.current_structure)
        if 'parameters' in kpath_results:
            self.out('seekpath_parameters', kpath_results['parameters'])

    def run_scf(self):
        """
        Run the SCF calculation
        """

        base_work = WorkflowFactory(self._base_wk_string)
        inputs = AttributeDict(self.exposed_inputs(base_work, namespace='scf'))
        inputs.metadata.call_link_label = 'scf'
        inputs.metadata.label = self.inputs.metadata.label + ' SCF'
        inputs.structure = self.ctx.current_structure

        # Turn off cleaning of the working directory
        clean = inputs.get('clean_workdir')
        if clean and clean.value is False:
            pass
        else:
            inputs.clean_workdir = orm.Bool(False)

        # Ensure that writing the CHGCAR file is on
        pdict = inputs.parameters.get_dict()
        if (pdict[OVERRIDE_NAMESPACE].get('lcharg') is False) or (pdict[OVERRIDE_NAMESPACE].get('LCHARG') is False):
            pdict[OVERRIDE_NAMESPACE]['lcharg'] = True
            inputs.parameters = orm.Dict(dict=pdict)
            self.report('Correction: setting LCHARG to True')

        # Take magmom from the context, in case that the magmom is rearranged in the primitive cell
        magmom = self.ctx.get('magmom')
        if magmom:
            inputs.parameters = nested_update_dict_node(inputs.parameters, {OVERRIDE_NAMESPACE: {'magmom': magmom}})

        running = self.submit(base_work, **inputs)
        self.report(f'Running SCF calculation {running}')
        self.to_context(workchain_scf=running)

    def verify_scf(self):
        """Inspect the SCF calculation"""
        scf_workchain = self.ctx.workchain_scf
        if not scf_workchain.is_finished_ok:
            self.report('SCF workchain finished with Error')
            return self.exit_codes.ERROR_SUB_PROC_SCF_FAILED

        # Store the charge density or remote reference
        if 'chgcar' in scf_workchain.outputs:
            self.ctx.chgcar = scf_workchain.outputs.chgcar
        else:
            self.ctx.chgcar = None
        self.ctx.restart_folder = scf_workchain.outputs.remote_folder
        self.report(f'SCF calculation {scf_workchain} completed')

    def run_bands_dos(self):
        """Run the bands and the DOS calculations"""
        base_work = WorkflowFactory(self._base_wk_string)

        # Use the SCF inputs as the base
        inputs = AttributeDict(self.exposed_inputs(base_work, namespace='scf'))
        inputs.structure = self.ctx.current_structure

        if self.ctx.restart_folder:
            inputs.restart_folder = self.ctx.restart_folder

        if self.ctx.chgcar:
            inputs.chgcar = self.ctx.chgcar

        if not (inputs.get('restart_folder') or inputs.get('chgcar')):
            raise RuntimeError('One of the restart_folder or chgcar must be set for non-scf calculations')

        running = {}

        only_dos = self.inputs.get('only_dos')

        if (only_dos is None) or (only_dos.value is False):
            if 'bands' in self.inputs:
                bands_input = AttributeDict(self.exposed_inputs(base_work, namespace='bands'))
            else:
                bands_input = AttributeDict({
                    'settings': orm.Dict(dict={'parser_settings': {
                        'add_bands': True
                    }}),
                    'parameters': orm.Dict(dict={'charge': {
                        'constant_charge': True
                    }}),
                })

            # Special treatment - combine the parameters
            parameters = inputs.parameters.get_dict()
            bands_parameters = bands_input.parameters.get_dict()

            if 'charge' in bands_parameters:
                bands_parameters['charge']['constant_charge'] = True
            else:
                bands_parameters['charge'] = {'constant_charge': True}

            nested_update(parameters, bands_parameters)

            # Apply updated parameters
            inputs.update(bands_input)
            inputs.parameters = orm.Dict(dict=parameters)

            # Check if add_bands
            settings = inputs.get('settings')
            essential = {'parser_settings': {'add_bands': True}}
            if settings is None:
                inputs.settings = orm.Dict(dict=essential)
            else:
                inputs.settings = nested_update_dict_node(settings, essential)

            # Swap with the default kpoints generated
            inputs.kpoints = self.ctx.bs_kpoints

            # Tag the calculation
            inputs.metadata.label = self.inputs.metadata.label + ' BS'
            inputs.metadata.call_link_label = 'bs'

            bands_calc = self.submit(base_work, **inputs)
            running['bands_workchain'] = bands_calc
            self.report(f'Submitted workchain {bands_calc} for band structure')

        # Do DOS calculation if dos input namespace is populated or a
        # dos_kpoints input is passed.
        if ('dos_kpoints_density' in self.inputs) or ('dos' in self.inputs):
            if 'dos' in self.inputs:
                dos_input = AttributeDict(self.exposed_inputs(base_work, namespace='dos'))
            else:
                dos_input = AttributeDict({
                    'settings': orm.Dict(dict={'add_dos': True}),
                    'parameters': orm.Dict(dict={'charge': {
                        'constant_charge': True
                    }}),
                })
            # Use the supplied kpoints density for DOS
            if 'dos_kpoints_density' in self.inputs:
                dos_kpoints = orm.KpointsData()
                dos_kpoints.set_cell_from_structure(self.ctx.current_structure)
                dos_kpoints.set_kpoints_mesh_from_density(self.inputs.dos_kpoints_density.value * 2 * np.pi)
                dos_input.kpoints = dos_kpoints

            # Special treatment - combine the parameters
            parameters = inputs.parameters.get_dict()
            dos_parameters = dos_input.parameters.get_dict()
            nested_update(parameters, dos_parameters)

            # Ensure we start from constant charge
            if 'charge' in dos_parameters:
                dos_parameters['charge']['constant_charge'] = True
            else:
                dos_parameters['charge'] = {'constant_charge': True}

            # Apply updated parameters
            inputs.update(dos_input)
            inputs.parameters = orm.Dict(dict=parameters)

            if 'dos' not in self.inputs:
                # kindly add `add_dos` if the `dos` input namespace is not
                # explicitly defined.
                settings = inputs.get('settings')
                essential = {'parser_settings': {'add_dos': True}}

                if settings is None:
                    inputs.settings = orm.Dict(dict=essential)
                else:
                    inputs.settings = nested_update_dict_node(settings, essential)

            # Set the label
            inputs.metadata.label = self.inputs.metadata.label + ' DOS'
            inputs.metadata.call_link_label = 'dos'

            dos_calc = self.submit(base_work, **inputs)
            running['dos_workchain'] = dos_calc
            self.report(f'Submitted workchain {dos_calc} for DOS')

        return self.to_context(**running)

    def inspect_bands_dos(self):
        """Inspect the bands and dos calculations"""

        exit_code = None

        if 'bands_workchain' in self.ctx:
            bands = self.ctx.bands_workchain
            if not bands.is_finished_ok:
                self.report(f'Bands calculation finished with error, exit_status: {bands}')
                exit_code = self.exit_codes.ERROR_SUB_PROC_BANDS_FAILED
            self.out(
                'band_structure',
                compose_labelled_bands(bands.outputs.bands, bands.inputs.kpoints),
            )
        else:
            bands = None

        if 'dos_workchain' in self.ctx:
            dos = self.ctx.dos_workchain
            if not dos.is_finished_ok:
                self.report(f'DOS calculation finished with error, exit_status: {dos.exit_status}')
                exit_code = self.exit_codes.ERROR_SUB_PROC_DOS_FAILED

            # Attach outputs
            self.out('dos', dos.outputs.dos)
            if 'projectors' in dos.outputs:
                self.out('projectors', dos.outputs.projectors)
        else:
            dos = None

        return exit_code

    def on_terminated(self):
        """
        Clean the remote directories of all called childrens
        """

        super().on_terminated()

        if self.inputs.clean_children_workdir.value != 'none':
            cleaned_calcs = []
            for called_descendant in self.node.called_descendants:
                if isinstance(called_descendant, orm.CalcJobNode):
                    try:
                        called_descendant.outputs.remote_folder._clean()  # pylint: disable=protected-access
                        cleaned_calcs.append(called_descendant.pk)
                    except (OSError, KeyError):
                        pass

            if cleaned_calcs:
                self.report(f"cleaned remote folders of calculations: {' '.join(map(str, cleaned_calcs))}")


@calcfunction
def seekpath_structure_analysis(structure, **kwargs):
    """Primitivize the structure with SeeKpath and generate the high symmetry k-point path through its Brillouin zone.
    This calcfunction will take a structure and pass it through SeeKpath to get the normalized primitive cell and the
    path of high symmetry k-points through its Brillouin zone. Note that the returned primitive cell may differ from the
    original structure in which case the k-points are only congruent with the primitive cell.
    The keyword arguments can be used to specify various Seekpath parameters, such as:
        with_time_reversal: True
        reference_distance: 0.025
        recipe: 'hpkot'
        threshold: 1e-07
        symprec: 1e-05
        angle_tolerance: -1.0
    Note that exact parameters that are available and their defaults will depend on your Seekpath version.
    """
    from aiida.tools import get_explicit_kpoints_path

    # All keyword arugments should be `Data` node instances of base type and so should have the `.value` attribute
    unwrapped_kwargs = {key: node.value for key, node in kwargs.items() if isinstance(node, orm.Data)}

    return get_explicit_kpoints_path(structure, **unwrapped_kwargs)


@calcfunction
def compose_labelled_bands(bands, kpoints):
    """
    Add additional information from the kpoints allow richer informations
    to be stored such as band structure labels.
    """
    new_bands = deepcopy(bands)
    new_bands.set_kpointsdata(kpoints)
    return new_bands


@calcfunction
def get_primitive_strucrture_and_scf_kpoints(structure):
    """
    This function dryruns a VASP calculation using the primitive structure obtained by performing seekpath analyses

    The input StructureData should be returned by an VaspRelaxWorkChain which will be used for dryun using local
    VASP and getting the explicity kpoints for SCF calculation.
    """
    # Locate the relaxation work
    from aiida.tools import get_explicit_kpoints_path

    from .common.dryrun import dryrun_relax_builder

    # Locate the relaxation work
    relax_work = structure.get_incoming(link_label_filter='relax__structure').one().node
    primitive = get_explicit_kpoints_path(structure)['primitive_structure']

    # Create an restart builder
    builder = relax_work.get_builder_restart()
    builder.structure = primitive

    # Dryrun and construct the SCF kpoints
    kpoint_weights = np.array(dryrun_relax_builder(builder)['kpoints_and_weights'])
    scf_kpoints = orm.KpointsData()
    scf_kpoints.set_kpoints(kpoint_weights[:, :3], weights=kpoint_weights[:, -1])
    return {'primitive': primitive, 'scf_kpoints': scf_kpoints}


class VaspHybridBandsWorkChain(VaspBandsWorkChain):
    """
    Bands workchain for hybrid calculations

    This workchain compute the bandstructure by adding band path segments as zero-weighted
    kpoints for self-consistent calculations. This is mainly for hybrid calculations, but can
    also be used for GGA calculations, although it would be not as efficient as the non-SCF
    method implemented in ``VaspBandsWorkChain``.

    In contrast to ``VaspBandsWorkChain`` this workflow requires and explicitly defined kpoints
    set for the ``scf.kpoints`` port. This can be obtained by parsing the ``IBZKPT`` file from
    and existing calculation or dryrun. Or by parsing the ``vasprun.xml`` file.

    If a relaxation workchain is run as part of the process, the ``kpoints`` output returned can
    be used for this purpose automatically.

    Only the `scf` namespace will be used for performing the calculation

    TODO:
     - Warn if the calculation is not actually a hybrid one
     - Automatic Kpoints from dryruns
    """

    @classmethod
    def define(cls, spec):
        """Initialise the WorkChain class"""
        super().define(spec)
        relax_work = WorkflowFactory(cls._relax_wk_string)
        base_work = WorkflowFactory(cls._base_wk_string)

        spec.input(
            'kpoints_per_split',
            valid_type=orm.Int,
            help='Number of kpoints per split, INCLUDING the weighted SCF kpoints.',
            required=True,
        )
        spec.input('structure', help='The input structure', valid_type=orm.StructureData)
        spec.expose_inputs(
            relax_work,
            namespace='relax',
            exclude=('structure',),
            namespace_options={
                'required': False,
                'populate_defaults': False,
                'help': 'Inputs for Relaxation workchain, if needed',
            },
        )
        spec.expose_inputs(
            base_work,
            namespace='scf',
            exclude=('structure',),
            namespace_options={
                'required': True,
                'populate_defaults': True,
                'help': 'Inputs for SCF workchain, mandatory',
            },
        )
        spec.input(
            'clean_children_workdir',
            valid_type=orm.Str,
            serializer=to_aiida_type,
            help='What part of the called children to clean',
            required=False,
            default=lambda: orm.Str('none'),
        )
        spec.outline(
            cls.setup,
            if_(cls.should_do_relax)(
                cls.run_relax,
                cls.verify_relax,
            ),
            if_(cls.should_generate_path)(cls.generate_path),
            cls.make_splitted_kpoints,  # Split the kpoints
            cls.run_scf_multi,  # Launch split calculation
            cls.inspect_and_combine_bands,  # Combined the band structure
        )

        spec.output(
            'primitive_structure',
            required=False,
            help='Primitive structure used for band structure calculations',
        )
        spec.output('band_structure', required=False, help='Computed band structure with labels')
        spec.output('seekpath_parameters', help='Parameters used by seekpath', required=False)

        spec.exit_code(501, 'ERROR_SUB_PROC_RELAX_FAILED', message='Relaxation workchain failed')
        spec.exit_code(502, 'ERROR_SUB_PROC_SCF_FAILED', message='SCF workchain failed')
        spec.exit_code(
            503,
            'ERROR_SUB_PROC_BANDS_FAILED',
            message='Band structure workchain failed',
        )
        spec.exit_code(504, 'ERROR_SUB_PROC_DOS_FAILED', message='DOS workchain failed')
        spec.exit_code(
            505,
            'ERROR_NO_VALID_SCF_KPOINTS_INPUT',
            message='Cannot found valid inputs for SCF kpoints',
        )
        spec.exit_code(
            601,
            'ERROR_INPUT_STRUCTURE_NOT_PRIMITIVE',
            message='The input structure is not the primitive one!',
        )

    def make_splitted_kpoints(self):
        """Split the kpoints"""
        # Fully specified band structure kpoints
        full_kpoints = self.ctx.bs_kpoints

        if 'kpoints' in self.inputs.scf:
            scf_kpoints = self.inputs.scf.kpoints
        # Relaxation workchain has kpoints output
        elif ('workchain_relax' in self.ctx and 'kpoints' in self.ctx['workchain_relax'].outputs):
            scf_kpoints = self.ctx.workchain_relax.outputs.kpoints
            self.report(f'Using output from <{self.ctx.workchain_relax}> for SCF kpoints.')
        # Parse from relaxation output
        elif 'workchain_relax' in self.ctx:
            # Try getting the kpoints from the retrieved folder
            scf_kpoints = extract_kpoints_from_calc(self.ctx.workchain_relax)
            self.report(f'Extracted SCF kpoints from retrieved vasprun.xml of <{self.ctx.workchain_relax}>.')
        else:
            self.report('No valid SCF kpoints is avaliable to use. Please define scf.kpoints explicitly!')
            return self.exit_codes.ERROR_NO_VALID_SCF_KPOINTS_INPUT  # pylint: disable=no-member

        # Number of kpoints per split, NOT including the SCF kpoints
        per_split = orm.Int(self.inputs.kpoints_per_split.value - scf_kpoints.get_kpoints().shape[0])
        kpoints_for_calc = split_kpoints(scf_kpoints, full_kpoints, per_split)
        self.ctx.kpoints_for_calc = kpoints_for_calc

    def run_scf_multi(self):
        """
        Launch multiple SCF calculations with zero-weighted kpoints for segments of the band structure
        """

        workflow_class = WorkflowFactory(self._base_wk_string)
        for key, value in self.ctx.kpoints_for_calc.items():
            idx = int(key.split('_')[-1])

            inputs = self.exposed_inputs(workflow_class, 'scf')

            # Ensure that the bands are parsed
            if 'settings' not in inputs:
                inputs.settings = orm.Dict(dict={'parser_settings': {'add_bands': True}})
            else:
                # Merge with 'parser_settings'
                inputs.settings = nested_update_dict_node(inputs.settings, {'parser_settings': {'add_bands': True}})

            # Swap the kpoints the the one with zero-weight parts
            inputs.kpoints = value
            inputs.metadata.label = self.inputs.metadata.label + f' SPLIT {idx}'
            inputs.metadata.call_link_label = f'bandstructure_split_{idx:03d}'
            inputs.structure = self.ctx.current_structure
            running = self.submit(workflow_class, **inputs)
            self.report(f'launching {workflow_class.__name__}<{running.pk}> for split #{idx}')
            self.to_context(workchains=append_(running))

    def inspect_and_combine_bands(self):
        """
        Inspect that all calculations have finished OK
        """
        workchains = self.ctx.workchains

        return_codes = [work.exit_status for work in workchains]
        if any(return_codes):
            self.report('At least one calculation did not have zero return code!')

        # Extract the bands information
        self.report(f'Extracting output bandstructure from {len(self.ctx.workchains)} workchains.')
        kwargs = {}
        for work in workchains:
            link_label = (work.get_incoming(link_type=LinkType.CALL_WORK).one().link_label)
            link_idx = int(link_label.split('_')[-1])
            kwargs[f'band_{link_idx:03d}'] = work.outputs.bands
            kwargs[f'kpoint_{link_idx:03d}'] = work.inputs.kpoints

        combined_bands = combine_bands_data(self.ctx.bs_kpoints, **kwargs)
        self.out('band_structure', combined_bands)


@calcfunction
def split_kpoints(scf_kpoints, band_kpoints, kpn_per_split):
    """
    Split the kpoints into multiple one and combined with SCF kpoints

    The kpoints for band structure calculation has zero weights
    """
    return _split_kpoints(scf_kpoints, band_kpoints, kpn_per_split)


def _split_kpoints(scf_kpoints: orm.KpointsData, band_kpoints: orm.KpointsData, kpn_per_split: orm.Int):
    """
    Split the kpoints into multiple one and combined with SCF kpoints

    The kpoints for band structure calculation has zero weights
    """
    scf_kpoints_array, scf_weights_array = scf_kpoints.get_kpoints(also_weights=True)
    band_kpn = band_kpoints.get_kpoints()
    nband_kpts = band_kpn.shape[0]
    nscf_kpts = scf_kpoints_array.shape[0]

    # Split the kpoints
    kpn_per_split = int(kpn_per_split)
    kpt_splits = [band_kpn[i:i + kpn_per_split] for i in range(0, nband_kpts, kpn_per_split)]

    splitted_kpoints = {}
    for isplit, skpts in enumerate(kpt_splits):
        kpt = orm.KpointsData()
        kpt_array = np.concatenate([scf_kpoints_array, skpts], axis=0)
        weights_array = np.zeros(kpt_array.shape[0])
        # Set the weights for SCF kpoints
        weights_array[:nscf_kpts] = scf_weights_array
        # Set kpoints and the weights
        kpt.set_kpoints(kpt_array, weights=weights_array)
        kpt.label = f'SPLIT {isplit:03d}'
        kpt.description = 'Splitted kpoints'
        splitted_kpoints[f'bs_kpoints_{isplit:03d}'] = kpt

    return splitted_kpoints


def dryrun_split_kpoints(
    structure: orm.StructureData,
    scf_kpoints: orm.KpointsData,
    kpn_per_split: orm.Int,
    kpoints_args=None,
    verbose=True,
):
    """
    Perform a "dryrun" for splitting the kpoints
    """
    from aiida.tools import get_explicit_kpoints_path

    if kpoints_args is None:
        kpoints_args = {}
    seekpath_results = get_explicit_kpoints_path(structure, **kpoints_args)
    explicit_kpoints = seekpath_results['explicit_kpoints']
    splitted = _split_kpoints(scf_kpoints, explicit_kpoints, kpn_per_split)
    if verbose:
        nseg = len(splitted)
        nkpts = [kpn.get_kpoints().shape[0] for kpn in splitted.values()]
        print(f'Splitted in to {nseg} segements with number of kpoints: {nkpts}')
    return seekpath_results, splitted


@calcfunction
def combine_bands_data(bs_kpoints, **kwargs):
    """
    Combine splitted bands and kpoints data

    The inputs should be supplied as keyword arguments like `band_001`, `kpoint_001` for the splitted
    kpoints and correspdonging bands data from each calculation.
    The `bs_kpoints` is the originally generated band structure path.

    Returns a `BandsData` by combining the zero-weighted bands from each calculation.
    """
    kpoints_list = [[node, int(key.split('_')[1])] for key, node in kwargs.items() if 'kpoint' in key]
    kpoints_list.sort(key=lambda x: x[1])
    kpoints_list = [item[0] for item in kpoints_list]

    bands_list = [[node, int(key.split('_')[1])] for key, node in kwargs.items() if 'band' in key]
    bands_list.sort(key=lambda x: x[1])
    bands_list = [item[0] for item in bands_list]

    return _combine_bands_data(bs_kpoints, kpoints_list, bands_list)


def _combine_bands_data(
    bs_kpoints: orm.KpointsData,
    kpoints_list: List[orm.KpointsData],
    bands_list: List[orm.BandsData],
):
    """
    Combine bands from splitted kpoints into a single bands node.

    The list of kpoints and bands must be sorted in the right order.
    """
    bands_array_combine = []
    occu_array_combine = []
    kpoints_combine = []

    for skpts, sbands in zip(kpoints_list, bands_list):
        kpt_array, weights_array = skpts.get_kpoints(also_weights=True)
        zero_weight_mask = weights_array == 0.0
        kpoints_combine.append(kpt_array[zero_weight_mask, :])

        bands_array = sbands.get_bands()
        if 'occupations' in sbands.get_arraynames():
            occ_array = sbands.get_array('occupations')
        else:
            occ_array = None

        # Bands array can have three or two dimensions, we have to handle it separately
        if bands_array.ndim == 3:
            bands_array_combine.append(bands_array[:, zero_weight_mask, :])
            if occ_array is not None:
                occu_array_combine.append(occ_array[:, zero_weight_mask, :])
        else:
            bands_array_combine.append(bands_array[zero_weight_mask, :])
            if occ_array is not None:
                occu_array_combine.append(occ_array[zero_weight_mask, :])

    # Concatenate arrays
    if bands_array.ndim == 3:
        band_array_full = np.concatenate(bands_array_combine, axis=1)
        if occu_array_combine:
            occu_array_full = np.concatenate(occu_array_combine, axis=1)
        else:
            occu_array_full = None
    else:
        band_array_full = np.concatenate(bands_array_combine, axis=0)
        if occu_array_combine:
            occu_array_full = np.concatenate(occu_array_combine, axis=0)
        else:
            occu_array_full = None

    # Sanity check all valid kpoints should combine into the original path
    all_kpoints = np.concatenate(kpoints_combine, axis=0)
    if not np.allclose(all_kpoints, bs_kpoints.get_kpoints()):
        raise ValueError('The k-path segements do not much the original path when combined!')

    # Compose the node
    band_data = orm.BandsData()
    band_data.set_kpointsdata(bs_kpoints)
    band_data.set_bands(band_array_full, occupations=occu_array_full)

    return band_data


def extract_kpoints_from_calc(calc):
    """
    Extract computed kpoints from a existing calculation
    """
    retrieved = calc.outputs.retrieved
    return extract_kpoints_from_retrieved(retrieved)


@calcfunction
def extract_kpoints_from_retrieved(retrieved):
    """
    Extract explicity kpoints from a finished calculation
    """
    return _extract_kpoints_from_retrieved(retrieved)


def _extract_kpoints_from_retrieved(retrieved):
    """
    Extract explicity kpoints from a finished calculation
    """

    with retrieved.base.repository.open('vasprun.xml', 'rb') as fh:
        parser = VasprunParser(handler=fh)

    vkpoints = parser.kpoints
    if vkpoints['mode'] != 'explicit':
        raise ValueError('Only explicity kpoints is supported!')

    kpoints_array = vkpoints['points']
    weights_array = vkpoints['weights']

    kpoints_data = orm.KpointsData()
    kpoints_data.set_kpoints(kpoints=kpoints_array, weights=weights_array)

    return kpoints_data
