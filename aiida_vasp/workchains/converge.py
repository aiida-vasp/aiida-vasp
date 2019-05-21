# pylint: disable=too-many-lines, too-many-locals, too-many-statements, too-many-public-methods, too-many-branches, attribute-defined-outside-init
"""
ConvergenceWorkChain.

Intended to be used to control convergence checks for plane-wave calculations.
"""
import copy
import numpy as np

from aiida.engine import WorkChain, append_, while_, if_, calcfunction
from aiida.common.extendeddicts import AttributeDict
from aiida.plugins import WorkflowFactory
from aiida.orm.nodes.data.array.bands import find_bandgap

from aiida_vasp.utils.aiida_utils import (get_data_class, get_data_node, displaced_structure, compressed_structure)
from aiida_vasp.utils.workchains import fetch_k_grid, prepare_process_inputs


class ConvergeWorkChain(WorkChain):
    """A workchain to perform convergence tests."""

    _verbose = False
    _next_workchain_string = 'vasp.relax'
    _next_workchain = WorkflowFactory(_next_workchain_string)

    _ALLOWED_CUTOFF_TYPES = {'energy': 0, 'forces': 1, 'vbm': 2, 'gap': 3}

    @classmethod
    def define(cls, spec):
        super(ConvergeWorkChain, cls).define(spec)
        spec.expose_inputs(cls._next_workchain, exclude=('kpoints', 'parameters', 'structure', 'settings'))
        spec.input('parameters', valid_type=get_data_class('dict'))
        spec.input('structure', valid_type=(get_data_class('structure'), get_data_class('cif')))
        spec.input('kpoints', valid_type=get_data_class('array.kpoints'), required=False)
        spec.input('settings', valid_type=get_data_class('dict'), required=False)
        spec.input(
            'encut',
            valid_type=get_data_class('float'),
            required=False,
            help="""
                   The plane-wave cutoff to be used during convergence tests in electron volts.
                   """)
        spec.input(
            'kgrid',
            valid_type=get_data_class('array'),
            required=False,
            help="""
                   The k-point grid to be used during convergence tests.
                   """)
        spec.input(
            'encut_start',
            valid_type=get_data_class('float'),
            required=False,
            default=get_data_node('float', 200.0),
            help="""
                   The plane-wave cutoff in electron volts.
                   """)
        spec.input(
            'encut_step',
            valid_type=get_data_class('float'),
            required=False,
            default=get_data_node('float', 50.0),
            help="""
                   The plane-wave cutoff step (increment) in electron volts.
                   """)
        spec.input(
            'encut_samples',
            valid_type=get_data_class('int'),
            required=False,
            default=get_data_node('int', 10),
            help="""
                   The number of plane-wave cutoff samples.
                   """)
        spec.input(
            'k_step',
            valid_type=get_data_class('float'),
            required=False,
            default=get_data_node('float', 0.1),
            help="""
                   The k-point step size in inverse AA.
                   """)
        spec.input(
            'k_spacing',
            valid_type=get_data_class('float'),
            required=False,
            default=get_data_node('float', 0.5),
            help="""
                   The default k-point spacing in inverse AA.
                   """)
        spec.input(
            'k_samples',
            valid_type=get_data_class('int'),
            required=False,
            default=get_data_node('int', 10),
            help="""
                   The number of k-point samples.
                   """)
        spec.input(
            'cutoff_type',
            valid_type=get_data_class('str'),
            required=False,
            default=get_data_node('str', 'energy'),
            help="""
                   The cutoff_type to check convergence against. Currently the following
                   options are accepted:
                   * energy
                   * gap
                   * vbm (not yet currently supported)
                   * forces
                   """)
        spec.input(
            'cutoff_value',
            valid_type=get_data_class('float'),
            required=False,
            default=get_data_node('float', 0.01),
            help="""
                   The cutoff value to be used. When the difference between two convergence
                   calculations are within this value for ``cutoff_type``, then it is
                   considered converged.
                   """)
        spec.input(
            'cutoff_value_r',
            valid_type=get_data_class('float'),
            required=False,
            default=get_data_node('float', 0.01),
            help="""
                   The relative cutoff value to be used. When the difference between two convergence
                   calculations are within this value for ``cutoff_type``, then it is
                   considered converged. However, in this case the cutoff value is the difference 
                   between `cutoff_type` for the input structure and an atomic displacement or a 
                   compression of the unitcell.
                   """)
        spec.input(
            'compress',
            valid_type=get_data_class('bool'),
            required=False,
            default=get_data_node('bool', False),
            help="""
                   If True, a convergence test of the compressed structure is also
                   performed. The difference of the ``cutoff_type`` values for each
                   calculations are evaluated and when the difference between these are
                   less than ``cutoff_value_r``, the calculation is considered converged.
                   The largest planw-wave cutoff and densest k-point grid are used.
                   """)
        spec.input(
            'displace',
            valid_type=get_data_class('bool'),
            required=False,
            default=get_data_node('bool', False),
            help="""
                   If True, a convergence test of the displaced structure is also
                   performed. The difference of the ``cutoff_type`` values for each
                   calculations are evaluated and when the difference between these are
                   less than ``cutoff_value_r``, the calculation is considered converged.
                   The largest planw-wave cutoff and densest k-point grid are used.
                   """)
        spec.input(
            'displacement_vector',
            valid_type=get_data_class('array'),
            required=False,
            default=default_array('array', np.array([1.0, 1.0, 1.0])),
            help="""
                   The displacement unit vector for the displacement test. Sets the direction
                   of displacement.
                   """)
        spec.input(
            'displacement_distance',
            valid_type=get_data_class('float'),
            required=False,
            default=get_data_node('float', 0.2),
            help="""
                   The displacement distance (L2 norm) for the displacement test in AA. Follows
                   the direction of ``displacement_vector``.
                   """)
        spec.input(
            'displacement_atom',
            valid_type=get_data_class('int'),
            required=False,
            default=get_data_node('int', 1),
            help="""
                   Which atom to displace? Index starts from 1 and follows the sequence for the
                   sites in the Aiida ``structure`` object.
                   """)
        spec.input(
            'volume_change',
            valid_type=get_data_class('array'),
            required=False,
            default=default_array('array', np.array([1.05, 1.05, 1.05])),
            help="""
                   The volume change in direct coordinates for each lattice vector.
                   """)
        spec.input(
            'converge_relax',
            valid_type=get_data_class('bool'),
            required=False,
            default=get_data_node('bool', False),
            help="""
                   If True, we relax for each convergence test.
                   """)
        spec.input(
            'total_energy_type',
            valid_type=get_data_class('str'),
            required=False,
            default=get_data_node('str', 'energy_no_entropy'),
            help="""
                   The energy type that is used when ``cutoff_type`` is set to `energy`.
                   Consult the options available in the parser for the current version.
                   """)
        spec.input(
            'testing',
            valid_type=get_data_class('bool'),
            required=False,
            default=get_data_node('bool', False),
            help="""
                   If True, we assume testing to be performed (e.g. dummy calculations).
                   """)

        spec.outline(
            cls.initialize,
            if_(cls.run_conv_calcs)(
                while_(cls.run_pw_conv_calcs)(
                    cls.init_pw_conv_calc,
                    cls.init_next_workchain,
                    cls.run_next_workchain,
                    cls.results_pw_conv_calc
                ),
                cls.analyze_pw_conv,
                while_(cls.run_kpoints_conv_calcs)(
                    cls.init_kpoints_conv_calc,
                    cls.init_next_workchain,
                    cls.run_next_workchain,
                    cls.results_kpoints_conv_calc
                ),
                cls.init_disp_conv,
                while_(cls.run_pw_conv_disp_calcs)(
                    cls.init_pw_conv_calc,
                    cls.init_next_workchain,
                    cls.run_next_workchain,
                    cls.results_pw_conv_calc
                ),
                if_(cls.analyze_pw_after_disp)(
                    cls.analyze_pw_conv,
                ),
                while_(cls.run_kpoints_conv_disp_calcs)(
                    cls.init_kpoints_conv_calc,
                    cls.init_next_workchain,
                    cls.run_next_workchain,
                    cls.results_kpoints_conv_calc
                ),
                cls.init_comp_conv,
                while_(cls.run_pw_conv_comp_calcs)(
                    cls.init_pw_conv_calc,
                    cls.init_next_workchain,
                    cls.run_next_workchain,
                    cls.results_pw_conv_calc
                ),
                if_(cls.analyze_pw_after_comp)(
                    cls.analyze_pw_conv,
                ),
                while_(cls.run_kpoints_conv_comp_calcs)(
                    cls.init_kpoints_conv_calc,
                    cls.init_next_workchain,
                    cls.run_next_workchain,
                    cls.results_kpoints_conv_calc
                ),
                cls.analyze_conv,
                cls.store_conv,
            ),
            cls.init_converged,
            cls.init_next_workchain,
            cls.run_next_workchain,
            cls.verify_next_workchain,
            cls.results,
            cls.finalize
        )  # yapf: disable

        spec.expose_outputs(cls._next_workchain)
        spec.output('output_convergence_data', valid_type=get_data_class('array'), required=False)
        spec.exit_code(199, 'ERROR_UNKNOWN',
            message='unknown error detected in the restart workchain')
        # spec.output('output_parameters', valid_type=get_data_class('dict'))
        # spec.output('remote_folder', valid_type=get_data_class('remote'))
        # spec.output('retrieved', valid_type=get_data_class('folder'))
        # spec.output('output_structure', valid_type=get_data_class('structure'), required=False)
        # spec.output('output_structure_relaxed', valid_type=get_data_class('structure'), required=False)

        # spec.output('output_kpoints', valid_type=get_data_class('array.kpoints'), required=False)
        # spec.output('output_trajectory', valid_type=get_data_class('array.trajectory'), required=False)
        # spec.output('output_chgcar', valid_type=get_data_class('vasp.chargedensity'), required=False)
        # spec.output('output_wavecar', valid_type=get_data_class('vasp.wavefun'), required=False)
        # spec.output('output_bands', valid_type=get_data_class('array.bands'), required=False)
        # spec.output('output_dos', valid_type=get_data_class('array'), required=False)
        # spec.output('output_occupancies', valid_type=get_data_class('array'), required=False)
        # spec.output('output_energies', valid_type=get_data_class('array'), required=False)
        # spec.output('output_projectors', valid_type=get_data_class('array'), required=False)
        # spec.output('output_dielectrics', valid_type=get_data_class('array'), required=False)
        # spec.output('output_born_charges', valid_type=get_data_class('array'), required=False)
        # spec.output('output_hessian', valid_type=get_data_class('array'), required=False)
        # spec.output('output_dynmat', valid_type=get_data_class('array'), required=False)
        # spec.output('output_final_forces', valid_type=get_data_class('array'), required=False)
        # spec.output('output_final_stress', valid_type=get_data_class('array'), required=False)

    def initialize(self):
        """Initialize."""
        self._init_context()
        self._init_conv()
        self._init_inputs()
        self._init_settings()

        return

    def _init_context(self):
        """Initialize context variables that are used during the logical flow of the BaseRestartWorkChain."""
        self._init_standard_context()
        self._init_converge_context()

        return

    def _init_standard_context(self):
        """Initialize the standard content of context."""

        self.ctx.exit_status = self.exit_codes.ERROR_UNKNOWN
        self.ctx.workchains = []
        self.ctx.inputs = AttributeDict()
        self.ctx.set_input_nodes = True

        return

    def _init_converge_context(self):
        """Initialize the converge part of the context."""
        self.ctx.converge = AttributeDict()
        self.ctx.converge.settings = AttributeDict()
        self._init_pw_context()
        self._init_kpoints_context()

        return

    def _init_inputs(self):
        """Initialize the inputs."""
        try:
            self._verbose = self.inputs.verbose.value
        except AttributeError:
            pass

    def _init_settings(self):
        """Initialize the settings."""
        # Make sure the parser settings at least contains 'add_bands' and the correct
        # output_params settings.
        if self.run_conv_calcs():
            dict_entry = {'add_bands': True, 'output_params': ['total_energies', 'maximum_force']}
            compress = False
            displace = False
            try:
                compress = self.inputs.compress.value
            except AttributeError:
                pass
            try:
                displace = self.inputs.displace.value
            except AttributeError:
                pass
            if compress or displace:
                dict_entry.update({'add_structure': True})
            if 'settings' in self.inputs:
                settings = AttributeDict(self.inputs.settings.get_dict())
                try:
                    settings.parser_settings.update(dict_entry)
                except AttributeError:
                    settings.parser_settings = dict_entry
            else:
                settings = AttributeDict({'parser_settings': dict_entry})
            self.ctx.inputs.settings = settings
        else:
            if 'settings' in self.inputs:
                self.ctx.inputs.settings = self.inputs.settings

        return

    def _init_pw_context(self):
        """Initialize plane wave cutoff variables and store in context."""
        settings = self.ctx.converge.settings
        self.ctx.running_pw = False
        self.ctx.pw_workchains = []
        self.ctx.converge.pw_data = None
        self.ctx.converge.run_pw_conv_calcs = False
        self.ctx.converge.run_pw_conv_calcs_org = False
        self.ctx.converge.encut_sampling = None
        self.ctx.converge.pw_iteration = 0
        settings.encut = None
        # Set supplied encut
        try:
            settings.encut = copy.deepcopy(self.inputs.encut.value)
        except AttributeError:
            pass
        # Check if encut is supplied in the parameters input, this takes presence over
        # the encut supplied in the inputs
        try:
            parameters_dict = self.inputs.parameters.get_dict()
            encut = parameters_dict.get('encut', None)
            settings.encut = encut
        except AttributeError:
            pass
        # We need a copy of the original encut as we will modify it
        self.ctx.converge.settings.encut_org = copy.deepcopy(settings.encut)

        return

    def _init_kpoints_context(self):
        """Initialize the k-point grid variables and store in context."""
        settings = self.ctx.converge.settings
        self.ctx.running_kpoints = False
        self.ctx.kpoints_workchains = []
        self.ctx.converge.k_data = None
        self.ctx.converge.run_kpoints_conv_calcs = False
        self.ctx.converge.run_kpoints_conv_calcs_org = False
        self.ctx.converge.kpoints_iteration = 0
        self.ctx.converge.k_sampling = None
        settings.kgrid = None
        # We need a special flag that lets us know that we have supplied
        # a k-point grid (e.g. then we do not have access to the grid sampling
        # etc. during user information etc.). Also, the user might want to run
        # plane wave cutoff tests with custom k-point grids. This takes
        # presence over a supplied `kgrid` setting.
        settings.supplied_kmesh = True
        try:
            self.inputs.kpoints
        except AttributeError:
            settings.supplied_kmesh = False
            try:
                settings.kgrid = np.array(self.inputs.kgrid.get_array('array'))
            except AttributeError:
                pass
        # We need a copy of the original kgrid as we will modify it
        if settings.kgrid is not None:
            self.ctx.converge.settings.kgrid_org = np.array(settings.kgrid)
        else:
            self.ctx.converge.settings.kgrid_org = None

        return

    def _init_conv(self):
        """Initialize the convergence tests."""

        # Fetch a temporary StructureData and Dict that we will use throughout,
        # overwrite previous inputs (they are still stored in self.inputs for later ref).
        # Since we cannot execute a calc (that seals the node on completion) we store
        # these in converge instead of input and copy them over when needed.
        if self.inputs.compress.value or self.inputs.displace.value:
            # Only copy if we are going to change the structure
            self.ctx.converge.structure = self.inputs.structure.clone()
        else:
            self.ctx.converge.structure = self.inputs.structure
        # Also create a dummy KpointsData in order to calculate the reciprocal
        # unit cell
        kpoints = get_data_class('array.kpoints')()
        kpoints.set_kpoints_mesh([1, 1, 1])
        kpoints.set_cell_from_structure(self.ctx.converge.structure)
        self.ctx.converge.kpoints = kpoints
        self._init_pw_conv()
        self._init_kpoints_conv()

        return

    def init_rel_conv(self):
        """Initialize the relative convergence tests."""

        # Most of the needed parameters are already set initially by `init_conv`. Here,
        # we only reset counters and clear workchain arrays to prepare for a new batch
        # of convergence tests.
        self.ctx.converge.pw_iteration = 0
        self.ctx.converge.kpoints_iteration = 0

        return

    def init_disp_conv(self):
        """Initialize the displacement convergence tests."""

        converge = self.ctx.converge
        settings = converge.settings
        if self.inputs.displace.value:
            # Make sure we reset the plane wave and k-point tests
            if converge.run_pw_conv_calcs_org:
                converge.run_pw_conv_calcs = True
            if converge.run_kpoints_conv_calcs_org:
                converge.run_kpoints_conv_calcs = True
            self.init_rel_conv()
            # Set the new displaced structure
            converge.structure = self._displace_structure()
            # Set extra information on verbose info
            converge.settings.inform_details = ', using a displaced structure'
        # Also, make sure the data arrays from previous convergence tests are saved
        # in order to be able to calculate the relative convergence
        # criterias later.
        converge.pw_data_org = copy.deepcopy(converge.pw_data)
        converge.k_data_org = copy.deepcopy(converge.k_data)
        # Emtpy arrays
        converge.pw_data = []
        converge.k_data = []
        # Finally, reset k-point grid if plane wave cutoff is not supplied
        if settings.encut_org is None:
            if not settings.supplied_kmesh and settings.kgrid_org is None:
                self._set_default_kgrid()

        return

    def init_comp_conv(self):
        """Initialize the compression convergence tests."""

        converge = self.ctx.converge
        settings = converge.settings
        if self.inputs.compress.value:
            # Make sure we reset the plane wave and k-point tests
            if converge.run_pw_conv_calcs_org:
                converge.run_pw_conv_calcs = True
            if converge.run_kpoints_conv_calcs_org:
                converge.run_kpoints_conv_calcs = True
            self.init_rel_conv()
            # Set the new compressed structure
            converge.structure = self._compress_structure()
            # Set extra information on verbose info
            converge.settings.inform_details = ', using a compressed structure'
        # Also, make sure the data arrays from previous convergence tests are saved
        # in order to be able to calculate the relative convergence criterias later.
        # If we jumped the displacement tests, we have already saved the original data.
        if self.inputs.displace.value:
            converge.pw_data_displacement = copy.deepcopy(converge.pw_data)
            converge.k_data_displacement = copy.deepcopy(converge.k_data)
        # Empty arrays
        converge.pw_data = []
        converge.k_data = []
        # Finally, reset k-point grid if plane wave cutoff is not supplied
        if settings.encut_org is None:
            if not settings.supplied_kmesh and settings.kgrid_org is None:
                self._set_default_kgrid()

        return

    def _init_pw_conv(self):
        """Initialize the plane wave convergence tests."""

        converge = self.ctx.converge
        settings = converge.settings
        supplied_kmesh = settings.supplied_kmesh
        encut_org = settings.encut_org
        kgrid_org = settings.kgrid_org
        encut_start = self.inputs.encut_start.value
        encut_step = self.inputs.encut_step.value
        encut_samples = self.inputs.encut_samples.value
        # Detect what kind of convergence tests that needs to be run.
        if encut_org is None:
            # No encut supplied, run plane wave convergence tests.
            converge.pw_data = []
            # Clone the input parameters if we have no encut,
            # we will eject this into the parameters as we go
            try:
                converge.parameters = self.inputs.parameters.clone()
            except AttributeError:
                converge.parameters = get_data_node('dict')
            if not supplied_kmesh and kgrid_org is None:
                self._set_default_kgrid()
            # Turn on plane wave convergene tests.
            converge.run_pw_conv_calcs = True
            converge.run_pw_conv_calcs_org = True
            # make encut test vector
            converge.encut_sampling = [encut_start + x * encut_step for x in range(encut_samples)]

        return

    def _init_kpoints_conv(self):
        """Initialize the kpoints convergence tests."""

        converge = self.ctx.converge
        settings = converge.settings
        kgrid_org = settings.kgrid_org
        supplied_kmesh = settings.supplied_kmesh
        if not supplied_kmesh and kgrid_org is None:
            converge.k_data = []
            # No kpoint grid supplied, run kpoints convergence tests.
            converge.run_kpoints_conv_calcs = True
            converge.run_kpoints_conv_calcs_org = True

            # Make kpoint test vectors.
            # Usually one expect acceptable convergence with a
            # step size of 0.1/AA, typically:
            # 8 AA lattice vector needs roughly 8 kpoints.
            # 4 AA lattice vector needs roughly 16 kpoints etc.
            # Start convergence test with a step size of 0.5/AA,
            # round values up.
            converge.k_sampling = [x * self.inputs.k_step for x in range(self.inputs.k_samples, 0, -1)]

        return

    def _set_default_kgrid(self):
        """Sets the default k-point grid for plane wave convergence tests."""
        converge = self.ctx.converge
        rec_cell = converge.kpoints.cell
        k_spacing = self.inputs.k_spacing.value
        kgrid = fetch_k_grid(rec_cell, k_spacing)
        converge.settings.kgrid = kgrid
        # Update grid.
        kpoints = get_data_class('array.kpoints')()
        kpoints.set_kpoints_mesh(kgrid)
        kpoints.set_cell_from_structure(converge.structure)
        converge.kpoints = kpoints

    def init_converged(self):
        """Prepare to run the final calculation."""
        # Make sure previous inputs are cleared
        self.ctx.inputs = AttributeDict()
        # Structure should be the same as the initial.
        self.ctx.inputs.structure = self.inputs.structure
        # Same with settings (now we do not do convergence, so any updates
        # from these routines to settings can be skipped)
        try:
            self.ctx.inputs.settings = self.inputs.settings
        except AttributeError:
            pass
        # The plane wave cutoff needs to be updated in the parameters to the set
        # value.
        converged_parameters_dict = self.inputs.parameters.get_dict()
        converged_parameters_dict.update({'encut': self.ctx.converge.settings.encut})
        self.ctx.inputs.parameters = get_data_node('dict', dict=converged_parameters_dict)
        # And finally, the k-point grid needs to be updated to the set value, but
        # only if a kpoint mesh was not supplied
        if not self.ctx.converge.settings.supplied_kmesh:
            kpoints = get_data_class('array.kpoints')()
            kpoints.set_kpoints_mesh(self.ctx.converge.settings.kgrid)
            kpoints.set_cell_from_structure(self.ctx.inputs.structure)
            self.ctx.inputs.kpoints = kpoints
        else:
            self.ctx.inputs.kpoints = self.inputs.kpoints

        self.ctx.running_kpoints = False
        self.ctx.running_pw = False
        if not self.inputs.testing.value:
            self.ctx.set_input_nodes = False

        # inform user
        if self._verbose:
            if not self.ctx.converge.settings.supplied_kmesh:
                self.report('executing a calculation with an assumed converged '
                            'plane wave cutoff of {encut} and a {kgrid0}x{kgrid1}x{kgrid2} '
                            'k-point grid'.format(
                                encut=self.ctx.converge.settings.encut,
                                kgrid0=self.ctx.converge.settings.kgrid[0],
                                kgrid1=self.ctx.converge.settings.kgrid[1],
                                kgrid2=self.ctx.converge.settings.kgrid[2]))
            else:
                self.report('executing a calculation with an assumed converged '
                            'plane wave cutoff of {encut} and a supplied k-point grid'.format(encut=self.ctx.converge.settings.encut))

        return

    def _set_input_nodes(self):
        """Replaces the ctx.input nodes from the previous calculations."""
        self.ctx.inputs.structure = self.ctx.converge.structure.clone()
        if self.ctx.converge.settings.encut_org is None or self.inputs.testing.value:
            self.ctx.inputs.parameters = self.ctx.converge.parameters.clone()
        else:
            self.ctx.inputs.parameters = self.inputs.parameters
        # Only the k-points if no mesh was supplied
        if not self.ctx.converge.settings.supplied_kmesh:
            self.ctx.inputs.kpoints = self.ctx.converge.kpoints.clone()
        else:
            self.ctx.inputs.kpoints = self.inputs.kpoints
        return

    def init_next_workchain(self):
        """Initialize the next workchain calculation."""

        try:
            self.ctx.inputs
        except AttributeError:
            raise ValueError('no input dictionary was defined in self.ctx.inputs')

        # Add exposed inputs
        self.ctx.inputs.update(self.exposed_inputs(self._next_workchain))

        # Check if we want to perform relaxation, if so, modify input
        try:
            relax = self.inputs.converge_relax.value
            if relax:
                self.ctx.inputs.relax = get_data_node('bool', True)
        except AttributeError:
            pass

        # If we are running tests, set the system flag in parameters to contain
        # information, such that it is possible to locate different runs
        if self.inputs.testing.value:
            self.report('TESTING')
            settings = self.ctx.converge.settings
            param_dict = self.inputs.parameters.get_dict()
            if not self.ctx.running_kpoints and not self.ctx.running_pw:
                # Converged run, so a special case
                if settings.encut_org is None and settings.supplied_kmesh:
                    location = 'test-case:test_converge_wc/pw'
                elif settings.encut_org is not None and not settings.supplied_kmesh:
                    location = 'test-case:test_converge_wc/kgrid'
                else:
                    location = 'test-case:test_converge_wc/both'
            else:
                if settings.encut_org is None and settings.supplied_kmesh:
                    location = 'test-case:test_converge_wc/pw/' + str(int(settings.encut))
                elif settings.encut_org is not None and not settings.supplied_kmesh:
                    location = 'test-case:test_converge_wc/kgrid/' + str(settings.kgrid[0]) + '_' + str(settings.kgrid[1]) + '_' + str(
                        settings.kgrid[2])
                else:
                    location = 'test-case:test_converge_wc/both/' + str(int(settings.encut)) + '_' + str(settings.kgrid[0]) + '_' + str(
                        settings.kgrid[1]) + '_' + str(settings.kgrid[2])
            param_dict['system'] = location
            self.ctx.converge.parameters = get_data_node('dict', dict=param_dict)

        # Set input nodes
        if self.ctx.set_input_nodes:
            self._set_input_nodes()

        # Make sure we do not have any floating dict (convert to Dict) in the input
        self.ctx.inputs = prepare_process_inputs(self.ctx.inputs)

    def run_next_workchain(self):
        """Run next workchain."""
        inputs = self.ctx.inputs
        running = self.submit(self._next_workchain, **inputs)
        if hasattr(running, 'pid'):
            self.report('launching {}<{}> '.format(self._next_workchain.__name__, running.pid))
        else:
            # Aiida < 1.0
            self.report('launching {}<{}> '.format(self._next_workchain.__name__, running.pk))

        if self.ctx.running_pw:
            return self.to_context(pw_workchains=append_(running))
        elif self.ctx.running_kpoints:
            return self.to_context(kpoints_workchains=append_(running))

        return self.to_context(workchains=append_(running))

    def run_pw_conv_calcs(self):
        """Should a new plane wave cutoff convergence calculation run?"""
        return self.ctx.converge.run_pw_conv_calcs

    def run_pw_conv_disp_calcs(self):
        """Should a new plane wave cutoff displacement convergence calculation run?"""

        return self.ctx.converge.run_pw_conv_calcs and self.inputs.displace.value

    def run_pw_conv_comp_calcs(self):
        """Should a new plane wave cutoff compression convergence calculation run?"""

        return self.ctx.converge.run_pw_conv_calcs and self.inputs.compress.value

    def run_kpoints_conv_calcs(self):
        """Should a new kpoints convergence calculation run?"""
        return self.ctx.converge.run_kpoints_conv_calcs

    def run_kpoints_conv_disp_calcs(self):
        """Should a new kpoints displacement convergence calculation run?"""

        return self.ctx.converge.run_kpoints_conv_calcs and self.inputs.displace.value

    def run_kpoints_conv_comp_calcs(self):
        """Should a new kpoints compression convergence calculation run?"""

        return self.ctx.converge.run_kpoints_conv_calcs and self.inputs.compress.value

    def init_pw_conv_calc(self):
        """Initialize a single plane wave convergence calculation."""

        # Update the plane wave cutoff
        encut = self.ctx.converge.encut_sampling[self.ctx.converge.pw_iteration]
        self.ctx.converge.settings.encut = encut
        parameters_dict = self.ctx.converge.parameters.get_dict()
        parameters_dict.update({'encut': self.ctx.converge.settings.encut})
        self.ctx.converge.parameters = get_data_node('dict', dict=parameters_dict)
        self.ctx.running_pw = True
        self.ctx.running_kpoints = False
        inform_details = self.ctx.converge.settings.get('inform_details')
        if inform_details is None:
            inform_details = ''

        # inform user
        if self._verbose:
            if self.ctx.converge.settings.supplied_kmesh:
                self.report('running plane wave convergence test on the supplied k-point '
                            'mesh for a plane wave cutoff {encut}'.format(encut=encut) + inform_details)
            else:
                self.report('running plane wave convergence test for k-point sampling '
                            'of {kgrid0}x{kgrid1}x{kgrid2} for a plane wave cutoff {encut}'.format(
                                kgrid0=self.ctx.converge.settings.kgrid[0],
                                kgrid1=self.ctx.converge.settings.kgrid[1],
                                kgrid2=self.ctx.converge.settings.kgrid[2],
                                encut=encut) + inform_details)

        return

    def results_pw_conv_calc(self):
        """Fetch and store the relevant convergence parameters for each plane wave calculation."""

        # Check if child workchain was successfull
        exit_status = self.ctx.pw_workchains[-1].exit_status
        if exit_status:
            self.report('This single convergence calculation has to be considered '
                        'failed as the exit status from the child {} is {}'.format(self._next_workchain, exit_status))

        # Update plane wave iteration index.
        self.ctx.converge.pw_iteration += 1
        # Check if the index has an entry, if not, do not perform further
        # calculations.
        try:
            self.ctx.converge.encut_sampling[self.ctx.converge.pw_iteration]
        except IndexError:
            self.ctx.converge.run_pw_conv_calcs = False

        try:
            workchain = self.ctx.pw_workchains[-1]
        except IndexError:
            self.report('the plane wave convergence calculation finished ' 'without returning a {}'.format(self._next_workchain.__name__))

        encut = self.ctx.converge.settings.encut
        if not exit_status:
            parameters = workchain.outputs.output_parameters.get_dict()
            # fetch total energy
            total_energy = parameters['total_energies'][self.inputs.total_energy_type.value]

            # fetch max force
            max_force = parameters['maximum_force']

            # fetch bands and occupations
            bands = workchain.outputs.output_bands

            # fetch band
            _, gap = find_bandgap(bands)
            if gap is None:
                gap = 0.0
            # Aiida cannot do VBM, yet, so set to zero for now
            max_valence_band = 0.0

            # add stuff to the converge context
            self.ctx.converge.pw_data.append([encut, total_energy, max_force, max_valence_band, gap])
        else:
            self.ctx.converge.pw_data.append([encut, None, None, None, None])

        return

    def init_kpoints_conv_calc(self):
        """Initialize a single k-point grid convergence calculation."""

        # Fetch k-point grid by using the distance between each point
        kstep = self.ctx.converge.k_sampling[self.ctx.converge.kpoints_iteration]
        rec_cell = self.ctx.converge.kpoints.cell
        kgrid = fetch_k_grid(rec_cell, kstep)
        self.ctx.converge.settings.kgrid = kgrid
        # Update grid.
        kpoints = get_data_class('array.kpoints')()
        kpoints.set_kpoints_mesh(kgrid)
        kpoints.set_cell_from_structure(self.ctx.converge.structure)
        self.ctx.converge.kpoints = kpoints
        self.ctx.running_kpoints = True
        self.ctx.running_pw = False
        inform_details = self.ctx.converge.settings.get('inform_details')
        if inform_details is None:
            inform_details = ''

        # inform user
        if self._verbose:
            self.report('running k-point convergence test for k-point sampling '
                        'of {}x{}x{} for a plane wave cutoff {encut}'.format(
                            kgrid[0], kgrid[1], kgrid[2], encut=self.ctx.converge.settings.encut) + inform_details)

        return

    def results_kpoints_conv_calc(self):
        """Fetch and store the relevant convergence parameters for each k-point grid calculation."""

        # Check if child workchain was successfull
        exit_status = self.ctx.kpoints_workchains[-1].exit_status
        if exit_status:
            self.report('This single convergence calculation has to be considered '
                        'failed as the exit status from a child workchain is not '
                        'zero, exit_status:{}'.format(exit_status))
        # Update kpoints iteration index
        self.ctx.converge.kpoints_iteration += 1
        # Check if the index has an entry, if not, do not perform further
        # calculations
        try:
            self.ctx.converge.k_sampling[self.ctx.converge.kpoints_iteration]
        except IndexError:
            self.ctx.converge.run_kpoints_conv_calcs = False

        try:
            workchain = self.ctx.kpoints_workchains[-1]
        except IndexError:
            self.report('the k-point grid convergence calculation ' 'finished without returning a {}'.format(self._next_workchain.__name__))

        kgrid = self.ctx.converge.settings.kgrid
        encut = self.ctx.converge.settings.encut
        if not exit_status:
            parameters = workchain.outputs.output_parameters.get_dict()
            # fetch total energy
            total_energy = parameters['total_energies'][self.inputs.total_energy_type.value]

            # fetch max force
            max_force = parameters['maximum_force']

            # fetch bands and occupations
            bands = workchain.outputs.output_bands
            # fetch band
            _, gap = find_bandgap(bands)
            if gap is None:
                gap = 0.0
            # Aiida cannot do VBM, yet, so set to zero for now
            max_valence_band = 0.0

            # add stuff to the converge context
            self.ctx.converge.k_data.append([kgrid[0], kgrid[1], kgrid[2], encut, total_energy, max_force, max_valence_band, gap])
        else:
            self.ctx.converge.k_data.append([kgrid[0], kgrid[1], kgrid[2], encut, None, None, None, None])

        return

    def analyze_pw_after_comp(self):
        """Return True if we are running compressed convergence tests."""
        return self.inputs.compress.value

    def analyze_pw_after_disp(self):
        """Return True if we are running displaced convergence tests."""
        return self.inputs.displace.value

    def analyze_pw_conv(self):
        """Analyze the plane wave convergence and store it if need be."""

        # Only analyze plane wave cutoff if the encut is not supplied
        if self.ctx.converge.settings.encut_org is None:
            encut = self._check_pw_converged()
            # Check if something went wrong
            if encut is None:
                self.report(
                    'We were not able to obtain a convergence of the plane wave cutoff '
                    'to the specified cutoff. This could also be caused by failures of '
                    'the calculations producing results for the convergence tests. Setting '
                    'the plane wave cutoff to the highest specified value: {encut} eV'.format(encut=self.ctx.converge.pw_data[-1][0]))
                self.ctx.converge.settings.encut = self.ctx.converge.pw_data[-1][0]
            else:
                self.ctx.converge.settings.encut = encut

        return

    def _set_encut_and_kgrid(self, encut, kgrid):
        """Sets the encut and kgrid (if mesh was not supplied)."""
        settings = self.ctx.converge.settings
        settings.encut = encut
        if not settings.supplied_kmesh:
            settings.kgrid = kgrid

    def analyze_conv(self):
        """Analyze convergence and store its parameters."""

        settings = self.ctx.converge.settings
        displace = self.inputs.displace.value
        compress = self.inputs.compress.value

        # Notify the user
        if self._verbose:
            self.report('All convergence tests are done.')

        if displace:
            encut_diff_displacement, kgrid_diff_displacement = self._analyze_conv_disp()
            self._set_encut_and_kgrid(encut_diff_displacement, kgrid_diff_displacement)

        if compress:
            # We have data sitting from the compression tests
            self.ctx.converge.pw_data_comp = self.ctx.converge.pw_data
            self.ctx.converge.k_data_comp = self.ctx.converge.k_data
            encut_diff_comp, kgrid_diff_comp = self._analyze_conv_comp()
            self._set_encut_and_kgrid(encut_diff_comp, kgrid_diff_comp)

        if displace and compress:
            encut_disp_comp, kgrid_disp_comp = self._analyze_conv_disp_comp(encut_diff_displacement, encut_diff_comp,
                                                                            kgrid_diff_displacement, kgrid_diff_comp)
            self._set_encut_and_kgrid(encut_disp_comp, kgrid_disp_comp)

        if not (displace or compress):
            encut, kgrid = self._analyze_conv()
            self._set_encut_and_kgrid(encut, kgrid)

        # Check if any we have None entries for encut or kgrid, which means something failed,
        # or that we where not able to reach the requested convergence.
        if settings.encut is None:
            self.report(
                'We were not able to obtain a convergence of the plane wave cutoff '
                'to the specified cutoff. This could also be caused by failures of '
                'the calculations producing results for the convergence tests. Setting '
                'the plane wave cutoff to the highest specified value: {encut} eV'.format(encut=self.ctx.converge.pw_data_org[-1][0]))
            settings.encut = self.ctx.converge.pw_data_org[-1][0]
        if not settings.supplied_kmesh and self.ctx.converge.settings.kgrid is None:
            self.report(
                'We were not able to obtain a convergence of the k-point grid '
                'to the specified cutoff. This could also be caused by failures of '
                'the calculations producing results for the convergence tests. Setting '
                'the k-point grid sampling to the highest specified value: {kgrid}'.format(kgrid=self.ctx.converge.k_data_org[-1][0:3]))
            settings.kgrid = self.ctx.converge.k_data_org[-1][0:3]

        return

    def _analyze_conv(self):
        """
        Analyze convergence using no displacements or compression.

        Note that, in the case of no displacements or compressions, the
        converged plane wave cutoff is already stored.
        """

        settings = self.ctx.converge.settings
        cutoff_type = self.inputs.cutoff_type.value
        cutoff_value = self.inputs.cutoff_value.value

        encut = settings.encut
        k_data = self.ctx.converge.k_data
        if self._verbose:
            self.report('No atomic displacements or compression were performed.' 'The convergence test suggests:')
        if settings.encut_org is None:
            if self._verbose:
                self.report('plane wave cutoff: {encut} eV.'.format(encut=encut))
        else:
            if self._verbose:
                self.report('plane wave cutoff: User supplied.')

        if not settings.supplied_kmesh:
            kgrid = self._check_kpoints_converged(k_data, cutoff_type, cutoff_value)
            if self._verbose:
                if kgrid is not None:
                    self.report('k-point grid: {kgrid0}x{kgrid1}x{kgrid2}'.format(kgrid0=kgrid[0], kgrid1=kgrid[1], kgrid2=kgrid[2]))
                else:
                    self.report('k-point grid: Failed')
        else:
            kgrid = None
            if self._verbose:
                self.report('k-point grid: User supplied')

        if self._verbose:
            self.report('for the convergence criteria {cutoff_type} and a cutoff of {cutoff_value}'.format(
                cutoff_type=cutoff_type, cutoff_value=cutoff_value))

        return encut, kgrid

    def _analyze_conv_disp_comp(self, encut_displacement, encut_comp, kgrid_displacement, kgrid_comp):  # noqa: MC0001
        """
        Analyze the convergence when both displacements and compression is performed.

        We take the maximum of the plane wave cutoff and the densest k-point grid as
        the recommended values.

        """

        cutoff_type = self.inputs.cutoff_type.value
        cutoff_value = self.inputs.cutoff_value_r.value
        # return the highest plane wave cutoff and densest grid (L2 norm)
        # of the two
        encut = max(encut_displacement, encut_comp)
        if self._verbose:
            self.report('The convergence tests, taking the highest required plane-wave and '
                        'k-point values for both the atomic displacement and compression '
                        'tests suggests:')

        if not self.ctx.converge.settings.supplied_kmesh:
            if np.sqrt(sum([x**2 for x in kgrid_displacement])) > np.sqrt(sum([x**2 for x in kgrid_comp])):
                kgrid = kgrid_displacement
            else:
                kgrid = kgrid_comp

        if self.ctx.converge.settings.encut_org is None and encut_displacement is not None and encut_comp is not None:
            if self._verbose:
                self.report('plane wave cutoff: {encut} eV'.format(encut=encut))
        elif self.ctx.converge.settings.encut_org is not None:
            if self._verbose:
                self.report('plane wave cutoff: User supplied')
        else:
            if self._verbose:
                self.report('plane wave cutoff: Failed')

        if not self.ctx.converge.settings.supplied_kmesh and kgrid_displacement is not None and kgrid_comp is not None:
            if self._verbose:
                self.report('k-point grid: {kgrid0}x{kgrid1}x{kgrid2}'.format(kgrid0=kgrid[0], kgrid1=kgrid[1], kgrid2=kgrid[2]))
        elif self.ctx.converge.settings.supplied_kmesh:
            if self._verbose:
                self.report('k-point grid: User supplied')
        else:
            if self._verbose:
                self.report('k-point grid: Failed')

        if self._verbose:
            self.report('for the convergence criteria '
                        '{cutoff_type} and a cutoff of {cutoff_value}.'.format(cutoff_type=cutoff_type, cutoff_value=cutoff_value))

        return encut, kgrid

    def _analyze_conv_disp(self):  # noqa: MC000
        """Analyze the convergence when atomic displacements are performed."""
        settings = self.ctx.converge.settings
        encut_org = settings.encut_org
        kgrid_org = settings.kgrid_org
        cutoff_type = self.inputs.cutoff_type.value
        cutoff_value = self.inputs.cutoff_value.value
        cutoff_value_r = self.inputs.cutoff_value_r.value
        pw_data_org = self.ctx.converge.pw_data_org
        k_data_org = self.ctx.converge.k_data_org
        pw_data_displacement = self.ctx.converge.pw_data_displacement
        encut_displacement = self._check_pw_converged(pw_data_displacement, cutoff_type, cutoff_value)
        if not settings.supplied_kmesh:
            k_data_displacement = self.ctx.converge.k_data_displacement
            kgrid_displacement = self._check_kpoints_converged(k_data_displacement, cutoff_type, cutoff_value)
        else:
            kgrid_diff_displacement = None
        # Calculate diffs for the plane wave cutoff
        if encut_org is None:
            pw_data = pw_data_displacement
            for index, _ in enumerate(pw_data):
                pw_data[index][1:] = [
                    pw_data_displacement[index][j + 1] - pw_data_org[index][j + 1] for j in range(len(pw_data_displacement[0]) - 1)
                ]
            encut_diff_displacement = self._check_pw_converged(pw_data, cutoff_type, cutoff_value_r)
        else:
            encut_diff_displacement = encut_org

        # Then for the k points
        if kgrid_org is None and not settings.supplied_kmesh:
            k_data = k_data_displacement
            for index, _ in enumerate(k_data_displacement):
                k_data[index][4:] = [
                    k_data_displacement[index][j + 4] - k_data_org[index][j + 4] for j in range(len(k_data_displacement[0]) - 4)
                ]
            kgrid_diff_displacement = self._check_kpoints_converged(k_data, cutoff_type, cutoff_value_r)
        if self._verbose:
            self.report('Performed atomic displacements.')
            self.report('The convergence test using the difference between ' 'the original and displaced dataset suggests:')
        if encut_org is None and encut_diff_displacement is not None and encut_displacement is not None:
            if self._verbose:
                self.report('plane wave cutoff: {encut_diff_displacement} '
                            '({encut_displacement} for the isolated displacement tests) eV'.format(
                                encut_diff_displacement=encut_diff_displacement, encut_displacement=encut_displacement))
        elif encut_org:
            if self._verbose:
                self.report('plane wave cutoff: User supplied')
        else:
            if self._verbose:
                self.report('plane wave cutoff: Failed')

        if not settings.supplied_kmesh and kgrid_diff_displacement is not None and kgrid_displacement is not None:
            if self._verbose:
                self.report('a k-point grid of {kgrid_diff_displacement0}x{kgrid_diff_displacement1}'
                            'x{kgrid_diff_displacement2} ({kgrids0}x{kgrids1}x{kgrids2} for the '
                            'isolated displacement tests)'.format(
                                kgrid_diff_displacement0=kgrid_diff_displacement[0],
                                kgrid_diff_displacement1=kgrid_diff_displacement[1],
                                kgrid_diff_displacement2=kgrid_diff_displacement[2],
                                kgrids0=kgrid_displacement[0],
                                kgrids1=kgrid_displacement[1],
                                kgrids2=kgrid_displacement[2]))
        elif settings.supplied_kmesh:
            if self._verbose:
                self.report('k-point grid: User supplied')
        else:
            if self._verbose:
                self.report('k-point grid: Failed')

        if self._verbose:
            self.report('for the convergence criteria {cutoff_type} and a cutoff '
                        'of {cutoff_value_r} ({cutoff_value} for the isolated displacement tests).'.format(
                            cutoff_type=cutoff_type, cutoff_value_r=cutoff_value_r, cutoff_value=cutoff_value))

        return encut_diff_displacement, kgrid_diff_displacement

    def _analyze_conv_comp(self):  # noqa: MC0001
        """Analize the relative convergence due to unit cell compression."""

        settings = self.ctx.converge.settings
        encut_org = settings.encut_org
        kgrid_org = settings.kgrid_org
        cutoff_type = self.inputs.cutoff_type.value
        cutoff_value = self.inputs.cutoff_value.value
        cutoff_value_r = self.inputs.cutoff_value_r.value
        pw_data_org = self.ctx.converge.pw_data_org
        k_data_org = self.ctx.converge.k_data_org
        pw_data_comp = self.ctx.converge.pw_data_comp
        encut_comp = self._check_pw_converged(pw_data_comp, cutoff_type, cutoff_value)
        if not settings.supplied_kmesh:
            k_data_comp = self.ctx.converge.k_data_comp
            kgrid_comp = self._check_kpoints_converged(k_data_comp, cutoff_type, cutoff_value)
        else:
            kgrid_diff_comp = None
        # Calculate diffs for encut
        if encut_org is None:
            pw_data = pw_data_comp
            for index, _ in enumerate(pw_data):
                pw_data[index][1:] = [pw_data_comp[index][j + 1] - pw_data_org[index][j + 1] for j in range(len(pw_data_comp[0]) - 1)]
            encut_diff_comp = self._check_pw_converged(pw_data, cutoff_type, cutoff_value_r)
        else:
            encut_diff_comp = encut_org
        # Then for the k points
        if kgrid_org is None and not settings.supplied_kmesh:
            k_data = k_data_comp
            for index, _ in enumerate(k_data_comp):
                k_data[index][4:] = [k_data_comp[index][j + 4] - k_data_org[index][j + 4] for j in range(len(k_data_comp[0]) - 4)]
            kgrid_diff_comp = self._check_kpoints_converged(k_data, cutoff_type, cutoff_value_r)
        if self._verbose:
            self.report('Performed compression.')
            self.report('The convergence test using the difference between the ' 'original and dataset with a volume change suggests:')
        if encut_org is None and encut_diff_comp is not None and encut_comp is not None:
            if self._verbose:
                self.report('plane wave cutoff: {encut_diff_comp} '
                            '({encut_comp} for the isolated compression tests) eV'.format(
                                encut_diff_comp=encut_diff_comp, encut_comp=encut_comp))
        elif encut_org:
            if self._verbose:
                self.report('plane wave cutoff: User supplied')
        else:
            if self._verbose:
                self.report('plane wave cutoff: Failed')
        if not settings.supplied_kmesh and kgrid_diff_comp is not None and kgrid_comp is not None:
            if self._verbose:
                self.report('k-point grid: {kgrid_diff_comp0}x{kgrid_diff_comp1}x{kgrid_diff_comp2} '
                            '({kgrid_comp0}x{kgrid_comp1}x{kgrid_comp2} for the isolated '
                            'compression tests)'.format(
                                kgrid_diff_comp0=kgrid_diff_comp[0],
                                kgrid_diff_comp1=kgrid_diff_comp[1],
                                kgrid_diff_comp2=kgrid_diff_comp[2],
                                kgrid_comp0=kgrid_comp[0],
                                kgrid_comp1=kgrid_comp[1],
                                kgrid_comp2=kgrid_comp[2]))
        elif settings.supplied_kmesh:
            if self._verbose:
                self.report('k-point grid: User supplied')
        else:
            if self._verbose:
                self.report('k-point grid: Failed')

        if self._verbose:
            self.report('for the convergence criteria {cutoff_type} and a cutoff '
                        'of {cutoff_value_r} ({cutoff_value} for the isolated compression tests).'.format(
                            cutoff_type=cutoff_type, cutoff_value_r=cutoff_value_r, cutoff_value=cutoff_value))

        return encut_diff_comp, kgrid_diff_comp

    @calcfunction
    def store_conv(self):
        """Set up the convergence data and put it in a data node."""
        convergence = get_data_class('array')()

        if self._verbose:
            self.report("attaching the node {}<{}> as '{}'".format(convergence.__class__.__name__, convergence.pk,
                                                                   'output_convergence_data'))

        # Store regular conversion data
        try:
            store_conv_data(convergence, 'pw_regular', self.ctx.converge.pw_data_org)
        except AttributeError:
            store_conv_data(convergence, 'pw_regular', self.ctx.converge.pw_data)

        try:
            store_conv_data(convergence, 'kpoints_regular', self.ctx.converge.k_data_org)
        except AttributeError:
            store_conv_data(convergence, 'kpoints_regular', self.ctx.converge.k_data)

        # Then possibly displacement
        try:
            store_conv_data(convergence, 'pw_displacement', self.ctx.converge.pw_data_displacement)
            store_conv_data(convergence, 'kpoints_displacement', self.ctx.converge.k_data_displacement)
        except AttributeError:
            pass

        # And finally for compression
        try:
            store_conv_data(convergence, 'pw_compression', self.ctx.converge.pw_data_comp)
            store_conv_data(convergence, 'kpoints_compression', self.ctx.converge.k_data_comp)
        except AttributeError:
            pass
        
        self.out('output_convergence_data', convergence)

        return

    def _check_pw_converged(self, pw_data=None, cutoff_type=None, cutoff_value=None):
        """
        Check if plane wave cutoffs are converged to the specified value.

        Returns
        -------
        encut : float
            The converged plane wave cutoff in eV

        """

        if pw_data is None:
            pw_data = self.ctx.converge.pw_data
        if cutoff_type is None:
            cutoff_type = self.inputs.cutoff_type.value
        if cutoff_value is None:
            cutoff_value = self.inputs.cutoff_value.value

        # Make sure we do not analyze entries that have a None entry
        pw_data = [elements for elements in pw_data if None not in elements]
        # Since we are taking deltas, make sure we have at least two entries,
        # otherwise return None
        if len(pw_data) < 2:
            return None
        # Analyze which encut to use further (cutoff_type sets which parameter)
        encut_okey = False
        index = 0
        criteria = self._ALLOWED_CUTOFF_TYPES[cutoff_type]
        for encut in range(1, len(pw_data)):
            delta = abs(pw_data[encut][criteria + 1] - pw_data[encut - 1][criteria + 1])
            if delta < cutoff_value:
                encut_okey = True
                index = encut
                break
        if not encut_okey:
            # if self._verbose:
            #     self.report('Could not obtain convergence for {cutoff_type} with a cutoff '
            #                 'parameter of {cutoff_value}'.format(cutoff_type=cutoff_type, cutoff_value=cutoff_value))
            return None

        return pw_data[index][0]

    def _check_kpoints_converged(self, k_data=None, cutoff_type=None, cutoff_value=None):
        """
        Check if the k-point grid are converged to the specified value.

        Returns
        -------
        kgrid : (3) list of int
            The converged k-point grid sampling in each direction.

        """
        if k_data is None:
            k_data = self.ctx.converge.k_data
        if cutoff_type is None:
            cutoff_type = self.inputs.cutoff_type.value
        if cutoff_value is None:
            cutoff_value = self.inputs.cutoff_value.value

        # Make sure we do not analyze entries that have a None entry
        k_data = [elements for elements in k_data if None not in elements]
        # Since we are taking deltas, make sure we have at least two entries,
        # otherwise return None
        if len(k_data) < 2:
            return None
        # now analyze which k-point grid to use
        k_cut_okey = False
        index = 0
        criteria = self._ALLOWED_CUTOFF_TYPES[cutoff_type]
        for k in range(1, len(k_data)):
            delta = abs(k_data[k][criteria + 4] - k_data[k - 1][criteria + 4])
            if delta < cutoff_value:
                k_cut_okey = True
                index = k
                break
        if not k_cut_okey:
            # self.report('Could not find a dense enough grid to obtain a {cutoff_type} '
            #             'cutoff of {cutoff_value})'.format(cutoff_type=cutoff_type, cutoff_value=cutoff_value))
            return None

        return k_data[index][0:3]

    def verify_next_workchain(self):
        """Verify and inherit exit status from child workchains."""

        workchain = self.ctx.workchains[-1]
        # Adopt exit status from last child workchain (supposed to be
        # successfull)
        next_workchain_exit_status = workchain.exit_status
        if not next_workchain_exit_status:
            self.exit_status = 0
        else:
            self.exit_status = next_workchain_exit_status
            self.report('The child {}<{}> returned a non-zero exit status, {}<{}> '
                        'inherits exit status {}'.format(workchain.__class__.__name__, workchain.pk, self.__class__.__name__, self.pid,
                                                         next_workchain_exit_status))

        return

    def results(self):
        """Attach the remaining output results."""

        if not self.exit_status:
            workchain = self.ctx.workchains[-1]
            self.out_many(self.exposed_outputs(workchain, self._next_workchain))

        return

    def finalize(self):
        """Finalize the workchain."""
        return self.exit_status

    def run_conv_calcs(self):
        """Determines if convergence calcs are to be run at all."""
        return self.run_kpoints_conv_calcs() or self.run_pw_conv_calcs()

    def _displace_structure(self):
        """Displace the input structure according to the supplied settings."""

        displacement_vector = self.inputs.displacement_vector.get_array('array')
        displacement_distance = self.inputs.displacement_distance.value
        displacement_atom = self.inputs.displacement_atom.value
        # Set displacement
        displacement = displacement_distance * displacement_vector

        # Displace and return new structure
        return displaced_structure(self.ctx.converge.structure, displacement, displacement_atom)

    def _compress_structure(self):
        """Compress the input structure according to the supplied settings."""

        volume_change = self.inputs.volume_change.get_array('array')
        # Apply compression and tension
        comp_structure = compressed_structure(self.ctx.converge.structure, volume_change)
        # Make sure we also reset the reciprocal cell
        kpoints = get_data_class('array.kpoints')()
        kpoints.set_kpoints_mesh([1, 1, 1])
        kpoints.set_cell_from_structure(comp_structure)
        self.ctx.converge.kpoints = kpoints

        return comp_structure


def default_array(name, array):
    """Used to set ArrayData for spec.input."""
    array_cls = get_data_node('array')
    array_cls.set_array(name, array)

    return array_cls


def store_conv_data(array, key, data):
    """Store convergence data in the array."""
    if data is not None:
        if data:
            array.set_array(key, np.array(data))

    return
