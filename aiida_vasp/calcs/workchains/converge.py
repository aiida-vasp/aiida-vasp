# pylint: disable=too-many-lines, too-many-locals, too-many-statements, too-many-public-methods, too-many-branches, attribute-defined-outside-init
"""
ConvergenceWorkChain.

Intended to be used to control convergence checks for plane-wave calculations.
"""
import copy
import numpy as np

from aiida.work.workchain import WorkChain, append_, while_, if_
from aiida.common.extendeddicts import AttributeDict
from aiida.orm import WorkflowFactory, Code

from aiida_vasp.utils.aiida_utils import (get_data_class, get_data_node, displaced_structure, compressed_structure, copy_structure,
                                          copy_parameter, copy_kpoints)
from aiida_vasp.utils.workchains import fetch_k_grid, init_input


class ConvergeWorkChain(WorkChain):
    """A workchain to perform convergence tests."""

    _verbose = True
    _next_workchain = WorkflowFactory('vasp.relax')

    # a few default values defined here for now
    _encut_start_default = 200
    _encut_step_default = 50
    _encut_samples_default = 3
    _k_step_default = 0.15
    _k_spacing_default = 0.5
    _k_samples_default = 3
    _cutoff_type_default = 'energy'
    _cutoff_value_default = 0.001
    _cutoff_value_r_default = 0.001
    _compress_default = False
    _displacement_default = False
    _displacement_vector_default = np.array([1.0, 1.0, 1.0])
    _displacement_distance_default = 0.2  # AA
    _displacement_atom_default = 1  # atom number starting from 1
    _volume_change_default = np.array([1.05, 1.05, 1.05])
    _relax_default = False

    _ALLOWED_CUTOFF_TYPES = {'energy': 0, 'forces': 1, 'vbm': 2, 'gap': 3}

    @classmethod
    def define(cls, spec):
        super(ConvergeWorkChain, cls).define(spec)
        spec.input('code', valid_type=Code)
        spec.input('structure', valid_type=(get_data_class('structure'), get_data_class('cif')))
        spec.input('potential_family', valid_type=get_data_class('str'))
        spec.input('potential_mapping', valid_type=get_data_class('parameter'))
        spec.input('incar', valid_type=get_data_class('parameter'))
        spec.input('options', valid_type=get_data_class('parameter'))
        spec.input('kpoints', valid_type=get_data_class('array.kpoints'), required=False)
        spec.input('settings', valid_type=get_data_class('parameter'), required=False)
        spec.input('restart.max_iterations', valid_type=get_data_class('int'), required=False)
        spec.input('restart.clean_workdir', valid_type=get_data_class('bool'), required=False)
        spec.input('verify.max_iterations', valid_type=get_data_class('int'), required=False)
        spec.input('verify.clean_workdir', valid_type=get_data_class('bool'), required=False)
        spec.input('relax.incar', valid_type=get_data_class('parameter'), required=False)
        spec.input('relax.perform', valid_type=get_data_class('bool'), required=False, default=get_data_node('bool', True))
        spec.input('relax.positions', valid_type=get_data_class('bool'), required=False, default=get_data_node('bool', True))
        spec.input('relax.shape', valid_type=get_data_class('bool'), required=False, default=get_data_node('bool', False))
        spec.input('relax.volume', valid_type=get_data_class('bool'), required=False, default=get_data_node('bool', False))
        spec.input('relax.convergence.on', valid_type=get_data_class('bool'), required=False, default=get_data_node('bool', False))
        spec.input('relax.convergence.absolute', valid_type=get_data_class('bool'), required=False, default=get_data_node('bool', False))
        spec.input('relax.convergence.max_iterations', valid_type=get_data_class('int'), required=False, default=get_data_node('int', 5))
        spec.input(
            'relax.convergence.shape.lengths', valid_type=get_data_class('float'), required=False,
            default=get_data_node('float', 0.1))  # in cartesian coordinates
        spec.input(
            'relax.convergence.shape.angles', valid_type=get_data_class('float'), required=False,
            default=get_data_node('float', 0.1))  # in degree in the cartesian system
        spec.input(
            'relax.convergence.volume', valid_type=get_data_class('float'), required=False,
            default=get_data_node('float', 0.01))  # in degree in the cartesian system
        spec.input(
            'relax.convergence.positions', valid_type=get_data_class('float'), required=False,
            default=get_data_node('float', 0.01))  # in degree in the cartesian system

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

        spec.output('output_parameters', valid_type=get_data_class('parameter'))
        spec.output('remote_folder', valid_type=get_data_class('remote'))
        spec.output('retrieved', valid_type=get_data_class('folder'))
        spec.output('output_structure', valid_type=get_data_class('structure'), required=False)
        spec.output('output_structure_relaxed', valid_type=get_data_class('structure'), required=False)
        spec.output('output_convergence_data', valid_type=get_data_class('array'), required=False)
        spec.output('output_kpoints', valid_type=get_data_class('array.kpoints'), required=False)
        spec.output('output_trajectory', valid_type=get_data_class('array.trajectory'), required=False)
        spec.output('output_chgcar', valid_type=get_data_class('vasp.chargedensity'), required=False)
        spec.output('output_wavecar', valid_type=get_data_class('vasp.wavefun'), required=False)
        spec.output('output_bands', valid_type=get_data_class('array.bands'), required=False)
        spec.output('output_dos', valid_type=get_data_class('array'), required=False)
        spec.output('output_occupations', valid_type=get_data_class('array'), required=False)
        spec.output('output_energies', valid_type=get_data_class('array'), required=False)
        spec.output('output_projectors', valid_type=get_data_class('array'), required=False)
        spec.output('output_dielectrics', valid_type=get_data_class('array'), required=False)
        spec.output('output_born_charges', valid_type=get_data_class('array'), required=False)
        spec.output('output_hessian', valid_type=get_data_class('array'), required=False)
        spec.output('output_dynmat', valid_type=get_data_class('array'), required=False)
        spec.output('output_final_forces', valid_type=get_data_class('array'), required=False)
        spec.output('output_final_stress', valid_type=get_data_class('array'), required=False)

    def initialize(self):
        """Initialize."""
        self._init_context()
        self._init_inputs()
        self._init_conv()

        return

    def _init_inputs(self):
        """Initialize the inputs."""
        self.ctx.inputs = init_input(self.inputs)

        return

    def _init_context(self):
        """Initialize context variables that are used during the logical flow of the BaseRestartWorkChain."""
        self._init_standard_context()
        self._init_converge_context()

        return

    def _init_standard_context(self):
        """Initialize the standard content of context."""

        self.ctx.workchains = []
        self.ctx.pw_workchains = []
        self.ctx.kpoints_workchains = []
        self.ctx.running_pw = False
        self.ctx.running_kpoints = False
        self.ctx.replace_nodes = True

        return

    def _init_converge_context(self):
        """Initialize the converge part of the context."""
        self.ctx.converge = AttributeDict()
        self.ctx.converge.settings = AttributeDict()
        self._init_pw_context()
        self._init_kpoints_context()
        self._init_param_context()
        self._check_context()

        return

    def _init_pw_context(self):
        """
        Initialize plane wave cutoff variables and store in context.

        Parameters
        ----------
        encut : float, optional
            The plane wave cutoff in eV. If not set, convergence tests on this
            parameter are intitiated.
        encut_start : float, optional
            The start value of `encut` during convergence tests.
        encut_step : float, optional
            The `encut` step size in eV.
        encut_samples : int
            The number of `encut` steps to perform.
        """
        converge_dict = self.ctx.converge.settings
        converge_dict.encut = None
        converge_dict.encut_start = None
        converge_dict.encut_step = None
        converge_dict.encut_samples = None
        if 'settings' in self.inputs:
            settings = self.inputs.settings.get_dict()
            if 'converge' in settings:
                converge = settings.get('converge')
                converge_dict.encut = converge.get('pw_energy_cutoff')
                converge_dict.encut_start = converge.get('pw_energy_cutoff_start')
                converge_dict.encut_step = converge.get('pw_energy_cutoff_step')
                converge_dict.encut_samples = converge.get('pw_energy_cutoff_samples')
        # Also check if encut is supplied in the incar input, this takes presence over
        # the encut supplied in settings
        if self.inputs.incar:
            incar_dict = self.inputs.incar.get_dict()
            encut = incar_dict.get('encut', None)
            converge_dict.encut = encut

        return

    def _init_kpoints_context(self):
        """
        Initialize the k-point grid variables and store in context.

        Parameters
        ----------
        kgrid : (3) list int, optional
            The k-point grid sampling along each axis, e.g. [3,3,3]. If not set,
            convergence test on the k-point grid is initiated.
        k_step : float, optional
            The k-point step size.
        k_samples : int, optional
            The number of k-point grid refinements to perform with
            the stepping of `k_step`.
        """
        converge_dict = self.ctx.converge.settings
        converge_dict.kgrid = None
        converge_dict.k_step = None
        converge_dict.k_spacing = None
        converge_dict.k_samples = None
        if 'settings' in self.inputs:
            settings = self.inputs.settings.get_dict()
            if 'converge' in settings:
                converge = settings.get('converge')
                converge_dict.kgrid = converge.get('kpoints')
                converge_dict.k_step = converge.get('k_step')
                converge_dict.k_spacing = converge.get('k_spacing')
                converge_dict.k_samples = converge.get('k_samples')
        # We need a special flag that lets us know that we have supplied
        # a k-point grid (e.g. then we do not have access to the grid sampling
        # etc. during user information etc.). Also, the user might want to run
        # plane wave cutoff tests with custom k-point grids. This takes
        # presence over a supplied `kgrid` setting.
        converge_dict.supplied_kmesh = True
        try:
            self.inputs.kpoints
        except AttributeError:
            converge_dict.supplied_kmesh = False

        return

    def _init_param_context(self):
        """
        Initialize general parameters that are used for the convergence tests.

        Parameters
        ----------
        displace : bool, optional
            If True, check also relative convergence by displacing an atom
            (`displacement_*` parameters in the supplied `converge` settings),
            else do not (defaults to True).
        compress : bool, optional
            If True, check also relative convergence by uniaxial change of the unit
            cell volume (`change` parameter in the `converge` settings), else do not
            (defaults to True).
        cutoff_type : ['ENERGY', 'FORCE', 'VBM', 'GAP']
            The type to check the convergence criteria against.
        cutoff_value : float
            The cutoff value for the convergence tests.
        cutoff_value_r : float
            The cutoff value for the relative convergence tests.
        displacement_vector : list
            A list of floats describing the direction of displacement for the atom that
            is to be displaced for the relative displacement convergence test.
        displacement_distance : float
            The vector distance to displacement the atom.
        displacement_atom : int
            The atom number (starting from 1) to displace during the relative
            convergence tests.
        volume_change : list of list
            A list of list of float that describes the relative change for each
            lattice vector for the relative compression displacement tests.
        """

        converge_dict = self.ctx.converge.settings
        converge_dict.compress = None
        converge_dict.displace = None
        converge_dict.cutoff_type = None
        converge_dict.cutoff_value = None
        converge_dict.cutoff_value_r = None
        converge_dict.displacement_vector = None
        converge_dict.displacement_distance = None
        converge_dict.displacement_atom = None
        converge_dict.volume_change = None
        converge_dict.relax = None
        if 'settings' in self.inputs:
            settings = self.inputs.settings.get_dict()
            if 'converge' in settings:
                converge = settings.get('converge')
                converge_dict.compress = converge.get('compress')
                converge_dict.displace = converge.get('displace')
                converge_dict.cutoff_type = converge.get('cutoff_type')
                converge_dict.cutoff_value = converge.get('cutoff_value')
                converge_dict.cutoff_value_r = converge.get('cutoff_value_r')
                converge_dict.displacement_vector = converge.get('displacement_vector')
                converge_dict.displacement_distance = converge.get('displacement_distance')
                converge_dict.displacement_atom = converge.get('displacement_atom')
                converge_dict.volume_change = converge.get('volume_change')
                converge_dict.relax = converge.get('relax')

        return

    def _check_context(self):
        """Check and verify context and set defaults."""

        self._check_context_pw()
        self._check_context_kpoints()
        self._check_context_relative()
        self._check_context_cutoffs()
        self._check_context_others()

        return

    def _check_context_others(self):
        """Check and verify general context parameters."""
        settings = self.ctx.converge.settings

        if settings.relax is None:
            if self._verbose:
                self.report('setting the default relaxation parameter for the convergence ' 'tests to {}'.format(self._relax_default))
                settings.relax = self._relax_default

    def _check_context_pw(self):
        """Check and verify plane wave related parameters."""
        settings = self.ctx.converge.settings

        # Check encut is supplied, if not, also set defaults for convergence
        # tests.
        if settings.encut is None:
            if settings.encut_start is None:
                if self._verbose:
                    self.report('setting the default minimum value of the plane wave '
                                'energy cutoff (pw_energy_cutoff_start) to {}'.format(self._encut_start_default))
                settings.encut_start = self._encut_start_default
            if settings.encut_step is None:
                if self._verbose:
                    self.report('setting the default step size of the plane wave energy '
                                'cutoff (pw_energy_cutoff_step) to {}'.format(self._encut_step_default))
                settings.encut_step = self._encut_step_default
            if settings.encut_samples is None:
                if self._verbose:
                    self.report('setting the default number of plane wave energy '
                                'samples (pw_energy_cutoff_samples) to {}'.format(self._encut_samples_default))
                settings.encut_samples = self._encut_samples_default

        return

    def _check_context_kpoints(self):
        """Check and verify k-point related parameters."""
        settings = self.ctx.converge.settings

        # Check if a k-point sampling is supplied, if not, also set defaults for
        # converge tests.
        if not settings.supplied_kmesh:
            if settings.k_step is None:
                if self._verbose:
                    self.report('setting the default value of the step size between '
                                'k-points (k_step) to {}'.format(self._k_step_default))
                    settings.k_step = self._k_step_default
            if settings.k_samples is None:
                if self._verbose:
                    self.report('setting the default value of the number of k-point samples '
                                'to perform in the convergence test (k_samples) to {}'.format(self._k_samples_default))
                settings.k_samples = self._k_samples_default
            if settings.k_spacing is None:
                if self._verbose:
                    self.report('setting the default value of the k-point spacing ' '(k_spacing) to {}'.format(self._k_spacing_default))
                settings.k_spacing = self._k_spacing_default

        return

    def _check_context_cutoffs(self):
        """Check and verify cutoff related parameters."""
        settings = self.ctx.converge.settings

        if settings.cutoff_type is None:
            if self._verbose:
                self.report('setting the default convergence parameter to the default ' '({})'.format(self._cutoff_type_default))
            settings.cutoff_type = self._cutoff_type_default
        try:
            self._ALLOWED_CUTOFF_TYPES[settings.cutoff_type]
        except KeyError:
            if self._verbose:
                self.report('the supplied convergence parameter is not supported, '
                            'setting it to the default ({})'.format(self._cutoff_type_default))
            settings.cutoff_type = self._cutoff_type_default
        if settings.cutoff_value is None:
            if self._verbose:
                self.report('setting the default convergence value to {}'.format(self._cutoff_value_default))
            settings.cutoff_value = self._cutoff_value_default
        if settings.cutoff_value_r is None:
            if self._verbose:
                self.report('setting the default relative convergence value to {}'.format(self._cutoff_value_r_default))
            settings.cutoff_value_r = self._cutoff_value_r_default

        return

    def _check_context_relative(self):
        """Check and verify existence of relative convergence tests."""
        settings = self.ctx.converge.settings

        if settings.compress is None:
            if self._verbose:
                self.report('turning the relative convergence test for compression on')
            settings.compress = True
        if settings.displace is None:
            if self._verbose:
                self.report('turning the relative convergence test for displacement on')
            settings.displace = True
        self._check_context_relative_st()

        return

    def _check_context_relative_st(self):
        """Check and verify parameters related to the relative convergence tests."""
        settings = self.ctx.converge.settings

        # Only care about displacement parameters if relative displacement
        # convergence tests are performed.
        if settings.displace:
            if settings.displacement_vector is None:
                if self._verbose:
                    self.report('setting the default value of the atomic displacement '
                                'displacement vector for relative convergence tests (displacement_vector) '
                                'to {}'.format(self._displacement_vector_default))
                settings.displacement_vector = self._displacement_vector_default
            if settings.displacement_distance is None:
                if self._verbose:
                    self.report('setting the default value of the atomic displacement '
                                'displacement distance for relative convergence tests (displacement_distance) '
                                'to {}'.format(self._displacement_distance_default))
                settings.displacement_distance = self._displacement_distance_default
            if settings.displacement_atom is None:
                if self._verbose:
                    self.report('setting the default value of the atom to be displaced '
                                'for the relative convergence tests (displacement_atom) '
                                'to {}'.format(self._displacement_atom_default))
                settings.displacement_atom = self._displacement_atom_default
        # Only care about volume parameter if relative compression convergence
        # tests are performed.
        if settings.compress:
            if settings.volume_change is None:
                if self._verbose:
                    self.report('setting the default value of the volume change '
                                'for the relative convergence tests (volume_change) '
                                'to [{}, {}, {}]'.format(self._volume_change_default[0], self._volume_change_default[1],
                                                         self._volume_change_default[2]))
                settings.volume_change = self._volume_change_default

        return

    def _init_conv(self):
        """Initialize the convergence tests."""

        # Fetch a temporary StructureData and ParameterData that we will use throughout,
        # overwrite previous inputs (they are still stored in self.inputs for later ref).
        # Since we cannot execute a calc (that seals the node on completion) we store
        # these in converge instead of input and copy them over when needed.
        self.ctx.converge.structure = copy_structure(self.ctx.inputs.structure)
        self.ctx.converge.incar = copy_parameter(self.ctx.inputs.incar)
        # Also create a dummy KpointsData in order to calculate the reciprocal
        # unit cell
        kpoints = get_data_class('array.kpoints')
        self.ctx.converge.kpoints = kpoints(kpoints_mesh=[1, 1, 1], cell_from_structure=self.ctx.inputs.structure)
        self._init_pw_conv()
        self._init_kpoints_conv()

        # If we do not want relaxations during convergence tests, disable
        # the inputs.relax.perform flag and enable it for the final
        # calculation.
        if not self.ctx.converge.settings.relax:
            self.ctx.inputs.relax.perform = get_data_node('bool', False)

        return

    def init_rel_conv(self):
        """Initialize the relative convergence tests."""

        # Most of the needed parameters are already set initially by `init_conv`. Here,
        # we only reset counters and clear workchain arrays to prepare for a new batch
        # of convergence tests.
        self.ctx.converge.pw_iteration = 0
        self.ctx.converge.pw_workchains = []
        self.ctx.converge.kpoints_iteration = 0
        self.ctx.converge.kpoints_workchains = []

        return

    def init_disp_conv(self):
        """Initialize the displacement convergence tests."""

        if self.ctx.converge.settings.displace:
            self.init_rel_conv()
            # Set the new displaced structure in the context input
            self.ctx.inputs.structure = self._displace_structure()
            # Set extra information on verbose info
            self.ctx.converge.settings.inform_details = ', using a displaced structure'
        # Also, make sure the data arrays from previous convergence tests are saved
        # in order to be able to calculate the relative convergence
        # criterias later.
        self.ctx.converge.pw_data_org = copy.deepcopy(self.ctx.converge.pw_data)
        self.ctx.converge.k_data_org = copy.deepcopy(self.ctx.converge.k_data)

        return

    def init_comp_conv(self):
        """Initialize the compression convergence tests."""

        if self.ctx.converge.settings.compress:
            self.init_rel_conv()
            # Set the new compressed structure
            self.ctx.inputs.structure = self._compress_structure()
            # Set extra information on verbose info
            self.ctx.converge.settings.inform_details = ', using a compressed structure'
        # Also, make sure the data arrays from previous convergence tests are saved
        # in order to be able to calculate the relative convergence criterias later.
        # If we jumped the displacement tests, we have already saved the original data.
        if self.ctx.converge.settings.displace:
            self.ctx.converge.pw_data_displacement = copy.deepcopy(self.ctx.converge.pw_data)
            self.ctx.converge.k_data_displacement = copy.deepcopy(self.ctx.converge.k_data)

        return

    def _init_pw_conv(self):
        """
        Initialize the plane wave convergence tests.

        Parameters
        ----------
        pw_iteration : int
            The current iteration for the plane wave cufoff convergence calculations.
        pw_conv_steps_reached : bool
            True if the last plane wave cutoff energy step has been reached.
        """
        settings = self.ctx.converge.settings
        encut = settings.encut
        kgrid = settings.kgrid
        supplied_kmesh = settings.supplied_kmesh
        # Set original encut and kgrid as we will overwrite these during the
        # course of the convergence tests.
        self.ctx.converge.settings.encut_org = copy.deepcopy(encut)
        self.ctx.converge.settings.kgrid_org = copy.deepcopy(kgrid)
        encut_org = self.ctx.converge.settings.encut_org
        encut_start = settings.encut_start
        encut_step = settings.encut_step
        encut_samples = settings.encut_samples
        rec_cell = self.ctx.converge.kpoints.cell
        k_spacing = self.ctx.converge.settings.k_spacing
        self.ctx.converge.pw_data = []
        self.ctx.converge.encut_sampling = None
        self.ctx.converge.pw_iteration = 0
        self.ctx.converge.run_pw_conv_calcs = False
        self.ctx.converge.pw_workchains = []

        # Detect what kind of convergence tests that needs to be run.
        if encut_org is None:
            # No encut supplied, run plane wave convergence tests.
            if not supplied_kmesh:
                # Set sensible k-point grid (k_spacing stepping in zone)
                # for plane wave cutoff tests.
                kgrid = fetch_k_grid(rec_cell, k_spacing)
                self.ctx.converge.settings.kgrid = kgrid
                # Update grid.
                kpoints = get_data_class('array.kpoints')
                self.ctx.converge.kpoints = kpoints(kpoints_mesh=kgrid, cell_from_structure=self.ctx.converge.structure)
            # Turn on plane wave convergene tests.
            self.ctx.converge.run_pw_conv_calcs = True
            # make encut test vector
            self.ctx.converge.encut_sampling = [encut_start + x * encut_step for x in range(encut_samples)]

        return

    def _init_kpoints_conv(self):
        """
        Initialize the kpoints convergence tests.

        Parameters
        ----------
        kpoints_iteration : int
            The current iteration for the kpoints convergence calculations.
        kpoints_conv_steps_reached : bool
            True if the last kpoints step has been reached.
        """
        settings = self.ctx.converge.settings
        kgrid_org = settings.kgrid_org
        supplied_kmesh = settings.supplied_kmesh
        self.ctx.converge.k_data = None
        self.ctx.converge.k_sampling = None
        self.ctx.converge.kpoints_iteration = 0
        self.ctx.converge.run_kpoints_conv_calcs = False
        self.ctx.converge.kpoints_workchains = []
        if kgrid_org is None and not supplied_kmesh:
            self.ctx.converge.k_data = []
            # No kpoint grid supplied, run kpoints convergence tests.
            self.ctx.converge.run_kpoints_conv_calcs = True

            # Make kpoint test vectors.
            # Usually one expect acceptable convergence with a
            # step size of 0.1/AA, typically:
            # 8 AA lattice vector needs roughly 8 kpoints.
            # 4 AA lattice vector needs roughly 16 kpoints etc.
            # Start convergence test with a step size of 0.5/AA,
            # round values up.
            self.ctx.converge.k_sampling = [0.6, 0.5, 0.4]
            # k_step = settings.k_step
            # k_samples = settings.k_samples
            # self.ctx.converge.k_sampling = [x * k_step for x in
            # range(k_samples, 0, -1)]

        return

    def init_converged(self):
        """Prepare to run the final calculation."""
        # Structure should be the same as the initial.
        self.ctx.inputs.structure = self.inputs.structure
        # The plane wave cutoff needs to be updated in the INCAR to the set
        # value.
        converged_incar_dict = self.inputs.incar.get_dict()
        converged_incar_dict.update({'encut': self.ctx.converge.settings.encut})
        self.ctx.inputs.incar = get_data_node('parameter', dict=converged_incar_dict)
        # And finally, the k-point grid needs to be updated to the set value, but
        # only if a kpoint mesh was not supplied
        if not self.ctx.converge.settings.supplied_kmesh:
            kpoints = get_data_class('array.kpoints')
            self.ctx.inputs.kpoints = kpoints(kpoints_mesh=self.ctx.converge.settings.kgrid, cell_from_structure=self.ctx.inputs.structure)

        self.ctx.running_kpoints = False
        self.ctx.running_pw = False
        self.ctx.replace_nodes = False

        # reset all the relaxation flag
        self.ctx.inputs.relax = self.inputs.relax

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

    def _replace_nodes(self):
        """Replaces the ctx.input nodes from the previous calculations."""
        self.ctx.inputs.structure = copy_structure(self.ctx.converge.structure)
        self.ctx.inputs.incar = copy_parameter(self.ctx.converge.incar)
        # Only the k-points if no mesh was supplied
        if not self.ctx.converge.settings.supplied_kmesh:
            self.ctx.inputs.kpoints = copy_kpoints(self.ctx.converge.kpoints, self.ctx.inputs.structure)

        return

    def init_next_workchain(self):
        """Initialize the next workchain calculation."""

        try:
            self.ctx.inputs
        except AttributeError:
            raise ValueError('no input dictionary was defined in self.ctx.inputs')

        if self.ctx.replace_nodes:
            self._replace_nodes()

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

        return self.ctx.converge.run_pw_conv_calcs and self.ctx.converge.settings.displace

    def run_pw_conv_comp_calcs(self):
        """Should a new plane wave cutoff compression convergence calculation run?"""

        return self.ctx.converge.run_pw_conv_calcs and self.ctx.converge.settings.compress

    def run_kpoints_conv_calcs(self):
        """Should a new kpoints convergence calculation run?"""
        return self.ctx.converge.run_kpoints_conv_calcs

    def run_kpoints_conv_disp_calcs(self):
        """Should a new kpoints displacement convergence calculation run?"""

        return self.ctx.converge.run_kpoints_conv_calcs and self.ctx.converge.settings.displace

    def run_kpoints_conv_comp_calcs(self):
        """Should a new kpoints compression convergence calculation run?"""

        return self.ctx.converge.run_kpoints_conv_calcs and self.ctx.converge.settings.compress

    def init_pw_conv_calc(self):
        """Initialize a single plane wave convergence calculation."""

        # Update the plane wave cutoff
        encut = self.ctx.converge.encut_sampling[self.ctx.converge.pw_iteration]
        self.ctx.converge.settings.encut = encut
        incar_dict = self.ctx.converge.incar.get_dict()
        incar_dict.update({'encut': self.ctx.converge.settings.encut})
        self.ctx.converge.incar = get_data_node('parameter', dict=incar_dict)
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
            self.ctx.pw_workchains[-1]
        except IndexError:
            self.report('the plane wave convergence calculation finished ' 'without returning a {}'.format(self._next_workchain.__name__))

        encut = self.ctx.converge.settings.encut
        if not exit_status:
            # fetch total energy
            energy = 0.0
            # energy =
            # workchain.out.output_energies.get_array('energy_no_entropy')

            # fetch force
            # forces = workchain.out.output_final_forces.get_array('forces')

            # locate maximum force
            # max_force = np.amax(forces)
            max_force = 0.0

            # fetch bands and occupations
            # bands = workchain.out.output_bands

            # And VBM and band gap (currently only do band gap
            # in order to use internal Aiida function)
            max_valence_band = 0.0
            #_, gap = find_bandgap(bands)
            gap = 0.0

            # add stuff to the converge context
            self.ctx.converge.pw_data.append([encut, energy, max_force, max_valence_band, gap])
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
        kpoints = get_data_class('array.kpoints')
        self.ctx.converge.kpoints = kpoints(kpoints_mesh=kgrid, cell_from_structure=self.ctx.inputs.structure)
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
            self.ctx.kpoints_workchains[-1]
        except IndexError:
            self.report('the k-point grid convergence calculation ' 'finished without returning a {}'.format(self._next_workchain.__name__))

        kgrid = self.ctx.converge.settings.kgrid
        encut = self.ctx.converge.settings.encut
        if not exit_status:
            # fetch total energy
            energy = 0.0
            # energy =
            # workchain.out.output_energies.get_array('energy_no_entropy')

            # fetch force
            # forces = workchain.out.output_final_forces.get_array('forces')

            # locate maximum force
            # max_force = np.amax(forces)
            max_force = 0.0

            # fetch bands and occupations
            # bands = workchain.out.output_bands

            # find VBM and band gap (currently only do band gap
            # in order to use internal Aiida function)
            max_valence_band = 0.0
            #_, gap = find_bandgap(bands)
            gap = 0.0
            # add stuff to the converge context
            self.ctx.converge.k_data.append([kgrid[0], kgrid[1], kgrid[2], encut, energy, max_force, max_valence_band, gap])
        else:
            self.ctx.converge.k_data.append([kgrid[0], kgrid[1], kgrid[2], encut, None, None, None, None])

        return

    def analyze_pw_after_comp(self):
        """Return True if we are running compressed convergence tests."""
        return self.ctx.converge.settings.compress

    def analyze_pw_after_disp(self):
        """Return True if we are running displaced convergence tests."""
        return self.ctx.converge.settings.displace

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
        displace = settings.displace
        compress = settings.compress
        cutoff_type = settings.cutoff_type
        cutoff_value = settings.cutoff_value

        # Notify the user
        if self._verbose:
            self.report('||||||||||||||||||||||||||||||||' 'All convergence tests are done.' '||||||||||||||||||||||||||||||||')

        # We know the recommended plane wave cutoff have been updated if no displacement
        # or compression tests have been performed, otherwise refetch
        # We have not yet set a recommendation for the k-point grid, need to fetch
        # the recommendation fromt the original convergence test
        if displace or compress:
            pw_data_org = self.ctx.converge.pw_data_org
            if pw_data_org is not None:
                encut = self._check_pw_converged(pw_data_org, cutoff_type, cutoff_value)
            else:
                encut = self.ctx.converge.settings.encut_org
            k_data_org = self.ctx.converge.k_data_org
            if not settings.supplied_kmesh:
                kgrid = self._check_kpoints_converged(k_data_org, cutoff_type, cutoff_value)
                if self._verbose:
                    self.report('The original convergence test suggest to use a plane wave cutoff '
                                'of {encut} eV and a k-point grid of {kgrid0}x{kgrid1}x{kgrid2} '
                                'for the convergence criteria {cutoff_type} and a cutoff of '
                                '{cutoff_value}. '
                                '||||||||||||||||||||||||||||||||'.format(
                                    encut=encut,
                                    kgrid0=kgrid[0],
                                    kgrid1=kgrid[1],
                                    kgrid2=kgrid[2],
                                    cutoff_type=cutoff_type,
                                    cutoff_value=cutoff_value))
            else:
                kgrid = None
                if self._verbose:
                    self.report('The original convergence test suggest to use a plane wave cutoff '
                                'of {encut} eV for the convergence criteria {cutoff_type} and a '
                                'cutoff of {cutoff_value}. The user supplied a k-point mesh.'
                                '||||||||||||||||||||||||||||||||'.format(encut=encut, cutoff_type=cutoff_type, cutoff_value=cutoff_value))

        if displace:
            if not compress:
                # We have data sitting from the displacement tests
                self.ctx.converge.pw_data_displacement = copy.deepcopy(self.ctx.converge.pw_data)
                self.ctx.converge.k_data_displacement = copy.deepcopy(self.ctx.converge.k_data)
            encut_diff_displacement, kgrid_diff_displacement = self._analyze_conv_disp(pw_data_org, k_data_org)
            self._set_encut_and_kgrid(encut_diff_displacement, kgrid_diff_displacement)

        if compress:
            # We have data sitting from the compression tests
            self.ctx.converge.pw_data_comp = copy.deepcopy(self.ctx.converge.pw_data)
            self.ctx.converge.k_data_comp = copy.deepcopy(self.ctx.converge.k_data)
            encut_diff_comp, kgrid_diff_comp = self._analyze_conv_comp(pw_data_org, k_data_org)
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
            self.report('We were not able to obtain a convergence of the plane wave cutoff'
                        'to the specified cutoff. This could also be caused by failures of '
                        'the calculations producing results for the convergence tests. Setting '
                        'the plane wave cutoff to the highest specified value: {encut} eV'.format(encut=self.ctx.converge.pw_data[-1][0]))
            settings.encut = self.ctx.converge.pw_data[-1][0:3]
        if not settings.supplied_kmesh and self.ctx.converge.settings.kgrid is None:
            self.report('We were not able to obtain a convergence of the k-point grid'
                        'to the specified cutoff. This could also be caused by failures of '
                        'the calculations producing results for the convergence tests. Setting '
                        'the k-point grid sampling to the highest specified value: {kgrid}'.format(kgrid=self.ctx.converge.k_data[-1][0]))
            settings.kgrid = self.ctx.converge.k_data[-1][0:3]

        return

    def _analyze_conv(self):
        """
        Analyze convergence using no displacements or compression.

        Note that, in the case of no displacements or compressions, the
        converged plane wave cutoff is already stored.
        """

        settings = self.ctx.converge.settings
        cutoff_type = settings.cutoff_type
        cutoff_value = settings.cutoff_value

        encut = settings.encut
        k_data = self.ctx.converge.k_data
        if not settings.supplied_kmesh:
            kgrid = self._check_kpoints_converged(k_data, cutoff_type, cutoff_value)
            if self._verbose:
                self.report('No atomic displacements or compression were performed. '
                            '||||||||||||||||||||||||||||||||'
                            'The convergence test suggest to use a plane wave cutoff of '
                            '{encut} eV and a k-point grid {kgrid0}x{kgrid1}x{kgrid2} for '
                            'the convergence criteria {cutoff_type} and a cutoff of {cutoff_value}'
                            '||||||||||||||||||||||||||||||||'.format(
                                encut=encut,
                                kgrid0=kgrid[0],
                                kgrid1=kgrid[1],
                                kgrid2=kgrid[2],
                                cutoff_type=cutoff_type,
                                cutoff_value=cutoff_value))
        else:
            kgrid = None
            if self._verbose:
                self.report('No atomic displacements or compression were performed. '
                            '||||||||||||||||||||||||||||||||'
                            'The convergence test suggest to use a plane wave cutoff of '
                            '{encut} eV for the convergence criteria {cutoff_type} and a '
                            'cutoff of {cutoff_value}. The user supplied a k-point mesh.'
                            '||||||||||||||||||||||||||||||||'.format(encut=encut, cutoff_type=cutoff_type, cutoff_value=cutoff_value))

        return encut, kgrid

    def _analyze_conv_disp_comp(self, encut_displacement, encut_comp, kgrid_displacement, kgrid_comp):
        """
        Analyze the convergence when both displacements and compression is performed.

        We take the maximum of the plane wave cutoff and the densest k-point grid as
        the recommended values.

        """

        settings = self.ctx.converge.settings
        cutoff_type = settings.cutoff_type
        cutoff_value = settings.cutoff_value_r
        # return the highest plane wave cutoff and densest grid (L2 norm)
        # of the two
        encut = max(encut_displacement, encut_comp)
        if not self.ctx.converge.settings.supplied_kmesh:
            if np.sqrt(sum([x**2 for x in kgrid_displacement])) > np.sqrt(sum([x**2 for x in kgrid_comp])):
                kgrid = kgrid_displacement
            else:
                kgrid = kgrid_comp
            if self._verbose:
                self.report('||||||||||||||||||||||||||||||||'
                            'The convergence tests, taking the highest required plane-wave and '
                            'k-point values for both the atomic displacement and compression '
                            'tests suggests to use a plane wave cutoff of '
                            '{encut} eV and a k-point grid of '
                            '{kgrid0}x{kgrid1}x{kgrid2} for the convergence criteria '
                            '{cutoff_type} and a cutoff of {cutoff_value} '
                            '||||||||||||||||||||||||||||||||'.format(
                                encut=encut,
                                kgrid0=kgrid[0],
                                kgrid1=kgrid[1],
                                kgrid2=kgrid[2],
                                cutoff_type=cutoff_type,
                                cutoff_value=cutoff_value))
        else:
            kgrid = None
            if self._verbose:
                self.report('||||||||||||||||||||||||||||||||'
                            'The convergence tests, taking the highest required plane-wave '
                            'for both the atomic displacement and compression '
                            'tests suggests to use a plane wave cutoff of '
                            '{encut} eV for the convergence criteria '
                            '{cutoff_type} and a cutoff of {cutoff_value}. The user supplied '
                            'a k-point mesh.'
                            '||||||||||||||||||||||||||||||||'.format(encut=encut, cutoff_type=cutoff_type, cutoff_value=cutoff_value))

        return encut, kgrid

    def _analyze_conv_disp(self, pw_data_org, k_data_org):
        """Analyze the convergence when atomic displacements are performed."""
        settings = self.ctx.converge.settings
        encut_org = settings.encut_org
        kgrid_org = settings.kgrid_org
        cutoff_type = settings.cutoff_type
        cutoff_value = settings.cutoff_value
        cutoff_value_r = settings.cutoff_value_r
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
        if not settings.supplied_kmesh:
            if self._verbose:
                self.report('Performed atomic displacements. '
                            '||||||||||||||||||||||||||||||||'
                            'The convergence test using the difference between the original '
                            'and displaced dataset suggest to use a plane wave cutoff of '
                            '{encut_diff_displacement} eV and a k-point grid of '
                            '{kgrid_diff_displacement0}x{kgrid_diff_displacement1}x'
                            '{kgrid_diff_displacement2} '
                            'for the convergence criteria {cutoff_type} and a cutoff of '
                            '{cutoff_value_r}.'
                            '(The displacement convergence test suggest to use a plane wave cutoff '
                            '{encut_displacement} eV and a k-point grid of '
                            '{kgrids0}x{kgrids1}x{kgrids2} '
                            'for the convergence criteria {cutoff_type} and a cutoff of {cutoff_value}.)'
                            '||||||||||||||||||||||||||||||||'.format(
                                encut_diff_displacement=encut_diff_displacement,
                                kgrid_diff_displacement0=kgrid_diff_displacement[0],
                                kgrid_diff_displacement1=kgrid_diff_displacement[1],
                                kgrid_diff_displacement2=kgrid_diff_displacement[2],
                                cutoff_type=cutoff_type,
                                cutoff_value_r=cutoff_value_r,
                                cutoff_value=cutoff_value,
                                encut_displacement=encut_displacement,
                                kgrids0=kgrid_displacement[0],
                                kgrids1=kgrid_displacement[1],
                                kgrids2=kgrid_displacement[2]))
        else:
            if self._verbose:
                self.report('Performed atomic displacements. '
                            '||||||||||||||||||||||||||||||||'
                            'The convergence test using the difference between the original '
                            'and displaced dataset suggest to use a plane wave cutoff of '
                            '{encut_diff_displacement} eV for the convergence criteria '
                            '{cutoff_type} and a cutoff of {cutoff_value_r}. The user supplied '
                            'a k-point grid.'
                            '(The displacement convergence test suggest to use a plane wave cutoff '
                            '{encut_displacement} eV using a cutoff of {cutoff_value}.)'
                            '||||||||||||||||||||||||||||||||'.format(
                                encut_diff_displacement=encut_diff_displacement,
                                cutoff_type=cutoff_type,
                                cutoff_value_r=cutoff_value_r,
                                cutoff_value=cutoff_value,
                                encut_displacement=encut_displacement))

        return encut_diff_displacement, kgrid_diff_displacement

    def _analyze_conv_comp(self, pw_data_org, k_data_org):
        """Analize the relative convergence due to unit cell compression."""

        settings = self.ctx.converge.settings
        encut_org = settings.encut_org
        kgrid_org = settings.kgrid_org
        cutoff_type = settings.cutoff_type
        cutoff_value = settings.cutoff_value
        cutoff_value_r = settings.cutoff_value_r
        # Do not be too confused here.
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
        if not settings.supplied_kmesh:
            if self._verbose:
                self.report('Performed compression. '
                            '||||||||||||||||||||||||||||||||'
                            'The convergence test using the difference between the '
                            'original and dataset with a volume change suggest to use a '
                            'plane wave cutoff of {encut_diff_comp} eV and a k-point grid '
                            'of {kgrid_diff_comp0}x{kgrid_diff_comp1}x{kgrid_diff_comp2} '
                            'for the convergence criteria {cutoff_type} and a cutoff of '
                            '{cutoff_value_r}.'
                            '(The convergence test with a volume change suggest to use a plane '
                            'wave cutoff of {encut_comp} eV and a k-point grid of '
                            '{kgrid_comp0}x{kgrid_comp1}x{kgrid_comp2} for the convergence criteria '
                            'using a cutoff of {cutoff_value}.) '
                            '||||||||||||||||||||||||||||||||'.format(
                                encut_diff_comp=encut_diff_comp,
                                kgrid_diff_comp0=kgrid_diff_comp[0],
                                kgrid_diff_comp1=kgrid_diff_comp[1],
                                kgrid_diff_comp2=kgrid_diff_comp[2],
                                cutoff_type=cutoff_type,
                                cutoff_value=cutoff_value,
                                cutoff_value_r=cutoff_value_r,
                                encut_comp=encut_comp,
                                kgrid_comp0=kgrid_comp[0],
                                kgrid_comp1=kgrid_comp[1],
                                kgrid_comp2=kgrid_comp[2]))
        else:
            if self._verbose:
                self.report('Performed compression. '
                            '||||||||||||||||||||||||||||||||'
                            'The convergence test using the difference between the '
                            'original and dataset with a volume change suggest to use a '
                            'plane wave cutoff of {encut_diff_comp} eV '
                            'for the convergence criteria {cutoff_type} and a cutoff of '
                            '{cutoff_value_r}. The user supplied a k-point mesh.'
                            '(The convergence test with a volume change suggest to use a plane '
                            'wave cutoff of {encut_comp} eV using a cutoff of {cutoff_value}.) '
                            '||||||||||||||||||||||||||||||||'.format(
                                encut_diff_comp=encut_diff_comp,
                                cutoff_type=cutoff_type,
                                cutoff_value=cutoff_value,
                                cutoff_value_r=cutoff_value_r,
                                encut_comp=encut_comp))

        return encut_diff_comp, kgrid_diff_comp

    def store_conv(self):
        """Store the obtained convergence data on nodes."""
        convergence = get_data_class('array')()

        if self._verbose:
            self.report("attaching the node {}<{}> as '{}'".format(convergence.__class__.__name__, convergence.pk,
                                                                   'output_convergence_data'))

        # Store regular conversion data
        try:
            convergence.set_array('pw_regular', np.array(self.ctx.converge.pw_data_org))
        except AttributeError:
            convergence.set_array('pw_regular', np.array(self.ctx.converge.pw_data))
        try:
            convergence.set_array('kpoints_regular', np.array(self.ctx.converge.k_data_org))
        except AttributeError:
            convergence.set_array('kpoints_regular', np.array(self.ctx.converge.k_data))

        # Then possibly displacement
        try:
            convergence.set_array('pw_displacement', np.array(self.ctx.converge.pw_data_displacement))
            convergence.set_array('kpoints_displacement', np.array(self.ctx.converge.k_data_displacement))
        except AttributeError:
            pass

        # And finally for compression
        try:
            convergence.set_array('pw_compression', np.array(self.ctx.converge.pw_data_comp))
            convergence.set_array('kpoints_compression', np.array(self.ctx.converge.k_data_comp))
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
            cutoff_type = self.ctx.converge.settings.cutoff_type
        if cutoff_value is None:
            cutoff_value = self.ctx.converge.settings.cutoff_value

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
                if self._verbose:
                    self.report('A plane wave cutoff of {encut} eV is considered converged '
                                'by observing the convergence of the {cutoff_type} to within a '
                                'difference of {cutoff_value}.'.format(
                                    encut=pw_data[encut][0], cutoff_type=cutoff_type, cutoff_value=cutoff_value))
                encut_okey = True
                index = encut
                break
        if not encut_okey:
            if self._verbose:
                self.report('Could not obtain convergence for {cutoff_type} with a cutoff '
                            'parameter of {cutoff_value}'.format(cutoff_type=cutoff_type, cutoff_value=cutoff_value))
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
            cutoff_type = self.ctx.converge.settings.cutoff_type
        if cutoff_value is None:
            cutoff_value = self.ctx.converge.settings.cutoff_value

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
                self.report('The k-point grid of {k_data0}x{k_data1}x{k_data2} is '
                            'considered converged by observing the convergence of '
                            'the {cutoff_type} to within a difference '
                            'of {cutoff_value}'.format(
                                k_data0=k_data[k][0],
                                k_data1=k_data[k][1],
                                k_data2=k_data[k][2],
                                cutoff_type=cutoff_type,
                                cutoff_value=cutoff_value))
                k_cut_okey = True
                index = k
                break
        if not k_cut_okey:
            self.report('Could not find a dense enough grid to obtain a {cutoff_type} '
                        'cutoff of {cutoff_value})'.format(cutoff_type=cutoff_type, cutoff_value=cutoff_value))
            return None

        return k_data[index][0:3]

    def verify_next_workchain(self):
        """Verify and inherit exit status from child workchains."""

        # Adopt exit status from last child workchain (supposed to be
        # successfull)
        next_workchain_exit_status = self.ctx.workchains[-1].exit_status
        if not next_workchain_exit_status:
            self.exit_status = 0
            return
        self.exit_status = next_workchain_exit_status
        self.report('The child {} returned a non-zero exit status, {} inherits exit status {}'.format(
            self._next_workchain, self.__class__.__name__, next_workchain_exit_status))
        return

    def results(self):
        """Attach the outputs specified in the output specification from the last completed calculation."""

        if not self.exit_status:
            self.report('{} completed'.format(self.__class__.__name__))

            workchain = self.ctx.workchains[-1]

            for name, _ in self.spec().outputs.iteritems():
                # if port.required and (name not in workchain.out or name not in self.out):
                #    self.report('the spec specifying the output {} as required '
                #                'but was not an output of {}<{}> or already stored '
                #                'in the output of this workchain'.
                #                format(name, self._next_workchain.__name__,
                #                       workchain.pk))
                if name in workchain.out:
                    node = workchain.out[name]
                    self.out(name, workchain.out[name])
                    if self._verbose:
                        self.report("attaching the node {}<{}> as '{}'".format(node.__class__.__name__, node.pk, name))

        return

    def finalize(self):
        """Finalize the workchain."""

        return self.exit_status

    def run_conv_calcs(self):
        """Determines if convergence calcs are to be run at all."""
        return self.run_kpoints_conv_calcs() or self.run_pw_conv_calcs()

    def _displace_structure(self):
        """Displace the input structure according to the supplied settings."""

        displacement_vector = self.ctx.converge.settings.displacement_vector
        displacement_distance = self.ctx.converge.settings.displacement_distance
        displacement_atom = self.ctx.converge.settings.displacement_atom
        # Set displacement
        displacement = displacement_distance * displacement_vector

        # Displace and return new structure
        return displaced_structure(self.ctx.inputs.structure, displacement, displacement_atom)

    def _compress_structure(self):
        """Compress the input structure according to the supplied settings."""

        volume_change = self.ctx.converge.settings.volume_change
        # Apply compression and tension
        comp_structure = compressed_structure(self.ctx.inputs.structure, volume_change)
        # Make sure we also reset the reciprocal cell
        kpoints = get_data_class('array.kpoints')
        self.ctx.converge.kpoints = kpoints(kpoints_mesh=[1, 1, 1], cell_from_structure=comp_structure)

        return comp_structure
