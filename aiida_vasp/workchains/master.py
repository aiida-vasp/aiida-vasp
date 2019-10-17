""" # noqa: D205
Verify workchain
----------------

Indented to be used to verify a calculation, perform corrections in inputs files and
restart depending on physical principles etc. E.g. issues that are outside the Calculators awereness,
or not currently checked in it. This workchain does currently nothing.
"""
# pylint: disable=attribute-defined-outside-init
from aiida.common.extendeddicts import AttributeDict
from aiida.engine import WorkChain, if_, append_
from aiida.plugins import WorkflowFactory
from aiida_vasp.utils.aiida_utils import get_data_class, get_data_node
from aiida_vasp.utils.workchains import prepare_process_inputs, compose_exit_code


class MasterWorkChain(WorkChain):
    """The master workchain that selects sub workchains to perform necessary calculations."""

    _verbose = False
    _base_workchains_string = 'vasp.converge'
    _bands_workchain_string = 'vasp.bands'
    _base_workchain = WorkflowFactory(_base_workchains_string)
    _bands_workchain = WorkflowFactory(_bands_workchain_string)

    @classmethod
    def define(cls, spec):
        super(MasterWorkChain, cls).define(spec)
        spec.expose_inputs(cls._base_workchain, exclude=['settings'])
        spec.input('settings', valid_type=get_data_class('dict'), required=False)
        spec.input('kpoints', valid_type=get_data_class('array.kpoints'), required=False)
        spec.input(
            'extract_bands',
            valid_type=get_data_class('bool'),
            required=False,
            default=get_data_node('bool', False),
            help="""
            Do you want to extract the band structure?
            """)
        spec.input(
            'extract_dos',
            valid_type=get_data_class('bool'),
            required=False,
            default=get_data_node('bool', False),
            help="""
            Do you want to extract the density of states?
            """)
        spec.input(
            'relax.perform',
            valid_type=get_data_class('bool'),
            required=False,
            default=get_data_node('bool', False),
            help="""
            Do you want to relax the structure?
            """)
        spec.input(
            'kpoints_distance',
            valid_type=get_data_class('float'),
            required=False,
            help="""
            The maximum distance between k-points in inverse AA.
            """)
        spec.outline(
            cls.initialize,
            cls.init_workchain,
            cls.run_next_workchain,
            cls.verify_next_workchain,
            if_(cls.extract_bands)(
                cls.init_bands,
                cls.init_workchain,
                cls.run_next_workchain,
                cls.verify_next_workchain
            ),
            cls.finalize
        )  # yapf: disable
        spec.expose_outputs(cls._bands_workchain)
        spec.exit_code(0, 'NO_ERROR', message='the sun is shining')
        spec.exit_code(420, 'ERROR_NO_CALLED_WORKCHAIN', message='no called workchain detected')
        spec.exit_code(500, 'ERROR_UNKNOWN', message='unknown error detected in the master workchain')

    def initialize(self):
        """Initialize."""
        self._init_context()
        self._init_inputs()
        self._init_settings()
        self._set_base_workchain()

    def _init_context(self):
        """Initialize context variables."""
        self.ctx.exit_code = self.exit_codes.ERROR_UNKNOWN  # pylint: disable=no-member
        self.ctx.inputs = AttributeDict()

    def _init_inputs(self):
        """Initialize inputs."""
        try:
            self._verbose = self.inputs.verbose.value
        except AttributeError:
            pass
        self._init_structure()
        self._init_kpoints()

    def _init_structure(self):
        """Initialize the structure."""
        self.ctx.inputs.structure = self.inputs.structure

    def _init_kpoints(self):
        """
        Initialize the kpoints.

        Since the reciprocal lattice is set in KpointsData, we need, if not performing
        convergence tests to access the single StructureData, update the
        supplied KpointsData with the right reciprocal lattice.

        """

        # First check if a k-point distance is supplied
        try:
            distance = self.inputs.kpoints_distance.value
        except AttributeError:
            # Then check if a KpointsData has been supplied
            try:
                kpoints = self.inputs.kpoints
            except AttributeError:
                # If neither, return, we need to run convergence tests
                return
            self.ctx.inputs.kpoints = kpoints
            return
        kpoints = get_data_class('array.kpoints')()
        kpoints.set_cell_from_structure(self.ctx.inputs.structure)
        kpoints.set_kpoints_mesh_from_density(distance)

    def _init_settings(self):
        """Initialize the settings."""
        # Make sure we parse the charge density from the previous run
        if 'settings' in self.inputs:
            settings = AttributeDict(self.inputs.settings.get_dict())
        else:
            settings = AttributeDict({'parser_settings': {}})
        if self.extract_bands():
            dict_entry = {'add_chgcar': True}
            try:
                settings.parser_settings.update(dict_entry)
            except AttributeError:
                settings.parser_settings = dict_entry
            dict_entry = {'ADDITIONAL_RETRIEVE_LIST': ['CHGCAR']}
            settings.update(dict_entry)
        self.ctx.inputs.settings = settings

    def _set_base_workchain(self):
        """Set the base workchain to be called."""
        self._next_workchain = self._base_workchain

    def init_bands(self):
        """Initialize the run to extract the band structure."""
        self._next_workchain = self._bands_workchain
        # Make sure the charge density is added from the previous run
        self.ctx.inputs.chgcar = self.ctx.workchains[-1].outputs['chgcar']
        # Remove parser extraction of the charge density file
        settings = self.ctx.inputs.settings.get_dict()
        try:
            del settings['ADDITIONAL_RETRIEVE_LIST']
        except AttributeError:
            pass
        try:
            del settings['parser_settings']['add_chgcar']
        except AttributeError:
            pass
        self.ctx.inputs.settings = settings

    def init_workchain(self):
        """Initialize the base workchain."""
        try:
            self.ctx.inputs
        except AttributeError:
            raise ValueError('No input dictionary was defined in self.ctx.inputs')

        # Add exposed inputs
        exclude = []
        if self.extract_bands:
            exclude = exclude + ['converge', 'relax', 'verify']
        self.ctx.inputs.update(self.exposed_inputs(self._next_workchain))

        # Make sure we do not have any floating dict (convert to Dict)
        self.ctx.inputs = prepare_process_inputs(self.ctx.inputs, namespaces=['relax', 'converge'])

    def run_next_workchain(self):
        """Run the next workchain."""
        inputs = self.ctx.inputs
        running = self.submit(self._next_workchain, **inputs)

        self.report('launching {}<{}> '.format(self._next_workchain.__name__, running.pk))

        return self.to_context(workchains=append_(running))

    def verify_next_workchain(self):
        """Inherit exit status from child workchains."""
        try:
            workchain = self.ctx.workchains[-1]
        except IndexError:
            self.report('There is no {} in the called workchain list.'.format(self._next_workchain.__name__))
            return self.exit_codes.ERROR_NO_CALLED_WORKCHAIN  # pylint: disable=no-member

        # Inherit exit status from last workchain (supposed to be
        # successfull)
        next_workchain_exit_status = workchain.exit_status
        next_workchain_exit_message = workchain.exit_message
        if not next_workchain_exit_status:
            self.ctx.exit_code = self.exit_codes.NO_ERROR  # pylint: disable=no-member
        else:
            self.ctx.exit_code = compose_exit_code(next_workchain_exit_status, next_workchain_exit_message)
            self.report('The called {}<{}> returned a non-zero exit status. '
                        'The exit status {} is inherited'.format(workchain.__class__.__name__, workchain.pk, self.ctx.exit_code))

        return self.ctx.exit_code

    def extract_bands(self):
        """Determines if we should extract the band structure."""
        return self.inputs.extract_bands.value

    def loop_structures(self):
        """Determine if we should continue to calculate structures."""
        return len(self.ctx.structures.get_list()) > 0

    def finalize(self):
        """Finalize the workchain."""

        workchain = self.ctx.workchains[-1]
        self.out_many(self.exposed_outputs(workchain, self._next_workchain))
