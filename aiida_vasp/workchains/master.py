# pylint: disable=attribute-defined-outside-init
"""
VerifyWorkChain.

Indented to be used to verify a calculation, perform corrections in inputs files and
restart depending on physical principles etc. E.g. issues that are outside the Calculators awereness,
or not currently checked in it. This workchain does currently nothing.
"""
from aiida.common.extendeddicts import AttributeDict
from aiida.engine.workchain import WorkChain, if_, append_
from aiida.orm import WorkflowFactory
from aiida_vasp.utils.aiida_utils import get_data_class, get_data_node
from aiida_vasp.utils.workchains import prepare_process_inputs


class MasterWorkChain(WorkChain):
    """The master workchain that selects sub workchains to perform necessary calculations."""

    _verbose = False
    _base_workchains_string = 'vasp.converge'
    _bands_workchain_string = 'vasp.bands'
    _base_workchains = WorkflowFactory(_base_workchains_string)
    _bands_workchain = WorkflowFactory(_bands_workchain_string)

    @classmethod
    def define(cls, spec):
        super(MasterWorkChain, cls).define(spec)
        spec.expose_inputs(cls._base_workchains, exclude=['extract_bands', 'settings'])
        spec.input(
            'extract_bands',
            valid_type=get_data_class('bool'),
            required=False,
            default=get_data_node('bool', False),
            help="""
            Do you want to extract the band structure?
            """)
        spec.input(
            'relax',
            valid_type=get_data_class('bool'),
            required=False,
            default=get_data_node('bool', False),
            help="""
            Do you want to relax the structure?
            """)
        spec.outline(
            cls.initialize,
            cls.init_workchain,
            cls.run_workchain,
            cls.verify_workchain,
            if_(cls.extract_bands)(
                cls.init_bands,
                cls.init_workchain,
                cls.run_workchain,
                cls.verify_workchain
            ),
            cls.finalize
        )  # yapf: disable

        # spec.output('output_parameters', valid_type=get_data_class('parameter'))
        # spec.output('remote_folder', valid_type=get_data_class('remote'))
        # spec.output('retrieved', valid_type=get_data_class('folder'))
        # spec.output('output_structure', valid_type=get_data_class('structure'), required=False)
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

        spec.expose_outputs(cls._bands_workchain)

    def initialize(self):
        """Initialize."""
        self._init_context()
        self._init_inputs()
        self._init_settings()
        self._set_base_workchain()

        return

    def _init_context(self):
        """Initialize context variables."""
        self.ctx.inputs = AttributeDict()

        return

    def _init_inputs(self):
        """Initialize inputs."""
        try:
            self._verbose = self.inputs.verbose.value
        except AttributeError:
            pass

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

        return

    def _set_base_workchain(self):
        """Set the base workchain to be called."""
        self._next_workchain = self._base_workchains

    def init_bands(self):
        """Initialize the run to extract the band structure."""
        self._next_workchain = self._bands_workchain
        # Make sure the charge density is added from the previous run
        self.ctx.inputs.chgcar = self.ctx.workchains[-1].out['output_chgcar']
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
        self.ctx.inputs.update(self.exposed_inputs(self._next_workchain))

        # Make sure we do not have any floating dict (convert to Dict)
        self.ctx.inputs = prepare_process_inputs(self.ctx.inputs)

    def run_workchain(self):
        """Run the workchain."""
        inputs = self.ctx.inputs
        running = self.submit(self._next_workchain, **inputs)

        if hasattr(running, 'pid'):
            self.report('launching {}<{}> '.format(self._next_workchain.__name__, running.pid))
        else:
            # Aiida < 1.0
            self.report('launching {}<{}> '.format(self._next_workchain.__name__, running.pk))

        return self.to_context(workchains=append_(running))

    def verify_workchain(self):
        """Inherit exit status from child workchains."""
        workchain = self.ctx.workchains[-1]
        # Adopt exit status from last child workchain (supposed to be
        # successfull)
        workchain_exit_status = workchain.exit_status
        if not workchain_exit_status:
            self.exit_status = 0
        else:
            self.exit_status = workchain_exit_status
            self.report('The child {}<{}> returned a non-zero exit status, {}<{}> '
                        'inherits exit status {}'.format(workchain.__class__.__name__, workchain.pk, self.__class__.__name__, self.pid,
                                                         workchain_exit_status))
        return

    def extract_bands(self):
        """Determines if we should extract the band structure."""
        return self.inputs.extract_bands.value

    def finalize(self):
        """Finalize the workchain."""

        if not self.exit_status:
            workchain = self.ctx.workchains[-1]
            self.out_many(self.exposed_outputs(workchain, self._next_workchain))
        return self.exit_status
