# This is an adapted copy of the quantum espresso base workchain.

from copy import deepcopy
from aiida.orm import Code
from aiida.orm.data.base import Bool, Int, Str
from aiida.orm.data.folder import FolderData
from aiida.orm.data.remote import RemoteData
from aiida.orm.data.parameter import ParameterData
from aiida.orm.data.structure import StructureData
from aiida.orm.data.array.bands import BandsData
from aiida.orm.data.array.kpoints import KpointsData
from aiida.orm.data.singlefile import SinglefileData
from aiida.orm.utils import CalculationFactory
from aiida.common.extendeddicts import AttributeDict
from aiida.work.run import submit
from aiida.work.workchain import ToContext, if_, while_

# imports of common workchain utils. these might be independent of the actual code
# and probably could be part of aiida-core?
from aiida_quantumespresso.common.exceptions import UnexpectedCalculationFailure
from aiida_quantumespresso.common.workchain.utils import ErrorHandlerReport
from aiida_quantumespresso.common.workchain.utils import register_error_handler
from aiida_quantumespresso.common.workchain.base.restart import BaseRestartWorkChain

# calculation specific utils, which have to be implemented for each plugin
from aiida_vasp.utils.pseudopotential import validate_and_prepare_pseudos_inputs


VASPCalculation = CalculationFactory('vasp.vasp')


class VASPBaseWorkChain(BaseRestartWorkChain):
    """
    Base workchain to launch a VASP calculation
    """
    _verbose = True
    _calculation_class = VASPCalculation
    _error_handler_entry_point = 'aiida_quantumespresso.workflow_error_handlers.pw.base'

    def __init__(self, *args, **kwargs):
        super(VASPBaseWorkChain, self).__init__(*args, **kwargs)

        
    @classmethod
    def define(cls, spec):
        super(VASPBaseWorkChain, cls).define(spec)
        spec.input('code', valid_type=Code)
        spec.input('structure', valid_type=StructureData)
        spec.input('kpoints', valid_type=KpointsData)
        spec.input('parameters', valid_type=ParameterData)
        spec.input_group('pseudos', required=False)
        spec.input('pseudo_family', valid_type=Str, required=False)
        spec.input('parent_folder', valid_type=RemoteData, required=False)
        spec.input('vdw_table', valid_type=SinglefileData, required=False)
        spec.input('settings', valid_type=ParameterData, required=False)
        spec.input('options', valid_type=ParameterData, required=False)
        spec.input('automatic_parallelization', valid_type=ParameterData, required=False)
        spec.outline(
            cls.setup,
            cls.validate_inputs,
            while_(cls.should_run_calculation)(
                cls.prepare_calculation,
                cls.run_calculation,
                cls.inspect_calculation,
            ),
            cls.results,
        )
        spec.output('output_band', valid_type=BandsData, required=False)
        spec.output('output_structure', valid_type=StructureData, required=False)
        spec.output('output_parameters', valid_type=ParameterData)
        spec.output('remote_folder', valid_type=RemoteData)
        spec.output('retrieved', valid_type=FolderData)

    def validate_inputs(self):
        """
        Define context dictionary 'inputs_raw' with the inputs for the PwCalculations as they were at the beginning
        of the workchain. Changes have to be made to a deep copy so this remains unchanged and we can always reset
        the inputs to their initial state. Inputs that are not required by the workchain will be given a default value
        if not specified or be validated otherwise.
        """
        self.ctx.inputs_raw = AttributeDict({
            'code': self.inputs.code,
            'structure': self.inputs.structure,
            'kpoints': self.inputs.kpoints,
            'parameters': self.inputs.parameters.get_dict()
        })
        

        if 'parent_folder' in self.inputs:
            self.ctx.inputs_raw.parent_folder = self.inputs.parent_folder
        else:
            self.ctx.inputs_raw.parent_folder = None

        if 'settings' in self.inputs:
            self.ctx.inputs_raw.settings = self.inputs.settings.get_dict()
        else:
            self.ctx.inputs_raw.settings = {}

        if 'options' in self.inputs:
            self.ctx.inputs_raw._options = self.inputs.options.get_dict()
        else:
            self.ctx.inputs_raw._options = {}

        if 'vdw_table' in self.inputs:
            self.ctx.inputs_raw.vdw_table = self.inputs.vdw_table

        # Options has to be specified
        if not any([key in self.inputs for key in ['options']]):
            self.abort_nowait('you have to specify the options input')
            return

        # We better make sure that the options satisfy minimum requirements
        num_machines = self.ctx.inputs_raw['_options'].get('resources', {}).get('num_machines', None)
        max_wallclock_seconds = self.ctx.inputs_raw['_options'].get('max_wallclock_seconds', None)

        if num_machines is None or max_wallclock_seconds is None:
            self.abort_nowait("no automatic_parallelization requested, but the options do not specify both '{}' and '{}'"
                .format('num_machines', 'max_wallclock_seconds'))

        # Validate the inputs related to pseudopotentials
        structure = self.inputs.structure
        pseudos = self.inputs.get('pseudos', None)
        pseudo_family = self.inputs.get('pseudo_family', None)

        try:
            self.ctx.inputs_raw.paw = validate_and_prepare_pseudos_inputs(structure, pseudos, pseudo_family)
        except ValueError as exception:
            self.abort_nowait('{}'.format(exception))

        # Assign a deepcopy to self.ctx.inputs which will be used by the BaseRestartWorkChain
        self.ctx.inputs = deepcopy(self.ctx.inputs_raw)


    def prepare_calculation(self):
        """
        Prepare the inputs for the next calculation
        """
        if self.ctx.restart_calc:
            self.ctx.inputs.parent_folder = self.ctx.restart_calc.out.remote_folder

    def _prepare_process_inputs(self, inputs):

        return super(VASPBaseWorkChain, self)._prepare_process_inputs(inputs)
