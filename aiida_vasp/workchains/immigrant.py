"""
Vasp immigrant workchain.

-------------------------
Workchain to import files of VASP calculation in a folder.
"""

# pylint: disable=attribute-defined-outside-init
from aiida.engine import while_
from aiida.plugins import WorkflowFactory, CalculationFactory
from aiida.orm import Code
from aiida.engine.processes.workchains.restart import BaseRestartWorkChain
from aiida_vasp.utils.aiida_utils import get_data_class, get_data_node


class VaspImmigrantWorkChain(BaseRestartWorkChain):
    """Import files of VASP calculation in the specified folder path."""

    _verbose = False
    _process_class = CalculationFactory('vasp.immigrant')
    _next_workchain_string = 'vasp.vasp'
    _next_workchain = WorkflowFactory(_next_workchain_string)

    @classmethod
    def define(cls, spec):
        super().define(spec)
        spec.input('code', valid_type=Code, required=True)
        spec.input('folder_path', valid_type=get_data_class('str'), required=True)
        spec.input('settings', valid_type=get_data_class('dict'), required=False)
        spec.input('options', valid_type=get_data_class('dict'), required=False)
        spec.input('potential_family', valid_type=get_data_class('str'), required=False)
        spec.input('potential_mapping', valid_type=get_data_class('dict'), required=False)
        spec.input('use_chgcar',
                   valid_type=get_data_class('bool'),
                   required=False,
                   default=lambda: get_data_node('bool', False),
                   help="""
            If True, WavefunData (of WAVECAR) is attached.
            """)
        spec.input('use_wavecar',
                   valid_type=get_data_class('bool'),
                   required=False,
                   default=lambda: get_data_node('bool', False),
                   help="""
            If True, WavefunData (of WAVECAR) is attached.
            """)
        spec.input('max_iterations',
                   valid_type=get_data_class('int'),
                   required=False,
                   default=lambda: get_data_node('int', 1),
                   help="""
            The maximum number of iterations to perform.
            """)
        spec.input('clean_workdir',
                   valid_type=get_data_class('bool'),
                   required=False,
                   default=lambda: get_data_node('bool', False),
                   help="""
            If True, clean the work dir upon the completion of a successfull calculation.
            """)
        spec.input('verbose',
                   valid_type=get_data_class('bool'),
                   required=False,
                   default=lambda: get_data_node('bool', False),
                   help="""
            If True, enable more detailed output during workchain execution.
            """)
        spec.exit_code(0, 'NO_ERROR', message='the sun is shining')
        spec.outline(
            cls.setup,
            cls.init_inputs,
            while_(cls.should_run_process)(
                cls.run_process,
                cls.inspect_process,
            ),
            cls.results,
        )  # yapf: disable
        spec.expose_outputs(cls._next_workchain)

    def init_inputs(self):  # pylint: disable=too-many-branches, too-many-statements
        """Set inputs of VaspImmigrant calculation

        Initial inputs (self.ctx.inputs) AttributeDict is obtained from
        VaspImmigrant.get_inputs_from_folder(), where the following items
        are already set in respective AiiDA data types,

            code
            settings.import_from_path
            metadata['options']
            parameters
            structure
            kpoints
            potential (optional)

        """

        kwargs = self._get_kwargs()
        self.ctx.inputs = self._process_class.get_inputs_from_folder(self.inputs.code, self.inputs.folder_path.value, **kwargs)
        if 'options' in self.inputs:
            self.ctx.inputs.metadata.options.update(self.inputs.options)
        if 'metadata' in self.inputs:
            label = self.inputs.metadata.get('label', '')
            description = self.inputs.metadata.get('description', '')
            self.ctx.inputs.metadata['label'] = label
            self.ctx.inputs.metadata['description'] = description
        return self.exit_codes.NO_ERROR  # pylint: disable=no-member

    def _get_kwargs(self):
        """kwargs dictionary for VaspImmigrant calculation is created."""
        kwargs = {'use_chgcar': False, 'use_wavecar': False}
        for key in ('use_chgcar', 'use_wavecar', 'settings', 'potential_family', 'potential_mapping'):
            if key in self.inputs:
                if key in ('settings', 'potential_mapping'):
                    kwargs[key] = self.inputs[key].get_dict()
                else:
                    kwargs[key] = self.inputs[key].value
        return kwargs
