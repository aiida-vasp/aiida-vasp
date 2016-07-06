from helper import WorkflowHelper
from aiida.orm import Workflow
from aiida.common import aiidalogger
from aiida.common.utils import classproperty

logger = aiidalogger.getChild('ScfWorkflow')


class ScfWorkflow(Workflow):

    '''AiiDA workflow to run a VASP scf calculation
    the resulting WAVECAR and CHGCAR can then be reused in further workflows.'''
    Helper = WorkflowHelper

    def __init__(self, **kwargs):
        self.helper = self.Helper(parent=self)
        super(ScfWorkflow, self).__init__(**kwargs)

    def get_calc_maker(self):
        maker = self.helper._get_calc_maker('vasp.scf')
        maker.add_settings(icharg=0, istart=0)
        return maker

    @Workflow.step
    def start(self):
        params = self.get_parameters()
        self.append_to_report(self.helper._wf_start_msg())
        maker = self.get_calc_maker()
        calc = maker.new()
        calc.description = params.get('description', '')
        calc.store_all()
        calc.set_extras(params.get('extras'))
        self.attach_calculation(calc)
        self.append_to_report(
            self.helper._calc_start_msg('scf VASP run', calc))
        self.next(self.end)

    @Workflow.step
    def end(self):
        calc = self.helper._get_first_step_calc(self.start)
        output_links = ['charge_density', 'wavefunctions']
        valid_out = self.helper._verify_calc_output(calc, output_links)
        if valid_out:
            self.add_result('calc', calc)
            self.append_to_report(
                'Added the scf calculation as a result')
        else:
            self.append_to_report(
                self.helper._calc_invalid_outs_msg(calc, output_links))
        self.next(self.exit)

    def set_params(self, params, **kwargs):
        self.helper._verify_params(params)
        super(ScfWorkflow, self).set_params(params, **kwargs)

    @classmethod
    def get_params_template(cls):
        '''returns a dictionary with the necessary keys to
        run this workflow and explanations to each key as values'''
        return cls.Helper.get_params_template(continuation=False)

    @classmethod
    def get_template(cls, *args, **kwargs):
        '''returns a JSON formatted string that could be stored
        in a file, edited, loaded and used as parameters to run
        this workflow.'''
        return cls.Helper.get_template(*args, wf_class=cls, **kwargs)
