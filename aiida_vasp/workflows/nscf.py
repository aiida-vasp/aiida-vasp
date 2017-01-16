from helper import WorkflowHelper
from aiida.orm import Workflow, Calculation


class NscfWorkflow(Workflow):

    '''
    AiiDA-VASP Workflow for continuing from an SCF Calculation.
    Can be used with or without the vasp2wannier90 interface.
    '''
    Helper = WorkflowHelper

    def __init__(self, **kwargs):
        self.helper = self.Helper(parent=self)
        super(NscfWorkflow, self).__init__(**kwargs)

    def get_calc_maker(self):
        params = self.get_parameters()
        maker = self.helper._get_calc_maker(
            'vasp.nscf', continue_from=Calculation.query(
                uuid=params['continue_from'])[0])
        nscf_settings = {'lwannier90': params['use_wannier'],
                         'icharg': 11}
        maker.rewrite_settings(**nscf_settings)
        return maker

    @Workflow.step
    def start(self):
        params = self.get_parameters()
        self.append_to_report(self.helper._wf_start_msg())
        maker = self.get_calc_maker()
        calc = maker.new()

        calc.description = params.get('desc', maker.label)
        calc.store_all()
        calc.set_extras(params.get('extras'))

        self.attach_calculation(calc)
        self.append_to_report(
            self.helper._calc_start_msg('NSCF Calculation', calc))
        self.next(self.end)

    @Workflow.step
    def end(self):
        params = self.get_parameters()
        calc = self.helper._get_first_step_calc(self.start)
        output_links = ['bands', 'dos']
        if params.get('use_wannier'):
            output_links += ['wannier_settings']
        valid = self.helper._verify_calc_output(calc, output_links)
        if valid:
            self.add_result('calc', calc)
            self.append_to_report(
                'Added the nscf calculation as a result')
        else:
            self.append_to_report(
                self.helper._calc_invalid_outs_msg(calc, output_links))
        self.next(self.exit)

    def _verify_param_kpoints(self, params):
        valid, log = self.helper._verify_kpoints(params)
        if params.get('use_wannier') and params.get('kpoints'):
            if not params['kpoints'].get('mesh'):
                log += ('{}: parameters: kpoints may only be given as a mesh '
                        'when using wannier.').format(self.__class__.__name__)
                valid = False
        return valid, log

    def set_params(self, params):
        self.helper._verify_params(params)
        super(NscfWorkflow, self).set_params(params)

    @classmethod
    def get_params_template(cls):
        '''returns a dictionary of keys and explanations how they
        can be used as parameters for this workflow.'''
        tmpl = cls.Helper.get_params_template(continuation=True)
        tmpl['use_wannier'] = ('True | False (if true, vasp_code must be '
                               'compiled with wannier interface')
        return tmpl

    @classmethod
    def get_template(cls, *args, **kwargs):
        '''returns a JSON formatted string that can be stored to a file,
        edited, loaded and used to run this Workflow.'''
        return cls.Helper.get_template(*args, wf_class=cls, **kwargs)
