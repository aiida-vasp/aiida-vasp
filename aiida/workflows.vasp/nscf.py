from helper import WorkflowHelper
from aiida.orm import Workflow, Calculation


class NscfWorkflow(Workflow):
    '''
    AiiDA-VASP Workflow for continuing from an SCF Calculation
    parameters are given using :py:func:set_params(parameter_dict).
    See below for a list of keys to use in parameter_dict.
    :key str vasp_code: the identifier for your vasp code,
        example: "code@computer"
    :key str queue: the name of the queue to run calculations in.
    :key dict resources: the aiida style resources specification. example:
        {'num_machines': 4, 'num_mpiprocs_per_machine': 16}
    :key dict kpoints: a dict containing one of the keys ['mesh', 'list', 'path']
        The value of that key must be suitable to create a KpointsData node.
    :key bool use_wannier: default=False, if True LWANNIER90 is used and the
        wannier input/output are retrieved as well.
        This requires the vasp code to be compiled with the VASP2WANNIER90
        flag.
    :key dict settings: optional, a dict with incar keys, overriding the ones
        of the original scf calc. Note: by default ICHARG is set to 11,
        LWANNIER90 depending on the 'use_wannier' key.
    :key str label: optional, overrides the default label given to the
        calculation.
    :key desc: optional, overrides the default description given to the
        calculation.
    :key dict continue_from: uuid of the scf calculation to continue from
    '''
    def __init__(self, **kwargs):
        self.helper = WorkflowHelper(parent=self)
        super(NscfWorkflow, self).__init__(**kwargs)

    def get_calc_maker(self):
        params = self.get_parameters()
        maker = self.helper._get_calc_maker(
            'vasp.nscf', continue_from=Calculation.query(uuid=params['continue_from'])[0])
        nscf_settings = {'lwannier90': params['use_wannier'],
                         'icharg': 11}
        maker.rewrite_settings(**nscf_settings)
        return maker

    @Workflow.step
    def start(self):
        params = self.get_parameters()
        maker = self.get_calc_maker()
        calc = maker.new()

        calc.description = params.get('desc', maker.label)
        calc.store_all()

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

    def get_params_template(self):
        tmpl = self.helper.get_params_template(continuation=True)
        tmpl['use_wannier'] = ('True | False (if true, vasp_code must be '
                              'compiled with wannier interface')
        return tmpl

    def get_template(self, *args, **kwargs):
        return self.helper.get_template(*args, **kwargs)
