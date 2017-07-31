from helper import WorkflowHelper
from aiida.orm import Workflow


class ProjectionsWorkflow(Workflow):

    '''
    AiiDA-VASP Workflow for continuing from an NSCF Calculation to get
    all the data necessary to run a wannier calculation.
    parameters are given using :py:func:set_params(parameter_dict).
    see py:func:get_params_template() for a list of parameters.
    '''
    Helper = WorkflowHelper

    def __init__(self, **kwargs):
        self.helper = self.Helper(parent=self)
        super(ProjectionsWorkflow, self).__init__(**kwargs)

    def set_params(self, params):
        self.helper._verify_params(params)
        super(ProjectionsWorkflow, self).set_params(params)

    def get_calc_maker(self):
        from aiida.orm import Calculation
        params = self.get_parameters()
        cont = Calculation.query(uuid=params['continue_from'])[0]
        maker = self.helper._get_calc_maker(
            'vasp.amn', copy_from=cont)

        nscf_parameters = {'lwannier90': True,
                         'icharg': 11}
        maker.rewrite_parameters(**nscf_parameters)

        cout = cont.get_outputs_dict()
        cinp = cont.get_inputs_dict()
        maker.wannier_parameters = cout.get(
            'wannier_parameters',
            cinp.get(
                'wannier_parameters',
                {}))
        return maker

    @Workflow.step
    def start(self):
        from aiida_vasp.utils.win import modify_wannier_parameters_inline
        from aiida.orm import DataFactory
        params = self.get_parameters()
        self.append_to_report(self.helper._wf_start_msg())
        maker = self.get_calc_maker()

        old_win = maker.wannier_parameters or {}
        mod_win = params.get('wannier_parameters', {})
        mod_win['projections'] = params['projections']

        nbands = maker.parameters.get('nbands')
        if nbands:
            mod_win['num_bands'] = nbands
        ParameterData = DataFactory('parameter')
        mod_win = ParameterData(dict=mod_win)
        mod_c, mod_d = modify_wannier_parameters_inline(
            original=old_win, modifications=mod_win)
        msg = ('added projections and overrides to wannier parameters.')
        self.append_to_report(msg)
        self.add_attribute('modify_wannier_parameters_uuid', mod_c.uuid)

        maker.wannier_parameters = mod_d['wannier_parameters']

        calc = maker.new()
        calc.description = params.get('desc', maker.label)
        calc.store_all()
        calc.set_extras(params.get('extras'))

        self.attach_calculation(calc)
        self.append_to_report(
            self.helper._calc_start_msg('Amn Calculation', calc))
        self.next(self.end)

    @Workflow.step
    def end(self):
        calc = self.helper._get_first_step_calc(self.start)
        output_links = ['wannier_data']
        valid = self.helper._verify_calc_output(calc, output_links)
        wdat_valid, wdat_log = self._verify_wannier_data(calc.out.wannier_data)
        if valid and wdat_valid:
            self.add_result('calc', calc)
            self.add_result('wannier_parameters', calc.inp.wannier_parameters)
            self.add_result('wannier_data', calc.out.wannier_data)
            self.append_to_report(
                'Added the nscf calculation as a result')
        elif not wdat_valid:
            self.append_to_report(wdat_log)
        else:
            self.append_to_report(
                self.helper._calc_invalid_outs_msg(calc, output_links))
        self.next(self.exit)

    def _verify_wannier_data(self, wdat):
        valid = False
        log = ''
        names = set(wdat.archive.getnames())
        required = {'wannier90.mmn',
                    'wannier90.eig',
                    'wannier90.amn'}
        if required.issubset(names):
            valid = True
        else:
            log += ('the retrieved wannier_data node does not contain a .amn'
                    'file. '
                    'something must have gone wrong, no output produced.')
        return valid, log

    def _verify_param_kpoints(self, params):
        valid, log = self.helper._verify_kpoints(params)
        if params.get('use_wannier'):
            if not params['kpoints'].get('mesh'):
                log += ('{}: parameters: kpoints may only be given as a mesh '
                        'when using wannier.')
                valid = False
        return valid, log

    @classmethod
    def get_params_template(cls):
        '''returns a dictionary of keys and explanations how they
        can be used as parameters for this workflow.'''
        tmpl = cls.Helper.get_params_template(continuation=True)
        tmpl['projections'] = ['XX : s; px; py; pz', 'YY: ...']
        tmpl['wannier_parameters'] = {
            '#_explanation': (
                'overrides the parameters taken from continue_from '
                'keys starting with # are ignored'),
            '#num_wann': 'int',
            '#use_bloch_phases': 'false | true',
            }
        return tmpl

    @classmethod
    def get_template(cls, *args, **kwargs):
        '''returns a JSON formatted string that can be stored to a file,
        edited, loaded and used to run this Workflow.'''
        return cls.Helper.get_template(*args, wf_class=cls, **kwargs)
