from base import WorkflowBase, Workflow


class ProjectionsWorkflow(WorkflowBase):
    '''
    AiiDA-VASP Workflow for continuing from an NSCF Calculation to get
    all the data necessary to run a wannier calculation.
    parameters are given using :py:func:set_params(parameter_dict).
    see py:func:get_params_template() for a list of parameters.
    '''
    def __init__(self, **kwargs):
        super(ProjectionsWorkflow, self).__init__(**kwargs)

    def get_calc_maker(self):
        from aiida.orm import Calculation
        params = self.get_parameters()
        cont = Calculation.query(uuid=params['continue_from'])[0]
        maker = self._get_calc_maker(
            'vasp.amn', copy_from=cont)

        nscf_settings = {'lwannier90': True,
                         'icharg': 11}
        maker.rewrite_settings(**nscf_settings)

        cout = cont.get_outputs_dict()
        cinp = cont.get_inputs_dict()
        maker.wannier_settings = cout.get('wannier_settings', cinp.get('wannier_settings', {}))
        return maker

    @Workflow.step
    def start(self):
        from aiida.tools.codespecific.vasp.win import modify_wannier_settings_inline
        params = self.get_parameters()
        maker = self.get_calc_maker()

        old_win = maker.wannier_settings or {}
        mod_win = params.get('wannier_settings', {})
        mod_win['projections'] = params['projections']
        mod_c, mod_d = modify_wannier_settings_inline(original=old_win, modifications=mod_win)

        maker.wannier_settings = mod_d['wannier_settings']

        calc = maker.new()
        calc.description = params.get('desc', maker.label)
        calc.store_all()

        self.attach_calculation(calc)
        self.append_to_report(
            self._calc_start_msg('Amn Calculation', calc))
        self.next(self.end)

    @Workflow.step
    def end(self):
        params = self.get_parameters()
        calc = self._get_first_step_calc(self.start)
        output_links = ['wannier_settings', 'wannier_data']
        params['use_wannier'] and output_links += ['wannier_settings']
        valid = self._verify_calc_output(calc, output_links)
        wdat_valid, wdat_log = self._verify_wannier_data(calc.out.wannier_data)
        if valid and wdat_valid:
            self.add_result('calc', calc)
            self.append_to_report(
                'Added the nscf calculation as a result')
        elif not wdat_valid:
            self.append_to_report(wdat_log)
        else:
            self.append_to_report(
                self._calc_invalid_outs_msg(calc, output_links))

    def _verify_wannier_data(self, wdat):
        valid = False
        log = ''
        names = wdat.archive.getnames()
        required = ['wannier90.mmn',
                    'wannier90.eig',
                    'wannier90.amn']
        if required in names:
            valid = True
        else:
            log += ('the retrieved wannier_data node does not contain a .amn file. '
                    'something must have gone wrong, no output produced.')
        return valid, log

    def _verify_kpoints(self, params):
        valid, log = super(ProjectionsWorkflow, self)._verify_kpoints(params)
        if params.get('use_wannier'):
            if not params['kpoints'].get('mesh'):
                log += ('{}: parameters: kpoints may only be given as a mesh '
                        'when using wannier.')
                valid = False
        return valid, log

    def get_params_template(self):
        tmpl = super(ProjectionsWorkflow, self).get_params_template(continuation=True)
        tmpl['projections'] = ['XX : s; px; py; pz', 'YY: ...']
        tmpl['wannier_settings'] = {'#_explanation': 'overrides the settings taken from continue_from',
                                    '#num_wann': 'int',
                                    '#use_bloch_phases': 'False | True',
                                    }
