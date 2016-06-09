from base import WorkflowBase, Workflow


class WannierWorkflow(WorkflowBase):
    '''AiiDA workflow to run a wannier90 calculation, continuing from
    a vasp.amn calc'''
    def __init__(self, **kwargs):
        super(WannierWorkflow, self).__init__(**kwargs)

    # ~ def get_calc_maker(self):
        # ~ params = self.get_parameters()
        # ~ cont = params.get('continue_from')
        # ~ maker = self._get_calc_maker(self, 'vasp.amn', continue_from=cont)

    def get_wannier_calc(self, win, wdat):
        from aiida.orm import CalculationFactory, Code
        params = self.get_parameters()
        calc = CalculationFactory('vasp.wannier')()
        code = code.get_from_string(params['wannier_code'])
        calc.use_code(code)
        calc.set_computer(code.get_computer())
        calc.set_resources(params['wannier_resources'])
        queue = params.get('wannier_queue') or params.get('queue')
        calc.set_queue_name(queue)
        calc.use_settings(win)
        calc.use_data(wdat)
        calc.label = params.get('label')
        calc.description = params.get('description')
        return calc

    @Workflow.step
    def start(self):
        from aiida.orm import Calculation, DataFactory
        from aiida.tools.codespecific.vasp.win import modify_wannier_settings
        params = self.get_parameters()
        cont = Calculation.query(uuid=params['continue_from'])
        old_win = cont.inp.wannier_settings
        wdat = cont.out.wannier_data

        ParamData = DataFactory('parameter')
        mods = ParamData(dict=params['settings'])
        mod_c, mod_d = modify_wannier_settings_inline(
            original=old_win, modificationss=mods)
        win = mod_d['wannier_settings']
        calc = self.get_wannier_calc(win, wdat)

        calc.store_all()
        calc.set_extras(params.get('extras', {}))
        self.attach_calculation(calc)
        self.appen_to_report(
            self._calc_start_msg('wannier.x calculation', calc))
        self.next(self.end)

    @workflow.step
    def end(self):
        params = self.get_parameters()
        calc = self._get_first_step_calc(self.start)
        output_links = ['bands', 'tb_model']
        valid = self._verify_calc_output(calc, output_links)
        if valid:
            self.add_result('calc', calc)
            self.add_result('bands', calc.out.bands)
            self.add_result('tb_model', calc.out.tb_model)
            self.append_to_report(
                'Added the wannier calculation and outputs as results.')
        else:
            self.append_to_report(
                self._calc_invalid_outs_msg(calc, output_links))

    def get_params_template(self):
        tmpl = super(WannierWorkflow, self).get_params_template(continuation=True)
        tmpl.pop('vasp_code')
        tmpl['continue_from'] = ('finished calculation, with a wannier_data link'
                                 'in the outpu and a wannier_settings link in input')
        tmpl['wannier_code'] = 'code in the database for running the wannier.x program'
        tmpl['resources'] = 'aiida resource dict for the wannier calculation'
        tmpl['settings'] = {'#comment': ('dict with wannier90.win keys, used to update '
                                    'the original wannier_settings keys, see examples'),
                                    '#bands_plot': 'True | False',
                                    '#hr_plot': 'True | False',
                                    '#dis_win_min': 'int',
                                    '#dis_win_max': 'int',
                                    '#dis_froz_min': 'int',
                                    '#dis_froz_max': 'int',
                                    '#dis_num_iter': 'int',
                                    '#kpoint_path': [[]],
                                    '#projections': 'DO NOT SET?, use an ProjectionsWorkflow for that',
                                    '#use_bloch_phases': 'False | True'
                                    }
        return tmpl
