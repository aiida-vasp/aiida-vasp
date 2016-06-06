from base import WorkflowBase, Workflow


class WannierWorkflow(WorkflowBase):
    '''AiiDA workflow to run a wannier90 calculation'''
    def __init__(self, **kwargs):
        super(WannierWorkflow, self).__init__(**kwargs)

    def get_calc_maker(self):
        params = self.get_parameters()
        cont = params.get('continue_from')
        maker = self._get_calc_maker(self, 'vasp.amn', continue_from=cont)

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

        calc.store_all()
        calc.set_extras(params.get('extras', {}))

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
                                    '#projections': 'DO NOT SET, use an ProjectionsWorkflow for that',
                                    '#use_bloch_phases': 'False | True'
                                    }
        return tmpl
