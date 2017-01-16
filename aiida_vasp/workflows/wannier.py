from aiida.orm import Workflow
from helper import WorkflowHelper


class WannierWorkflow(Workflow):

    '''AiiDA workflow to run a wannier90 calculation, continuing from
    a vasp.amn calc'''
    Helper = WorkflowHelper

    def __init__(self, **kwargs):
        self.helper = self.Helper(parent=self)
        super(WannierWorkflow, self).__init__(**kwargs)

    def get_wannier_calc(self, win, wdat):
        from aiida.orm import CalculationFactory, Code
        params = self.get_parameters()
        calc = CalculationFactory('vasp.wannier')()
        code = Code.get_from_string(params['wannier_code'])
        calc.use_code(code)
        calc.set_computer(code.get_computer())
        calc.set_resources(params['resources'])
        queue = params.get('queue') or params.get('queue')
        calc.set_queue_name(queue)
        calc.use_settings(win)
        calc.use_data(wdat)
        calc.label = params.get('label')
        calc.description = params.get('description')
        return calc

    @Workflow.step
    def start(self):
        from aiida.orm import Calculation, DataFactory, load_node
        from aiida_vasp.utils.win import modify_wannier_settings_inline
        params = self.get_parameters()
        self.append_to_report(self.helper._wf_start_msg())
        cont = params['continue_from']
        if isinstance(cont, dict):
            old_win = load_node(uuid=cont['settings'])
            wdat = load_node(uuid=cont['data'])
        else:
            cont = Calculation.query(uuid=cont)[0]
            old_win = cont.inp.wannier_settings
            wdat = cont.out.wannier_data

        ParamData = DataFactory('parameter')
        mods = ParamData(dict=params['settings'])
        mod_c, mod_d = modify_wannier_settings_inline(
            original=old_win, modifications=mods)
        win = mod_d['wannier_settings']
        calc = self.get_wannier_calc(win, wdat)

        calc.store_all()
        calc.set_extras(params.get('extras', {}))
        self.attach_calculation(calc)
        self.append_to_report(
            self.helper._calc_start_msg('wannier.x calculation', calc))
        self.next(self.end)

    @Workflow.step
    def end(self):
        params = self.get_parameters()
        wset = params['settings']
        calc = self.helper._get_first_step_calc(self.start)
        output_links = ['bands', 'tb_model']
        valid = self.helper._verify_calc_output(calc, output_links)
        if valid:
            self.add_result('calc', calc)
            if calc.get_outputs_dict().get('bands'):
                self.add_result('bands', calc.out.bands)
            if calc.get_outputs_dict().get('tb_model'):
                self.add_result('tb_model', calc.out.tb_model)
            self.append_to_report(
                'Added the wannier calculation and outputs as results.')
        else:
            self.append_to_report(
                self.helper._calc_invalid_outs_msg(calc, output_links))
        self.next(self.exit)

    @classmethod
    def get_params_template(cls):
        '''returns a dictionary of keys and explanations how they
        can be used as parameters for this workflow.'''
        tmpl = cls.Helper.get_params_template(continuation=True)
        tmpl.pop('vasp_code')
        tmpl.pop('kpoints')
        tmpl['continue_from'] = ('finished calculation, with a wannier_data link'
                                 'in the output and a wannier_settings link in input\n'
                                 'or a dict with keys [settings, data], and uuids for values')
        tmpl['wannier_code'] = 'code in the database for running the wannier.x program'
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

    @classmethod
    def get_template(cls, *args, **kwargs):
        '''returns a JSON formatted string that can be stored to a file,
        edited, loaded and used to run this Workflow.'''
        return cls.Helper.get_template(*args, wf_class=cls, **kwargs)
