from base import WorkflowBase, Workflow


class ScfWorkflow(WorkflowBase):
    '''AiiDA workflow to run a VASP scf calculation
    to reuse WAVECAR and CHGCAR'''
    def __init__(self, **kwargs):
        super(ScfWorkflow, self).__init__(**kwargs)

    def get_calc_maker(self):
        maker = self._get_calc_maker('vasp.scf')
        maker.add_settings(icharg=0, istart=0)
        return maker

    def get_params_template(self):
        return super(ScfWorkflow, self).get_params_template(continuation=False)

    @Workflow.step
    def start(self):
        params = self.get_parameters()
        maker = self.get_calc_maker()
        calc = maker.new()
        calc.description = params.get('description', '')
        calc.store_all()
        calc.set_extras(params.get('extras'))
        self.attach_calculation(calc)
        self.append_to_report(self._calc_start_msg('scf VASP run', calc))
        self.next(self.end)

    @Workflow.step
    def end(self):
        calc = self._get_first_step_calc(self.start)
        output_links = ['charge_density', 'wavefunctions']
        valid_out = self._verify_calc_output(calc, output_links)
        if valid_out:
            self.add_result('calc', calc)
            self.append_to_report(
                'Added the scf calculation as a result')
        else:
            self.append_to_report(
                self._calc_invalid_outs_msg(calc, output_links))
