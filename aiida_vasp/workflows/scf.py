"""AiiDA - VASP workflow to run a VASP scf calculation"""
from aiida.common import aiidalogger
from aiida.orm import Workflow

from .helper import WorkflowHelper

LOGGER = aiidalogger.getChild('ScfWorkflow')


class ScfWorkflow(Workflow):
    """
    AiiDA workflow to run a VASP scf calculation

    the resulting WAVECAR and CHGCAR can then be reused in further workflows.
    """
    Helper = WorkflowHelper

    def __init__(self, **kwargs):
        self.helper = self.Helper(parent=self)
        super(ScfWorkflow, self).__init__(**kwargs)

    def get_calc_maker(self):
        maker = self.helper._get_calc_maker('vasp.scf')  # pylint: disable=protected-access
        maker.add_parameters(icharg=0, istart=0)
        return maker

    @Workflow.step
    # pylint: disable=protected-access
    def start(self):
        """Submit the calculation"""
        params = self.get_parameters()
        self.append_to_report(self.helper._wf_start_msg())
        maker = self.get_calc_maker()
        calc = maker.new()
        calc.description = params.get('description', '')
        calc.store_all()
        calc.set_extras(params.get('extras'))
        self.attach_calculation(calc)
        self.append_to_report(self.helper._calc_start_msg('scf VASP run', calc))
        self.next(self.end)

    @Workflow.step
    # pylint: disable=protected-access
    def end(self):
        """Set results"""
        calc = self.helper._get_first_step_calc(self.start)
        output_links = ['charge_density', 'wavefunctions']
        valid_out = self.helper._verify_calc_output(calc, output_links)
        if valid_out:
            self.add_result('calc', calc)
            self.append_to_report('Added the scf calculation as a result')
        else:
            self.append_to_report(self.helper._calc_invalid_outs_msg(calc, output_links))
        self.next(self.exit)

    def set_params(self, params, force=False):
        self.helper._verify_params(params)  # pylint: disable=protected-access
        super(ScfWorkflow, self).set_params(params, force=force)

    @classmethod
    def get_params_template(cls):
        """Gives a dictionary: keys are workflow parameters, values are explanations."""
        return cls.Helper.get_params_template(continuation=False)

    @classmethod
    def get_template(cls, *args, **kwargs):
        """
        Returns a JSON formatted template string.

        The template can be stored in a file, edited, loaded and used as parameters to run this workflow.
        """
        return cls.Helper.get_template(*args, wf_class=cls, **kwargs)
