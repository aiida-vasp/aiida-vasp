#~ import aiida.common
from aiida.common import aiidalogger
from aiida.orm.workflow import Workflow

logger = aiidalogger.getChild('Bandstructure')

class Bandstructure(Workflow):
    '''
    AiiDA workflow to get bandstructure of a material using VASP.

    Necessary steps are:
        * selfconsistent run with few kpoints
        * nonselfconsistent run with more kpoints and CHGCAR from
          selfconsistent run
    '''
    def __init__(self, **kwargs):
        self.queue = kwargs.get('queue')
        super(Bandstructure, self).__init__(**kwargs)

    def get_vasp_calc(self, sc=True):
        from aiida.orm import Code

        vaspcalc = Code.get_from_string('vasp.asevasp')
        vaspcalc.set_queue_name(self.queue)
        if sc:
            vaspcalc.set_max_wallclock_seconds(60)
            vaspcalc.set_withmpi(False)
            vaspcalc.set_resources({'num_machines':1, 'num_mpiprocs_per_machine':1})

    @Workflow.step
    def start(self):
        params = self.get_parameters()
        calc = get_vasp_calc(sc=True)
        incar = params['incar'].copy()
        incar.remove(
        calc.use_settings(params['incar'])
