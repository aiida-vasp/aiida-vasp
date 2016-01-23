#~ import aiida.common
from aiida.common import aiidalogger
from aiida.orm.workflow import Workflow
from aiida.orm import DataFactory, load_node
from aiida.orm.calculation.job.vasp.maker import VaspMaker

logger = aiidalogger.getChild('Bandstructure')
ParameterData = DataFactory('parameter')

class Bandstructure(Workflow):
    '''
    AiiDA workflow to get bandstructure of a material using VASP.

    Necessary steps are:
        * selfconsistent run with few kpoints
        * nonselfconsistent run with more kpoints and CHGCAR from
          selfconsistent run
    '''
    def __init__(self, **kwargs):
        super(Bandstructure, self).__init__(**kwargs)

    def get_calc_maker(self):
        '''
        return an VaspMaker instance.
        '''
        from aiida.orm import Code
        params = self.get_parameters()
        maker = VaspMaker(params['structure'])
        maker.add_settings(icharg=0, istart-0, lorbit=11, lsorbit=True,
                sigma=0.05, ismear=0, gga='PE', gga_compat=False)
        maker.code = Code.get_from_string(params['code'])
        maker.computer = maker.code.get_computer()
        maker.label = params.get('name')
        maker.queue = params['queue']
        return maker

    @Workflow.step
    def start(self):
        '''
        prepare, store and attach the selfconsistent run to get the charge density.
        '''
        params = self.get_parameters()
        maker = self.get_calc_maker()

        maker.kpoints.set_kpoints_mesh([8, 8, 8])
        calc = calc.new()
        calc.description = '{}: selfconsistent run for {}'.format(params.get('name'), formula)
        calc.store_all()
        self.attach_calculation(calc)
        self.append_to_report('{}: starting selfconsistent run, PK={}, uuid={}'.format(
            params.get('name'), calc.pk, calc.uuid))
        self.add_attributes({'scstep': {'pk': calc.pk, 'uuid': calc.uuid}})
        self.next(self.bandrun)

    @Workflow.step
    def bandrun(self):
        '''
        prepare, store and attach the non-selfconsistent run to get the band structure.
        '''
        params = self.get_parameters()
        # get chgcar from previous step
        scstep = self.get_attributes()['scstep']
        prev = self.get_step_calculations(self.start)
        sccalc = prev.get(uuid=scstep['uuid'])
        self.append_to_report('{}: retrieved sc step calculation'.format(params.get('name')))

        ichgdens = sccalc.out.chgcar

        # set up band structure step
        maker = self.get_calc_maker()
        maker.rewrite_settings(icharg=11)
        maker.kpoints.set_kpoints_path()
        
        maker._charge_density = ichgdens

        formula = self.get_formula()
        calc.description = '{}: bandstructure run for {}'.format(params.get('name'), formula)
        calc.set_max_wallclock_seconds(1800)
        calc.set_withmpi(True)
        nbands = incar.get('nbands') or incar.get('NBANDS')
        resources = {}
        resources['num_machines'] = nbands / 10 + (nbands%10 and 1 or 0)
        resources['tot_num_mpiprocs'] = nbands
        calc.set_resources(resources)
        calc.store_all()

        self.attach_calculation(calc)
        self.append_to_report('{}: starting band structure run, PK={}, uuid={}'.format(
            params.get('name'), calc.pk, calc.uuid))
        self.add_attributes({'bsstep': {'pk': calc.pk, 'uuid': calc.uuid}})
        self.add_result('band_run', calc)
        self.next(self.exit)
