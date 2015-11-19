#~ import aiida.common
from aiida.common import aiidalogger
from aiida.orm.workflow import Workflow
from aiida.orm import DataFactory

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

        self.queue = kwargs.get('queue')
        self.name = kwargs.get('name') or 'unnamed'
        super(Bandstructure, self).__init__(**kwargs)

        def get_vasp_calc(self, sc=True):
            '''
            return an asevasp code instance.
            kwargs:
                sc=True: if True, sets scheduler parameters as for a selfconsistent run.
                         aka: max walltime = 60, no mpi, single core.
            '''
            from aiida.orm import Code

            vaspcalc = Code.get_from_string('vasp.asevasp')
            vaspcalc.label = self.name
            if self.queue:
                vaspcalc.set_queue_name(self.queue)
            if sc:
                vaspcalc.set_max_wallclock_seconds(60)
                vaspcalc.set_withmpi(False)
                vaspcalc.set_resources({'num_machines':1, 'num_mpiprocs_per_machine':1})

        def get_formula(self):
            params = self.get_parameters()
            return params['structure'].get_formula()

        @Workflow.step
        def start(self):
            '''
            prepare, store and attach the selfconsistent run to get the charge density.
            '''
            params = self.get_parameters()
            calc = get_vasp_calc(sc=True)
            
            calc.use_settings(params['incar'])
            calc.use_structure(params['structure'])
            calc.use_potentials(params['potentials'])
            calc.use_kpoints(params['sc_kpoints'])
            calc.use_chgcar(ParameterData(dict={}))
            formula = self.get_formula()
            calc.description = '{}: selfconsistent run for {}'.format(self.name, formula)
            calc.store_all()
            self.attach_calculation(calc)
            self.append_to_report('{}: starting selfconsistent run, PK={}, uuid={}'.format(
                self.name, calc.pk, calc.uuid))
            self.add_attributes({'scstep': {'pk': calc.pk, 'uuid': calc.uuid}})
            self.next(self.bandrun)

        def bandrun(self):
            '''
            prepare, store and attach the non-selfconsistent run to get the band structure.
            '''
            # get chgcar from previous step
            scstep = self.get_attributes()['scstep']
            prev = self.get_step_calculations(self.start)
            if prev:
                sccalc = [i for i in prev if prev.uuid == scstep.uuid][0]
                self.append_to_report('{}: retrieved sc step calculation'.format(self.name))

            ichg = sccalc.out.chgcar.copy()

            # set up band structure step
            params = self.get_parameters()
            calc = get_vasp_calc(sc=False)

            incar = params['incar'].get_dict()
            incar['icharg'] = 11
            settings = ParameterData(dict=incar)
            calc.use_settings(settings)
            calc.use_structure(params['structure'])
            calc.use_potentials(params['potentials'])
            calc.use_kpoints(params['bs_kpoints'])
            calc.use_chgcar(ichg)

            formula = self.get_formula()
            calc.description = '{}: bandstructure run for {}'.format(self.name, formula)
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
                self.name, calc.pk, calc.uuid))
            self.add_attributes({'bsstep': {'pk': calc.pk, 'uuid': calc.uuid}})
            self.set_result('band_run', calc)
            self.next(self.exit)
