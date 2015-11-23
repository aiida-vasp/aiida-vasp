#~ import aiida.common
from aiida.common import aiidalogger
from aiida.orm.workflow import Workflow
from aiida.orm import DataFactory, load_node

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

        #~ self.queue = kwargs.get('queue')
        #~ params.get('name') = kwargs.get('name') or 'unnamed'
        #~ self.infiles = kwargs['infiles']
        super(Bandstructure, self).__init__(**kwargs)

    def get_structure(self, sdict):
        SData = DataFactory('structure')
        st = SData(cell = sdict['cell'])
        for atom in sdict['atoms']:
            st.append_atom(position = atom['position'], symbols = atom['symbols'])
        return st

    def get_vasp_calc(self, sc=True):
        '''
        return an asevasp code instance.
        kwargs:
            sc=True: if True, sets scheduler parameters as for a selfconsistent run.
                        aka: max walltime = 60, no mpi, single core.
        '''
        from aiida.orm import Code

        params = self.get_parameters()
        vaspcode = Code.get_from_string('asevasp')
        vaspcalc = vaspcode.new_calc()
        vaspcalc.label = params.get('name')
        if params.get('queue'):
            vaspcalc.set_queue_name(params['queue'])
        if sc:
            vaspcalc.set_max_wallclock_seconds(600)
            vaspcalc.set_withmpi(False)
            vaspcalc.set_resources({'num_machines':1, 'num_mpiprocs_per_machine':1})
        return vaspcalc

    def get_formula(self):
        params = self.get_parameters()
        return load_node(params['structure']).get_formula()

    @Workflow.step
    def start(self):
        '''
        prepare, store and attach the selfconsistent run to get the charge density.
        '''
        params = self.get_parameters()
        calc = self.get_vasp_calc(sc=True)

        calc.use_settings(load_node(params['incar']))
        calc.use_structure(load_node(params['structure']))
        calc.use_potentials(load_node(params['potentials']))
        calc.use_kpoints(load_node(params['sc_kpoints']))
        calc.use_chgcar(ParameterData(dict={}))
        formula = self.get_formula()
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
        if prev:
            #~ sccalc = [i for i in prev if i.uuid == scstep['uuid']][0]
            sccalc = prev.get(uuid=scstep['uuid'])
            self.append_to_report('{}: retrieved sc step calculation'.format(params.get('name')))

        ichg = sccalc.out.chgcar.copy()

        # set up band structure step
        calc = self.get_vasp_calc(sc=False)

        incar = params['incar'].get_dict()
        incar['icharg'] = 11
        settings = ParameterData(dict=incar)
        calc.use_settings(settings)
        calc.use_structure(load_node(params['structure']))
        calc.use_potentials(load_node(params['potentials']))
        calc.use_kpoints(load_node(params['sc_kpoints']))
        #~ structure = ParameterData(dict=params['structure'])
        #~ calc.use_structure(self.get_structure(params['structure']))
        #~ potentials = ParameterData(dict=params['potentials'])
        #~ calc.use_potentials(potentials)
        #~ kpoints = ParameterData(dict=params['bs_kpoints'])
        #~ calc.use_kpoints(kpoints)
        calc.use_chgcar(ichg)

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
        self.set_result('band_run', calc)
        self.next(self.exit)
