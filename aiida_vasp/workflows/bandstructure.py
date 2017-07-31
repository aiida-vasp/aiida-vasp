# ~ import aiida.common
from aiida.common import aiidalogger
from aiida.orm.workflow import Workflow
from aiida.orm import DataFactory
from aiida_vasp.calcs.maker import VaspMaker

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
        maker = VaspMaker(structure=params['structure'])
        maker.add_parameters(icharg=0, istart=0, lorbit=11, lsorbit=True,
                           sigma=0.05, ismear=0, gga='PE', gga_compat=False)
        maker.code = Code.get_from_string(params['code'])
        # ~ maker.computer = maker.code.get_computer()
        maker.label = params.get('name')
        maker.queue = params['queue']
        return maker

    @Workflow.step
    def start(self):
        '''
        prepare, store and attach the selfconsistent run
        to get the charge density.
        '''
        params = self.get_parameters()
        maker = self.get_calc_maker()

        maker.kpoints.set_kpoints_mesh([8, 8, 8])
        maker.label += ': sc run'
        calc = maker.new()
        calc.description = '{}: selfconsistent run'.format(params.get('name'))
        calc.set_resources({'num_machines': 8, 'num_mpiprocs_per_machine': 2})
        calc.store_all()
        self.attach_calculation(calc)
        self.append_to_report(
            '{}: starting selfconsistent run, PK={}, uuid={}'.format(
                params.get('name'), calc.pk, calc.uuid))
        self.add_attributes({'scstep': {'pk': calc.pk, 'uuid': calc.uuid}})
        self.next(self.bandrun)

    @Workflow.step
    def bandrun(self):
        '''
        prepare, store and attach the non-selfconsistent run
        to get the band structure.
        '''
        params = self.get_parameters()
        # get chgcar from previous step
        scstep = self.get_attributes()['scstep']
        prev = self.get_step_calculations(self.start)
        sccalc = prev.get(uuid=scstep['uuid'])
        self.append_to_report(
            '{}: retrieved sc step calculation'.format(
                params.get('name')))

        # set up band structure step
        maker = self.get_calc_maker()
        maker._init_from(sccalc)
        maker.set_kpoints_path(value=params.get('kp_path'))
        maker.label += ': bands run'
        calc = maker.new()

        calc.description = '{}: bandstructure run'.format(params.get('name'))
        calc.set_resources({'num_machines': 8, 'num_mpiprocs_per_machine': 2})
        calc.store_all()

        self.attach_calculation(calc)
        self.append_to_report(
            '{}: starting band structure run, PK={}, uuid={}'.format(
                params.get('name'), calc.pk, calc.uuid))
        self.add_attributes({'bsstep': {'pk': calc.pk, 'uuid': calc.uuid}})
        self.add_result('band_run', calc)
        self.next(self.get_results)

    @Workflow.step
    def get_results(self):
        params = self.get_parameters()
        prev = self.get_step_calculations(self.bandrun)
        bandcalc = prev.get(uuid=self.get_attribute('bsstep')['uuid'])
        self.add_to_report(
            '{}: calculations done, gathering results'.format(
                params.get('name')))
        efermi = bandcalc.out.results.get_dict()['efermi']
        self.add_to_report(
            '{}: E_fermi = {}'.format(params.get('name'), efermi))
        self.add_result('efermi', bandcalc.out.results.get_dict()['efermi'])
        self.next(self.exit)
