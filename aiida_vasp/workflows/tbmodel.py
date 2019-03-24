"""AiiDA - VASP workflow to get tight binding model of a material using VASP and wannier90"""
from aiida.common import aiidalogger
from aiida.common.exceptions import NotExistent
from aiida.plugins import DataFactory, Group, CalculationFactory, Code
from aiida.orm.workflow import Workflow

from aiida_vasp.calcs.maker import VaspMaker

LOGGER = aiidalogger.getChild('Tbmodel')
PARAMETER_CLS = DataFactory('parameter')


class TbmodelWorkflow(Workflow):
    """
    AiiDA workflow to get tight binding model of a material using VASP and wannier90.

    Necessary steps are:
        * selfconsistent run with few kpoints
        * nonselfconsistent run with lwannier90
        * lwannier90 run with projections block
        * wannier.x run with hr_plot = True
    """

    def __init__(self, **kwargs):
        super(TbmodelWorkflow, self).__init__(**kwargs)
        self.group = None
        try:
            self.group = Group.get_from_string('tbmodel')
        except NotExistent:
            self.group = Group(name='tbmodel')
            self.group.store()

    def get_calc_maker(self):
        """Return an VaspMaker instance."""
        params = self.get_parameters()
        maker = VaspMaker(structure=params['structure'])
        maker.add_parameters(icharg=0, istart=0, lorbit=11, lsorbit=True, sigma=0.05, ismear=0, gga='PE', gga_compat=False)
        maker.code = Code.get_from_string(params['vasp'])
        maker.label = params.get('name')
        return maker

    @Workflow.step
    def start(self):
        """Prepare, store and attach the selfconsistent run to get the charge density."""
        params = self.get_parameters()
        maker = VaspMaker(structure=params['structure'], calc_cls='vasp.scf')
        maker.rewrite_parameters(**params['parameters'])
        kpoints = params['kpmesh']
        maker.kpoints.set_kpoints_mesh(kpoints)
        maker.code = Code.get_from_string(params['vasp'])
        maker.queue = params['queue']
        maker.resources['num_machines'] = 4
        maker.resources['num_mpiprocs_per_machine'] = 2
        maker.label = params.get('name') + ': sc run'
        calc = maker.new()
        calc.description = '{}: selfconsistent run'.format(params.get('name'))
        calc.store_all()
        calc.set_extras({'experiment': 'tight-binding'})
        self.group.add_nodes(calc)
        self.attach_calculation(calc)
        self.append_to_report('{}: starting selfconsistent run, PK={}, uuid={}'.format(params.get('name'), calc.pk, calc.uuid))
        self.add_attributes({'scstep': {'pk': calc.pk, 'uuid': calc.uuid}})
        self.next(self.winrun)

    @Workflow.step
    def winrun(self):
        """Prepare, store and attach the non-selfconsistent run to get the wannier input file"""
        params = self.get_parameters()
        # get chgcar from previous step
        scstep = self.get_attributes()['scstep']
        prev = self.get_step_calculations(self.start)
        sccalc = prev.get(uuid=scstep['uuid'])
        self.append_to_report('{}: retrieved sc step calculation'.format(params.get('name')))

        # set up band structure step
        maker = VaspMaker(calc_cls='vasp.nscf', continue_from=sccalc)
        maker.add_parameters(lwannier90=True, icharg=11)
        maker.label = params.get('name') + ': win run'
        calc = maker.new()

        calc.description = '{}: win run'.format(params.get('name'))
        calc.store_all()
        calc.set_extras({'experiment': 'tight-binding'})

        self.group.add_nodes(calc)
        self.attach_calculation(calc)
        self.append_to_report('{}: starting win run, PK={}, uuid={}'.format(params.get('name'), calc.pk, calc.uuid))
        self.add_attributes({'winstep': {'pk': calc.pk, 'uuid': calc.uuid}})
        self.next(self.amnrun)

    @Workflow.step
    def amnrun(self):
        """Launch the Vasp2Wannier calculation to get the .amn file"""
        params = self.get_parameters()
        winstep = self.get_attributes()['winstep']
        prev = self.get_step_calculations(self.winrun)
        wincalc = prev.get(uuid=winstep['uuid'])
        self.append_to_report('{}: retrieved nscf calculation'.format(params.get('name')))

        maker = VaspMaker(calc_cls='vasp.amn', copy_from=wincalc)
        maker.wannier_parameters = wincalc.out.wannier_parameters.copy()
        num_bands = maker.wannier_parameters.get_dict()['num_wann']
        maker.wannier_parameters.update_dict({'num_bands': num_bands})
        maker.wannier_parameters.update_dict(params['win'])
        maker.label = params.get('name') + ': amn run'
        calc = maker.new()

        calc.description = '{}: amn run'.format(params.get('name'))
        calc.store_all()
        calc.set_extras({'experiment': 'tight-binding'})

        self.group.add_nodes(calc)
        self.attach_calculation(calc)
        self.append_to_report('{}: starting amn run, PK={}, uuid={}'.format(params.get('name'), calc.pk, calc.uuid))
        self.add_attributes({'amnstep': {'pk': calc.pk, 'uuid': calc.uuid}})
        self.next(self.wannrun)

    @Workflow.step
    def wannrun(self):
        """Launch the Wannier90 calculation that returns the final tight-binding model and set results"""
        params = self.get_parameters()
        amnstep = self.get_attributes()['amnstep']
        prev = self.get_step_calculations(self.amnrun)
        amncalc = prev.get(uuid=amnstep['uuid'])
        self.append_to_report('{}: retrieved amn calculation'.format(params.get('name')))

        calc = CalculationFactory('vasp.wannier')()
        code = Code.get_from_string(params['wannier_x'])
        calc.use_code(code)
        calc.set_computer(code.get_computer())
        calc.use_parameters(amncalc.inp.wannier_parameters)
        calc.use_data(amncalc.out.wannier_data)
        calc.label = params.get('name') + ': wannier run'
        calc.set_resources({'num_machines': 1})
        calc.set_queue_name(amncalc.get_queue_name())
        calc.description = calc.label

        calc.store_all()
        calc.set_extras({'experiment': 'tight-binding'})

        self.group.add_nodes(calc)
        self.attach_calculation(calc)
        self.add_result('wannier_run', calc)
        self.append_to_report('{}: starting wannier run, PK={}, uuid={}'.format(params.get('name'), calc.pk, calc.uuid))
        self.add_attributes({'wannstep': {'pk': calc.pk, 'uuid': calc.uuid}})
        self.next(self.exit)
