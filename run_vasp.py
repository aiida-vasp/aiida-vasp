from aiida.work.run import run
from aiida.orm import CalculationFactory, Code
from aiida.orm.data.base import Str, Int
from aiida.orm.data.array.kpoints import KpointsData
from aiida.orm.data.parameter import ParameterData
from aiida.workflows.user.base import VASPBaseWorkChain

VaspCalculation = CalculationFactory('vasp.vasp')

options = {
	'resources': {
		'num_machines': 1,
		'tot_num_mpiprocs': 1,
	},
	'max_wallclock_seconds': 1800,
}

kpoints = KpointsData()
kpoints.set_kpoints_mesh([1, 1, 1])

inputs = {
	'code': Code.get_from_string('VASP.5.4.4@Raichu'),
	'structure': load_node(888),
	'kpoints': kpoints,
	'parameters': ParameterData(dict={}),
	'settings': ParameterData(dict={}),
	'pseudo_family': Str('vasp-pbe'),
        'options' : ParameterData( dict = { 
                      'max_wallclock_seconds' : 3600,
                      'max_memory_kb': 10000000,
                      'resources' : { 'num_machines': 1
                                    },
                    }),
        'max_iterations' : Int(1),
}

process = VaspCalculation.process()
# running = run(process, **inputs)
running = run(VASPBaseWorkchain, **inputs)
