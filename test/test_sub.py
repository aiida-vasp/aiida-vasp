from aiida_vasp.calcs.vasp import VaspCalculation
from aiida.orm import DataFactory, Code
from aiida.common.folders import Folder
import pymatgen as pmg
import sys

#~ vc = VaspCalculation()

incar = {'SYSTEM': 'TestSystem',
         'ediff': 1E-5,
         'GGA_COMPAT': False,
         }
poscar = pmg.io.vasp.Poscar.from_file(sys.argv[1]).as_dict()
potcar = {}
kpoints = {}

ParameterData = DataFactory('parameter')
code = Code.get_from_string('vasp')

vc = code.new_calc()

vc.use_settings(ParameterData(dict=incar))
vc.use_structure(ParameterData(dict=poscar))
vc.use_potentials(ParameterData(dict=potcar))
vc.use_kpoints(ParameterData(dict=kpoints))

vc.label = 'Test Vasp Calculation'
vc.description = 'Test Vasp Calculation Plugin (with invalid data)'
vc.set_max_wallclock_seconds(1) #should not run at all
vc.set_withmpi(True)
vc.set_resources({'num_machines':2, 'num_mpiprocs_per_machine':12})
vc.set_queue_name('bogus_queue_to_prevent_accidental_queuing')

fldr = Folder(sys.argv[2])
vc.submit_test(fldr)
