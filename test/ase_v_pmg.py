__doc__ = '''sets up two selfconsistent vasp calculations, and runs them on monch to compare ase vs pmg structure'''
from aiida.orm import DataFactory, Code
#~ from pymatgen import Lattice
from numpy import array, linalg
import ase
import sys
import os

ParameterData = DataFactory('parameter')
StructureData = DataFactory('structure')
#~ SinglefileData = DataFactory('singlefile')

# ---------- set up the incar info ----------
incar = ParameterData(dict={
    'SYSTEM': 'InAs',
    'EDIFF': 1e-5,
    'LORBIT': 11,
    'LSORBIT': True,
    'GGA_COMPAT': False,
    'ISMEAR': 0,
    'SIGMA': 0.05,
    'GGA': 'PE',
    'ENCUT': 280.00 * ase.units.eV,
    'MAGMOM': 6*[0.0],
    'NBANDS': 24,
})

# ---------- set up the Structure ----------
larray = array([[0,.5,.5],
                [.5,0,.5],
                [.5,.5,0]])
alat = 6.058
def carthesian_from_direct(structure, coords):
    rc = structure.get_ase().get_reciprocal_cell()
    return linalg.solve(rc, array(coords).T).T

def append_with_direct(structure, position, symbols):
    cart = carthesian_from_direct(structure, position)
    structure.append_atom(position=cart, symbols=symbols)

structure = StructureData(cell=larray*alat)
append_with_direct(structure, position=[0,0,0], symbols='In')
append_with_direct(structure, position=[.25,.25,.25], symbols='As')

# ---------- set up potentials ----------
#~ potcar = Potcar.from_file(sys.argv[2]).as_dict()
from os.path import abspath
#~ potcar = abspath(sys.argv[2])
potcar = ParameterData(dict={
    'potpaw_path': abspath(sys.argv[2]),
    'special_symbols' : {'In': '_d'},
})

# ---------- set kpoints info ----------
kpoints = ParameterData(dict={
    'comment': 'fcc InAs - selfconsistent',
    'generation_style': 'Gamma',
    'kpoints': (8,8,8),
    'usershift': [0,0,0],
})

code = Code.get_from_string('asevasp')
calc = code.new_calc()

# ---------- set inputs ----------
calc.use_settings(incar)
calc.use_structure(structure)
#~ calc.use_potentials(ParameterData(dict=potcar))
calc.use_potentials(potcar)
calc.use_kpoints(kpoints)

calc.label = 'Vasp Selfconsistent Test Run'
calc.description = 'Test the Vasp Plugin with a simple selfconsistent run for fcc InAs'
calc.set_max_wallclock_seconds(60)
calc.set_withmpi(False)
calc.set_resources({'num_machines':1, 'num_mpiprocs_per_machine':1})
calc.set_queue_name('dphys_compute')

if sys.argv[1] == '--submit':
    calc.store_all()
    print "created calculation; with uuid='{}' and PK={}".format(calc.uuid,calc.pk)
    calc.submit()
if sys.argv[1] == '--test':
    from aiida.common.folders import Folder
    calc.submit_test(Folder(sys.argv[3]))

