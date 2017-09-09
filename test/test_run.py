__doc__ = '''sets up a selfconsistent vasp calculation and tries to run it on monch'''
from aiida.orm import DataFactory, Code
from pymatgen import Lattice
from pymatgen.io.vasp import Lattice, Structure, Poscar, Potcar, Kpoints
from numpy import array
import sys

incar = {
    'SYSTEM': 'InAs',
    'EDIFF': 1e-5,
    'LORBIT': 11,
    'LSORBIT': '.True.',
    'GGA_COMPAT': '.False.',
    'ISMEAR': 0,
    'SIGMA': 0.05,
    'GGA': 'PE',
    'ENCUT': '280.00 eV',
    'MAGMOM': '6*0.0',
    'NBANDS': 24,
}

larray = array([[0, .5, .5], [.5, 0, .5], [.5, .5, 0]])
lattice = Lattice(larray * 6.058)
species = ['In', 'As']
coords = [[0, 0, 0], [.25, .25, .25]]
structure = Structure(lattice, species, coords, coords_are_cartesian=False)
poscar = Poscar(structure, 'fcc InAs').as_dict()

#~ potcar = Potcar.from_file(sys.argv[2]).as_dict()
from os.path import abspath
potcar = abspath(sys.argv[2])

kpoints = {
    'comment': 'fcc InAs - selfconsistent',
    'generation_style': 'Gamma',
    'kpoints': [[8, 8, 8]],
    'usershift': [0, 0, 0],
}

code = Code.get_from_string('vasp')
calc = code.new_calc()

ParameterData = DataFactory('parameter')
SinglefileData = DataFactory('singlefile')
calc.use_parameters(ParameterData(dict=incar))
calc.use_structure(ParameterData(dict=poscar))
#~ calc.use_potentials(ParameterData(dict=potcar))
calc.use_potentials(SinglefileData(file=potcar))
calc.use_kpoints(ParameterData(dict=kpoints))

calc.label = 'Vasp Selfconsistent Test Run'
calc.description = 'Test the Vasp Plugin with a simple selfconsistent run for fcc InAs'
calc.set_max_wallclock_seconds(60)
calc.set_withmpi(False)
calc.set_resources({'num_machines': 1, 'num_mpiprocs_per_machine': 1})
calc.set_queue_name('dphys_compute')

if sys.argv[1] == '--submit':
    calc.store_all()
    print "created calculation; with uuid='{}' and PK={}".format(
        calc.uuid, calc.pk)
    calc.submit()
if sys.argv[1] == '--test':
    from aiida.common.folders import Folder
    calc.submit_test(Folder(sys.argv[3]))
