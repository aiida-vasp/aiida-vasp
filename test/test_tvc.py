from aiida import load_dbenv
load_dbenv()

from aiida.orm import CalculationFactory, DataFactory, Code, Computer
from aiida.common.folders import Folder
from aiida.orm.calculation.job.vasp import cif
import sys
import os

Tvc = CalculationFactory('vasp.base.TentativeVaspCalc')
Par = DataFactory('parameter')
Str = DataFactory('structure')
Paw = DataFactory('vasp.potpaw')
Kp = DataFactory('array.kpoints')

cifname = sys.argv[1]
cifnode = cif.cif_from_file(cifname)
equivalent_cifs = cif.filter_cifs_for_structure(
    cif.get_cifs_with_name(os.path.basename(cifname)),
    cifnode.get_ase())
if equivalent_cifs:
    cifnode = equivalent_cifs[0]

kp = Kp()
kp.set_kpoints_mesh([4,4,4])

paw_in = Paw.load_paw(family='LDA', symbol='In_d')[0]
paw_as = Paw.load_paw(family='LDA', symbol='As')[0]

tvc = Tvc()
tvc.use_code(Code.get_from_string('asevasp@monch'))

tvc.use_incar(Par(dict={'a': 200, 'b': True, 'c': 3*[6.]}))
tvc.use_poscar(cif.cif_to_structure(cifnode=cifnode))
tvc.use_potcar(paw_in, kind='In')
tvc.use_potcar(paw_as, kind='As')
tvc.use_kpcar(kp)

tvc.set_computer(Computer.get('monch'))
tvc.set_resources({
    'num_machines': 1,
    'num_mpiprocs_per_machine': 1})

tvc.submit_test(Folder('../submit_test'), subfolder_name='tvc')
print tvc.elements
