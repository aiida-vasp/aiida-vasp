"""Test TentativeVaspCalculation: test used in development, this is now outdated."""

import sys
import os

from aiida import load_dbenv
load_dbenv()

from aiida.orm import CalculationFactory, DataFactory, Code, Computer
from aiida.orm.querytool import QueryTool
from aiida.common.folders import Folder
from aiida_vasp.utils.io import cif
from aiida_vasp.utils.default_paws import DEFAULT_LDA

Tvc = CalculationFactory('vasp.base.TentativeVaspCalc')
Par = DataFactory('parameter')
Str = DataFactory('structure')
Paw = DataFactory('vasp.potpaw')
Kp = DataFactory('array.kpoints')

cifname = sys.argv[1]
cifnode = cif.cif_from_file(cifname)
equivalent_cifs = cif.filter_cifs_for_structure(
    cif.get_cifs_with_name(os.path.basename(cifname)), cifnode.get_ase())
if equivalent_cifs:
    cifnode = equivalent_cifs[0]

structure = cif.cif_to_structure(cifnode=cifnode)

kp = Kp()
kp.set_cell_from_structure(structure)
kp.set_kpoints_mesh([4, 4, 4])

tvc = Tvc()
tvc.use_code(Code.get_from_string('asevasp@monch'))

nions = len(structure.sites)
nbands = 9 * nions + 8
if nbands % 8:
    nbands += 8 - (nbands % 8)
tvc.use_incar(
    Par(dict={
        'gga': 'PE',
        'gga_compat': False,
        'encut': 280,
        'ediff': 1e-5,
        'ismear': 0,
        'lorbit': 11,
        'magmom': 3 * nions * [0.],
        'lsorbit': True,
        'nbands': nbands,
        'sigma': 0.05,
        'system': os.path.basename(cifname),
        'npar': 8,
    }))
tvc.use_structure(structure)

for kind in structure.get_kind_names():
    tvc.use_paw(Paw.load_paw(family='LDA', symbol=kind)[0], kind=kind)

tvc.use_kpoints(kp)

tvc.set_computer(Computer.get('monch'))
tvc.set_queue_name('dphys_compute')
tvc.set_resources({'num_machines': nbands / 8, 'num_mpiprocs_per_machine': 8})
tvc.set_max_memory_kb(16000000)

tvc.set_parser_name('vasp.vasp5')

tag = sys.argv[2]
q = QueryTool()
q.set_class(Tvc)
q.add_extra_filter('tag', '=', tag)
ql = map(lambda c: c.get_extra('test-nr'), q.run_query())

last_tn = max(ql)

tvc.store_all()

tvc.set_extras({'tag': tag, 'test-nr': last_tn + 1})
print 'pk: {}, test-nr: {}'.format(tvc.pk, last_tn + 1)
tvc.submit()
