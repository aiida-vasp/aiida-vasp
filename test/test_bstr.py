__doc__ = '''starts a bandstructure workflow and runs it on monch'''
from aiida.orm import DataFactory, WorkflowFactory
from numpy import array, linalg
import sys

BStr = WorkflowFactory('vasp.bandstructure.Bandstructure')
PDat = DataFactory('parameter')
SDat = DataFactory('structure')

# ---------- set up incar ----------
incar = PDat(dict={
    'system': 'InAs',
    'ediff': 1e-5,
    'lorbit': 11,
    'lsorbit': True,
    'gga_compat': False,
    'ismear': 0,
    'sigma': .05,
    'gga': 'PE',
    'encut': 280,
    'magmom': 6*[0.0],
    'nbands': 24
})

# ---------- set up structure ----------
def carthesian_from_direct(structure, coords):
    rc = structure.get_ase().get_reciprocal_cell()
    return linalg.solve(rc, array(coords).T).T

def append_with_direct(structure, position, symbols):
    cart = carthesian_from_direct(structure, position)
    structure.append_atom(position=cart, symbols=symbols)

larray = array([[0,.5,.5],[.5,0,.5],[.5,.5,0]])
alat = 6.058
structure = SDat(cell=larray*alat)
append_with_direct(structure, position=3*[0], symbols='In')
append_with_direct(structure, position=3*[.25], symbols='As')
strucdict = {'cell': larray*alat,
             'atoms': [{'position': carthesian_from_direct(structure, 3*[0]),
                       'symbols': 'In'},
                      {'position': carthesian_from_direct(structure, 3*[.25]),
                       'symbols': 'As'}]
             }

# ---------- set up potentials ----------
from os.path import abspath
potcar = PDat(dict={
    'potpaw_path': abspath(sys.argv[1]),
    'special_symbols': {'In': '_d'},
})

# ---------- set kpoints info ----------
kpoints_sc = PDat(dict={
    'generation_style': 'Gamma',
    'kpoints': 3*[8],
    'usershift': 3*[0]
})

kpoints_bs = PDat(dict={
    'reciprocal': True,
    'intersections': 50,
    'kpoints': [
        [1.,0.,0.],
        [0.,0.,0.],
#
        [0.,0.,0.],
        [.5,.5,.5]
    ]
})

# ---------- store input nodes ----------
incar.store()
structure.store()
potcar.store()
kpoints_sc.store()
kpoints_bs.store()

# ---------- setup workflow ----------
wfpar = {'incar': incar.pk,
         'structure': structure.pk,
         'potentials': potcar.pk,
         'sc_kpoints': kpoints_sc.pk,
         'bs_kpoints': kpoints_bs.pk,
         'name': 'bstr_workflow_test',
         'queue': 'dphys_compute',
         }
wofl = BStr(params=wfpar)
wofl.start()
