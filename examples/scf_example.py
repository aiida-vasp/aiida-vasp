#!runaiida
from aiida.orm import Code, Group
from aiida_vasp.calcs.maker import VaspMaker
from aiida.common.exceptions import NotExistent
import sys

usage = '''runaiida scf_example.py <code@comp> <cif/POSCAR file>'''


def add_to_group(node):
    try:
        g = Group.get_from_string('examples')
    except NotExistent:
        g = Group(name='examples')
        g.store()
    g.add_nodes(node)

if __name__ == '__main__':
    if not len(sys.argv) == 3:
        print usage
        sys.exit()
    mkr = VaspMaker(
        structure=sys.argv[2],
        calc_cls='vasp.scf',
        paw_map={'Ga': 'Ga', 'As': 'As'},
        paw_family='pbe',
    )
    mkr.add_settings(
        gga='PE',
        gga_compat=False,
        ediff=1e-5,
        lorbit=11,
        ismear=0,
        sigma=.05,
        encut=280
    )
    mkr.set_kpoints_mesh([2, 2, 2])
    mkr.queue = 'dphys_compute'
    mkr.resources['num_machines'] = 2
    mkr.resources['num_mpiprocs_per_machine'] = 4
    mkr.code = Code.get_from_string(sys.argv[1])
    mkr.label = 'scf_example'
    calc = mkr.new()
    calc.store_all()
    calc.submit()
    add_to_group(calc)
    print 'submitted SCF : PK = {pk}'.format(pk=calc.pk)
