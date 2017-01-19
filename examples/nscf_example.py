#!runaiida
from aiida.orm import Code, Group, load_node
from aiida_vasp.calcs.maker import VaspMaker
from aiida.common.exceptions import NotExistent
import sys

usage = '''runaiida scf_example.py <code@comp> <scf pk>
    scf pk can be chosen from the following list (if you've run a scf calc before:
'''

def get_group():
    try:
        g = Group.get_from_string('examples')
    except NotExistent:
        g = Group(name='examples')
        g.store()
    return g

def add_to_group(node):
    g = get_group()
    g.add_nodes(node)

def has_scf_output(node):
    out = node.get_outputs_dict()
    search_for = {
        'charge_density',
        'wavefunctions'
    }
    return search_for.issubset(out)

def get_eligible_scf():
    from aiida.orm.querytool import QueryTool
    from aiida.orm import Calculation
    q = QueryTool()
    q.set_class(Calculation)
    return filter(has_scf_output, q.run_query())

def list_item(node=None):
    it = '\t{pk} {label} {structure}'
    if node == None:
        return it.format(pk='PK', label='Label', structure='Structure')
    return it.format(pk=node.pk, label=node.label,
                     structure=node.inp.structure.get_ase().get_chemical_formula())

scflist = '\n'.join(map(list_item, get_eligible_scf()))

if __name__ == '__main__':
    if not len(sys.argv) == 3:
        print(usage)
        print(list_item())
        print(scflist)
        sys.exit()
    mkr = VaspMaker(continue_from=load_node(pk=sys.argv[2]), calc_cls='vasp.nscf')
    mkr.add_settings(icharg=11)
    mkr.set_kpoints_path()
    mkr.resources['num_machines'] = 2
    mkr.resources['num_mpiprocs_per_machine'] = 4
    mkr.code = Code.get_from_string(sys.argv[1])
    mkr.label = 'nscf_example'
    calc = mkr.new()
    calc.store_all()
    calc.submit()
    add_to_group(calc)
    print 'submitted NSCF : PK = {pk}'.format(pk=calc.pk)
