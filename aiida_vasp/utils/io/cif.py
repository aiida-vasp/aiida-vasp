from aiida.orm.calculation.inline import optional_inline
from aiida.orm.querytool import QueryTool
from aiida.orm import DataFactory
import os


def cif_from_file(ciffile):
    cifabs = os.path.abspath(ciffile)
    Cif = DataFactory('cif')
    cif = Cif()
    cif.set_file(cifabs)
    return cif


def cif_to_structure(cifnode=None):
    Structure = DataFactory('structure')
    structure = Structure()
    structure.set_ase(cifnode.get_ase())
    return structure


@optional_inline
def cif_to_structure_inline(cif=None):
    Structure = DataFactory('structure')
    structure = Structure()
    structure.set_ase(cif.ase)
    return {'structure': structure}


def get_or_create_structure(cif=None, use_first=False):
    cts = filter(cts_filter, cif.get_outputs())
    if len(cts) > 1 and not use_first:
        raise Exception('more than one cif->str calculations found')
    elif len(cts) == 1:
        return cts[0].out.structure
    else:
        return cif_to_structure_inline(cif=cif, store=True)['structure']


def cts_filter(node):
    fn = node.get_attr('function_name', None)
    return fn == 'cif_to_structure_inline'


def get_cifs_with_name(filename):
    q = QueryTool()
    q.set_class(DataFactory('cif'))
    q.add_attr_filter('filename', '=', filename)
    return q.run_query()


def filter_cifs_for_structure(cif_seq, structure):
    ase = None
    if hasattr(structure, 'get_ase'):
        ase = structure.get_ase()
    else:
        ase = structure
    return filter(lambda s: s.get_ase() == ase,
                  cif_seq)
