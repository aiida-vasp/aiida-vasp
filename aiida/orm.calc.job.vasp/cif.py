from aiida.orm import DataFactory, JobCalculation
from aiida.orm.data.structure import has_ase
from aiida.orm.calculation.inline import optional_inline
from aiida.orm.querytool import QueryTool
import ase
import os

def cif_from_file(ciffile):
    cifabs = os.path.abspath(ciffile)
    Cif = DataFactory('cif')
    cif = Cif()
    cif.set_file(cifabs)
    return cif

@optional_inline
def cif_to_structure(cifnode=None):
    Structure = DataFactory('structure')
    structure = Structure()
    structure.set_ase(cifnode.get_ase())
    return structure

def get_cifs_with_name(filename):
    q = QueryTool()
    q.set_class(DataFactory('cif'))
    q.add_attr_filter('filename', '=', filename)
    return q.run_query()

def filter_cifs_for_structure(cif_seq, structure):
    Str = DataFactory('structure')
    Cif = DataFactory('cif')
    ase = None
    if hasattr(structure, 'get_ase'):
        ase = structure.get_ase()
    else:
        ase = structure
    return filter(
            lambda s : s.get_ase() == ase,
            cif_seq)
