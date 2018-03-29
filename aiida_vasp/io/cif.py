#"""Cif data utilities"""
#import os

#from aiida.orm.calculation.inline import optional_inline
#from aiida.orm.querybuilder import QueryBuilder
#from aiida.orm import DataFactory

#def cif_from_file(cif_file):
#    """Create a CifData node from a cif file"""
#    cifabs = os.path.abspath(cif_file)
#    cif_cls = DataFactory('cif')
#    cif = cif_cls()
#    cif.set_file(cifabs)
#    return cif

#def cif_to_structure(cifnode=None):
#    structure_cls = DataFactory('structure')
#    structure = structure_cls()
#    structure.set_ase(cifnode.get_ase())
#    return structure

#@optional_inline
#def cif_to_structure_inline(cif=None):
#    """Create a structure from a CIF file in an inline calculation."""
#    structure_cls = DataFactory('structure')
#    structure = structure_cls()
#    structure.set_ase(cif.ase)
#    return {'structure': structure}

#def get_or_create_structure(cif=None, use_first=False):
#    """Find a structure for a cif file in the DB or else create (store) a new one."""
#    cts = filter(cts_filter, cif.get_outputs())
#    if len(cts) > 1 and not use_first:
#        raise Exception('more than one cif->str calculations found')
#    elif len(cts) == 1:
#        return cts[0].out.structure
#    else:
#        return cif_to_structure_inline(cif=cif, store=True)['structure']  # pylint: disable=unexpected-keyword-arg

#def cts_filter(node):
#    function_name = node.get_attr('function_name', None)
#    return function_name == 'cif_to_structure_inline'

#def get_cifs_with_name(filename):
#    query_builder = QueryBuilder()
#    query_builder.append(tag='node', cls=DataFactory('cif'))
#    query_builder.add_filter('node', {'filename': filename})
#    return query_builder.all()

#def filter_cifs_for_structure(cif_seq, structure):
#    """Return all cif files in the sequence which match the given structure."""
#    ase = None
#    if hasattr(structure, 'get_ase'):
#        ase = structure.get_ase()
#    else:
#        ase = structure
#    return [s for s in cif_seq if s.get_ase() == ase]
