"""Unittests for PoscarIo"""
# pylint: disable=unused-import,redefined-outer-name,unused-argument,unused-wildcard-import,wildcard-import
import pytest
import sys

from aiida_vasp.utils.fixtures import *
from aiida_vasp.utils.fixtures.testdata import data_path
from aiida_vasp.io.poscar import PoscarParser

@pytest.mark.parametrize(['vasp_structure'], [('str',)], indirect=True)
def test_poscar_io(fresh_aiida_env, vasp_structure_poscar):
    poscario = vasp_structure_poscar
    assert poscario.potentials_order == ['In', 'As', 'In_d', 'As']

@pytest.mark.parametrize(['vasp_structure'], [('str-Al',)], indirect=True)
def test_parse_poscar(fresh_aiida_env, vasp_structure):
    """Parse a reference POSCAR file with the PoscarParser and 
    compare the result to a reference structure.

    """
    
    path = data_path('poscar', 'POSCAR')
    parser = PoscarParser(file_path=path)
    result = parser.get_quantity('structure', {})
    structure = vasp_structure
    assert result['structure'].cell == structure.cell
    assert result['structure'].get_site_kindnames() == structure.get_site_kindnames()
    assert result['structure'].sites[2].position == structure.sites[2].position
    
@pytest.mark.parametrize(['vasp_structure'], [('str-Al',)], indirect=True)
def test_parse_poscar_write(fresh_aiida_env, vasp_structure, tmpdir):
    """Parse (write) a reference POSCAR file with the PoscarParser, 
    parse (read), and compare to reference structure.

    """

    structure = vasp_structure
    parser = PoscarParser(data = structure)
    result = parser.get_quantity('structure', {})
    assert result['structure'].cell == structure.cell
    assert result['structure'].get_site_kindnames() == structure.get_site_kindnames()
    assert result['structure'].sites[2].position == structure.sites[2].position
    temp_file = str(tmpdir.join('POSCAR'))
    parser.write(temp_file)
    parser = PoscarParser(file_path=temp_file)
    result_reparse = parser.get_quantity('structure', {})
    assert result_reparse['structure'].cell == structure.cell
    assert result_reparse['structure'].get_site_kindnames() == \
        structure.get_site_kindnames()
    assert result_reparse['structure'].sites[2].position == \
        structure.sites[2].position

@pytest.mark.parametrize(['vasp_structure'], [('str-Al',)], indirect=True)
def test_parse_poscar_reparse(fresh_aiida_env, vasp_structure, tmpdir):
    """Parse (read) a reference POSCAR file with the PoscarParser, parse(write), 
    parse (read), and compare to reference structure.

    """
    
    path = data_path('poscar', 'POSCAR')
    parser = PoscarParser(file_path=path)
    result = parser.get_quantity('structure', {})
    temp_file = str(tmpdir.join('POSCAR'))
    parser.write(temp_file)
    parser = PoscarParser(file_path=temp_file)
    result_reparse = parser.get_quantity('structure', {})
    structure = vasp_structure
    assert result_reparse['structure'].cell == structure.cell
    assert result_reparse['structure'].get_site_kindnames() == \
        structure.get_site_kindnames()
    assert result_reparse['structure'].sites[2].position == \
        structure.sites[2].position

# The following tests can be enabled when Aiida gets updated
# (hopefully) with an element X.
    
# def test_parse_poscar_silly_read(fresh_aiida_env):
#     """Parse a reference POSCAR with silly elemental names with the 
#     PoscarParser and compare the result to a reference structure.

#     """

#     path = data_path('poscar', 'POSCARSILLY')
#     parser = PoscarParser(file_path=path)
#     result = parser.get_quantity('structure', {})
#     names = result['structure'].get_site_kindnames()
#     assert names == ['Hamburger', 'Pizza']
#     symbols = result['structure'].get_symbols_set()
#     assert symbols == set(['X', 'X'])
    
# @pytest.mark.parametrize(['vasp_structure'], [('str-InAs',)], indirect=True)
# def test_parse_poscar_silly_write(fresh_aiida_env, vasp_structure, tmpdir):
#     """Parse a reference POSCAR with silly elemental names with the 
#     PoscarParser and compare the result to a reference structure.

#     """

#     parser = PoscarParser(data=vasp_structure)
#     result = parser.get_quantity('structure', {})
#     names = result['structure'].get_site_kindnames()
#     assert names == ['Hamburger', 'Pizza']
#     symbols = result['structure'].get_symbols_set()
#     assert symbols == set(['As', 'In'])
#     temp_file = str(tmpdir.join('POSCAR'))
#     parser.write(temp_file)
#     parser = PoscarParser(file_path=temp_file)
#     result_reparse = parser.get_quantity('structure', {})
#     names = result_reparse['structure'].get_site_kindnames()
#     assert names == ['Hamburger', 'Pizza']
#     symbols = result_reparse['structure'].get_symbols_set()
#     assert symbols == set(['X', 'X'])    
    

@pytest.mark.parametrize(['vasp_structure'], [('str',)], indirect=True)
def test_parse_poscar_undercase(fresh_aiida_env, vasp_structure, tmpdir):
    """Parse a reference POSCAR with potential elemental names with the 
    PoscarParser and compare the result to a reference structure.

    """

    parser = PoscarParser(data=vasp_structure)
    result = parser.get_quantity('structure', {})
    names = result['structure'].get_site_kindnames()
    assert names == ['In', 'As', 'As', 'In_d', 'In_d', 'As']
    symbols = result['structure'].get_symbols_set()
    assert symbols == set(['As', 'In'])
    temp_file = str(tmpdir.join('POSCAR'))
    parser.write(temp_file)
    parser = PoscarParser(file_path=temp_file)
    result_reparse = parser.get_quantity('structure', {})
    names = result_reparse['structure'].get_site_kindnames()
    assert names == ['In', 'As', 'As', 'In_d', 'In_d', 'As']
    symbols = result_reparse['structure'].get_symbols_set()
    assert symbols == set(['As', 'In'])
