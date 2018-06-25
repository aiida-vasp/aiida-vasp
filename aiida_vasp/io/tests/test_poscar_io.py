"""Unittests for PoscarIo"""
# pylint: disable=unused-import,redefined-outer-name,unused-argument,unused-wildcard-import,wildcard-import
import pytest

from aiida_vasp.utils.fixtures import *
from aiida_vasp.utils.fixtures.environment import aiida_version
from aiida_vasp.utils.fixtures.testdata import data_path
from aiida_vasp.io.poscar import PoscarParser


@pytest.mark.parametrize(['vasp_structure'], [('str-Al',)], indirect=True)
def test_parse_poscar(fresh_aiida_env, vasp_structure):
    """
    Parse a reference POSCAR file.

    Using the PoscarParser and compare the result to a reference
    structure.

    """

    path = data_path('poscar', 'POSCAR')
    parser = PoscarParser(file_path=path)
    result = parser.get_quantity('poscar-structure', {})
    structure = vasp_structure
    assert result['poscar-structure'].cell == structure.cell
    assert result['poscar-structure'].get_site_kindnames() == structure.get_site_kindnames()
    assert result['poscar-structure'].sites[2].position == structure.sites[2].position


@pytest.mark.parametrize(['vasp_structure'], [('str-Al',)], indirect=True)
def test_parse_poscar_write(fresh_aiida_env, vasp_structure, tmpdir):
    """
    Parse (write) a reference POSCAR file.

    Using the PoscarParser, parse (read), and compare to reference
    structure.

    """

    structure = vasp_structure
    parser = PoscarParser(data=structure)
    result = parser.get_quantity('poscar-structure', {})
    assert result['poscar-structure'].cell == structure.cell
    assert result['poscar-structure'].get_site_kindnames() == structure.get_site_kindnames()
    assert result['poscar-structure'].sites[2].position == structure.sites[2].position
    temp_file = str(tmpdir.join('POSCAR'))
    parser.write(temp_file)
    parser = PoscarParser(file_path=temp_file)
    result_reparse = parser.get_quantity('poscar-structure', {})
    assert result_reparse['poscar-structure'].cell == structure.cell
    assert result_reparse['poscar-structure'].get_site_kindnames() == \
        structure.get_site_kindnames()
    assert result_reparse['poscar-structure'].sites[2].position == \
        structure.sites[2].position


@pytest.mark.parametrize(['vasp_structure'], [('str-Al',)], indirect=True)
def test_parse_poscar_reparse(fresh_aiida_env, vasp_structure, tmpdir):
    """
    Parse (read) a reference POSCAR file.

    Using the PoscarParser, parse(write), parse (read), and compare
    to reference structure.

    """

    path = data_path('poscar', 'POSCAR')
    parser = PoscarParser(file_path=path)
    _ = parser.get_quantity('poscar-structure', {})
    temp_file = str(tmpdir.join('POSCAR'))
    parser.write(temp_file)
    parser = PoscarParser(file_path=temp_file)
    result_reparse = parser.get_quantity('poscar-structure', {})
    structure = vasp_structure
    assert result_reparse['poscar-structure'].cell == structure.cell
    assert result_reparse['poscar-structure'].get_site_kindnames() == \
        structure.get_site_kindnames()
    assert result_reparse['poscar-structure'].sites[2].position == \
        structure.sites[2].position


@pytest.mark.xfail(aiida_version() < (1, 0), reason="Element X only present in Aiida >= 1.x")
def test_parse_poscar_silly_read(fresh_aiida_env):
    """
    Parse (read) a reference POSCAR with silly elemental names.

    Using the PoscarParser and compare the result to a reference
    structure.

    """

    path = data_path('poscar', 'POSCARSILLY')
    parser = PoscarParser(file_path=path)
    result = parser.get_quantity('poscar-structure', {})
    names = result['poscar-structure'].get_site_kindnames()
    assert names == ['Hamburger', 'Pizza']
    symbols = result['poscar-structure'].get_symbols_set()
    assert symbols == set(['X', 'X'])


@pytest.mark.xfail(aiida_version() < (1, 0), reason="Element X only present in Aiida >= 1.x")
@pytest.mark.parametrize(['vasp_structure'], [('str-InAs',)], indirect=True)
def test_parse_poscar_silly_write(fresh_aiida_env, vasp_structure, tmpdir):
    """
    Parse (read, write) a reference POSCAR with silly elemental names.

    Using the PoscarParser and compare the result to a reference structure.

    """

    parser = PoscarParser(data=vasp_structure)
    result = parser.get_quantity('poscar-structure', {})
    names = result['poscar-structure'].get_site_kindnames()
    assert names == ['Hamburger', 'Pizza']
    symbols = result['poscar-structure'].get_symbols_set()
    assert symbols == set(['As', 'In'])
    temp_file = str(tmpdir.join('POSCAR'))
    parser.write(temp_file)
    parser = PoscarParser(file_path=temp_file)
    result_reparse = parser.get_quantity('poscar-structure', {})
    names = result_reparse['poscar-structure'].get_site_kindnames()
    assert names == ['Hamburger', 'Pizza']
    symbols = result_reparse['poscar-structure'].get_symbols_set()
    assert symbols == set(['X', 'X'])


@pytest.mark.parametrize(['vasp_structure'], [('str',)], indirect=True)
def test_parse_poscar_undercase(fresh_aiida_env, vasp_structure, tmpdir):
    """
    Parse a reference POSCAR.

    With potential elemental names using the PoscarParser and compare
    the result to a reference structure.

    """

    parser = PoscarParser(data=vasp_structure)
    result = parser.get_quantity('poscar-structure', {})
    names = result['poscar-structure'].get_site_kindnames()
    assert names == ['In', 'As', 'As', 'In_d', 'In_d', 'As']
    symbols = result['poscar-structure'].get_symbols_set()
    assert symbols == set(['As', 'In'])
    temp_file = str(tmpdir.join('POSCAR'))
    parser.write(temp_file)
    parser = PoscarParser(file_path=temp_file)
    result_reparse = parser.get_quantity('poscar-structure', {})
    names = result_reparse['poscar-structure'].get_site_kindnames()
    assert names == ['In', 'As', 'As', 'In_d', 'In_d', 'As']
    symbols = result_reparse['poscar-structure'].get_symbols_set()
    assert symbols == set(['As', 'In'])
