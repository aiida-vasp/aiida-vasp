"""Test the POSCAR parser."""
# pylint: disable=unused-import,redefined-outer-name,unused-argument,unused-wildcard-import,wildcard-import, import-outside-toplevel
import pytest

from aiida_vasp.utils.fixtures import *
from aiida_vasp.utils.fixtures.testdata import data_path
from aiida_vasp.utils.aiida_utils import aiida_version, cmp_version
from aiida_vasp.parsers.file_parsers.poscar import PoscarParser


@pytest.mark.parametrize(['vasp_structure'], [('str-Al',)], indirect=True)
def test_parse_poscar(fresh_aiida_env, vasp_structure):
    """
    Parse a reference POSCAR file.

    Using the PoscarParser and compare the result to a reference
    structure.

    """

    path = data_path('poscar', 'POSCAR')
    parser = PoscarParser(file_path=path)
    result = parser.structure
    structure = vasp_structure

    assert result.cell == structure.cell
    assert result.get_site_kindnames() == structure.get_site_kindnames()
    assert result.sites[2].position == structure.sites[2].position


@pytest.mark.parametrize(['vasp_structure'], [('str-Al',)], indirect=True)
def test_parse_poscar_write(fresh_aiida_env, vasp_structure, tmpdir):
    """
    Parse (write) a reference POSCAR file.

    Using the PoscarParser, and compare to reference
    structure.

    """

    structure = vasp_structure
    parser = PoscarParser(data=structure)

    temp_file = str(tmpdir.join('POSCAR'))
    parser.write(temp_file)

    parser = PoscarParser(file_path=temp_file)
    result_reparse = parser.structure

    assert result_reparse.cell == structure.cell
    assert result_reparse.get_site_kindnames() == \
        structure.get_site_kindnames()
    assert result_reparse.sites[2].position == \
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

    temp_file = str(tmpdir.join('POSCAR'))
    parser.write(temp_file)

    parser = PoscarParser(file_path=temp_file)
    result = parser.structure

    structure = vasp_structure
    assert result.cell == structure.cell
    assert result.get_site_kindnames() == structure.get_site_kindnames()
    assert result.sites[2].position == structure.sites[2].position


@pytest.mark.xfail(aiida_version() < cmp_version('1.0.0a1'), reason='Element X only present in Aiida >= 1.x')
def test_parse_poscar_silly_read(fresh_aiida_env):
    """
    Parse (read) a reference POSCAR with silly elemental names.

    Using the PoscarParser and compare the result to a reference
    structure.

    """

    path = data_path('poscar', 'POSCARSILLY')
    parser = PoscarParser(file_path=path)
    result = parser.structure

    names = result.get_site_kindnames()
    assert names == ['Hamburger', 'Pizza']
    symbols = result.get_symbols_set()
    assert symbols == set(['X', 'X'])


@pytest.mark.xfail(aiida_version() < cmp_version('1.0.0a1'), reason='Element X only present in Aiida >= 1.x')
@pytest.mark.parametrize(['vasp_structure'], [('str-InAs',)], indirect=True)
def test_parse_poscar_silly_write(fresh_aiida_env, vasp_structure, tmpdir):
    """
    Parse (read, write) a reference POSCAR with silly elemental names.

    Using the PoscarParser and compare the result to a reference structure.

    """

    parser = PoscarParser(data=vasp_structure)
    result = parser.get_quantity('poscar-structure')
    names = result.get_site_kindnames()
    assert names == ['Hamburger', 'Pizza']
    symbols = result.get_symbols_set()
    assert symbols == set(['As', 'In'])

    temp_file = str(tmpdir.join('POSCAR'))
    parser.write(temp_file)

    parser = PoscarParser(file_path=temp_file)
    result_reparse = parser.structure

    names = result_reparse.get_site_kindnames()
    assert names == ['Hamburger', 'Pizza']
    symbols = result_reparse.get_symbols_set()
    assert symbols == set(['X', 'X'])


@pytest.mark.parametrize(['vasp_structure'], [('str',)], indirect=True)
def test_parse_poscar_undercase(fresh_aiida_env, vasp_structure, tmpdir):
    """
    Parse a reference POSCAR.

    With potential elemental names using the PoscarParser and compare
    the result to a reference structure.

    """

    parser = PoscarParser(data=vasp_structure)
    result = parser.get_quantity('poscar-structure')
    names = result.get_site_kindnames()
    assert names == ['In', 'As', 'As', 'In_d', 'In_d', 'As']
    symbols = result.get_symbols_set()
    assert symbols == set(['As', 'In'])
    temp_file = str(tmpdir.join('POSCAR'))
    parser.write(temp_file)
    parser = PoscarParser(file_path=temp_file)
    result_reparse = parser.structure
    names = result_reparse.get_site_kindnames()
    assert names == ['In', 'As', 'As', 'In_d', 'In_d', 'As']
    symbols = result_reparse.get_symbols_set()
    assert symbols == set(['As', 'In'])


@pytest.mark.parametrize(['vasp_structure'], [('str-Al',)], indirect=True)
def test_consistency_with_parsevasp(fresh_aiida_env, vasp_structure):
    """
    Compare the poscar-dict returned by parsevasp to the dict created by the PoscarParser.

    This tests purpose is to give a warning if we are overriding keys in parsevasps poscar-dict.
    """
    from aiida_vasp.parsers.file_parsers.poscar import parsevasp_to_aiida
    from parsevasp.poscar import Poscar

    path = data_path('poscar', 'POSCAR')
    poscar = Poscar(file_path=path, prec=12, conserve_order=True)

    poscar_dict = poscar.get_dict(direct=False)
    result_dict = parsevasp_to_aiida(poscar)['poscar-structure']
    compare_objects(poscar_dict, result_dict)


def compare_objects(obj_a, obj_b):
    """Compare two potentially nested objects assuming they have the same structure."""
    import numpy as np

    if isinstance(obj_a, dict):
        for key in obj_a:
            compare_objects(obj_a[key], obj_b[key])
            return

    if isinstance(obj_a, (list, tuple, np.ndarray)):
        for item_a, item_b in zip(obj_a, obj_b):
            compare_objects(item_a, item_b)
            return

    assert obj_a == obj_b
