"""Test the POSCAR parser."""
# pylint: disable=unused-import,redefined-outer-name,unused-argument,unused-wildcard-import,wildcard-import, import-outside-toplevel
import pytest
import numpy as np

from aiida_vasp.utils.fixtures import *
from aiida_vasp.utils.fixtures.testdata import data_path
from aiida_vasp.utils.aiida_utils import aiida_version, cmp_version
from aiida_vasp.parsers.content_parsers.poscar import PoscarParser
from aiida_vasp.utils.aiida_utils import get_data_class


@pytest.mark.parametrize(['poscar_parser'], [('poscar',)], indirect=True)
def test_parse_poscar(poscar_parser):
    """Load a reference POSCAR parser.

    We check that it parses and provides the correct content.

    """

    # The structure for the POSCAR parser should have the key `poscar-structure`
    result = poscar_parser.get_quantity('poscar-structure')
    compare_poscar_content(result)


@pytest.mark.parametrize(['poscar_parser'], [('poscar',)], indirect=True)
def test_parse_poscar_write(poscar_parser, tmpdir):
    """Load a reference POSCAR parser and check that the write functionality works.

    Here we make sure the write function of the content parser works.

    """

    # Write the loaded structure to file
    temp_path = str(tmpdir.join('POSCAR'))
    poscar_parser.write(temp_path)

    # Load the written structure using a new content parser instance
    content = None
    with open(temp_path, 'r', encoding='utf8') as handler:
        content = handler.readlines()
    ref_content = [
        '# Compound: Al4. Old comment: Al\n', '  1.000000000000\n', '  4.040000000000   0.000000000000   0.000000000000\n',
        '  0.000000000000   4.040000000000   0.000000000000\n', '  0.000000000000   0.000000000000   4.040000000000\n', 'Al\n', '    4\n',
        'Direct\n', '  0.000000000000   0.000000000000   0.000000000000\n', '  0.000000000000   0.500000000000   0.500000000000\n',
        '  0.500000000000   0.000000000000   0.500000000000\n', '  0.500000000000   0.500000000000   0.000000000000\n', '\n',
        '  0.000000000000   0.000000000000   0.000000000000\n', '  0.000000000000   0.000000000000   0.000000000000\n',
        '  0.000000000000   0.000000000000   0.000000000000\n', '  0.000000000000   0.000000000000   0.000000000000\n'
    ]
    assert content == ref_content


@pytest.mark.parametrize(['vasp_structure'], [('str-Al',)], indirect=True)
def test_parse_poscar_data(vasp_structure, tmpdir):
    """Load a reference AiiDA StructureData and check that the parser can
    initialize using the data.

    Using the StructureData sitting in the initialized parser it should
    write that content to a POSCAR file when issuing write which is also tested,
    file is reloaded and content checked.

    """

    # Initialize parser with an existing reference StructureData
    poscar_parser = PoscarParser(data=vasp_structure)

    # Check that get_quantity return the same StructureData instance
    assert vasp_structure == poscar_parser.get_quantity('key_does_not_matter')

    # Write the loaded StructureData to file, which behind the scenes convert it
    # to a POSCAR format
    temp_path = str(tmpdir.join('POSCAR'))
    poscar_parser.write(temp_path)

    # Load the written structure using a new content parser instance
    parser = None
    with open(temp_path, 'r', encoding='utf8') as handler:
        parser = PoscarParser(handler=handler)
    result = parser.get_quantity('poscar-structure')
    # When we start from StructureData we do not have any velocity or predictor
    # values present, so let us not compare them. The reference POSCAR has no
    # velocities, but predictors, so we can not compare those directly.
    compare_poscar_content(result, predictors=False)


@pytest.mark.parametrize(['vasp_structure'], [('str-Al',)], indirect=True)
def test_consistency_with_parsevasp(vasp_structure):
    """Compare the poscar-dict returned by parsevasp to the dict created by the PoscarParser.

    This tests purpose is to give a warning if we are overriding keys in parsevasps poscar-dict.
    """
    from aiida_vasp.parsers.content_parsers.poscar import parsevasp_to_aiida
    from parsevasp.poscar import Poscar

    path = data_path('poscar', 'POSCAR')
    poscar = Poscar(file_path=path, prec=12, conserve_order=True)

    poscar_dict = poscar.get_dict(direct=False)
    result_dict = parsevasp_to_aiida(poscar)
    compare_objects(poscar_dict, result_dict)


def compare_objects(obj_a, obj_b):
    """Compare two potentially nested objects assuming they have the same structure."""
    if isinstance(obj_a, dict):
        for key in obj_a:
            compare_objects(obj_a[key], obj_b[key])
            return

    if isinstance(obj_a, (list, tuple, np.ndarray)):
        for item_a, item_b in zip(obj_a, obj_b):
            compare_objects(item_a, item_b)
            return

    assert obj_a == obj_b


def compare_poscar_content(result, velocities=True, predictors=True):
    """Compare the content of a supplied get_quantity('poscar-structure') with reference data."""
    unitcell = np.array([[4.04, 0., 0.], [0., 4.04, 0.], [0., 0., 4.04]])
    keys = ['comment', 'unitcell', 'sites']
    species = [{
        'specie': 'Al',
        'position': np.array([0., 0., 0.]),
        'selective': [True, True, True],
        'velocities': None,
        'predictors': np.array([0., 0., 0.]),
        'direct': False,
        'symbol': 'Al',
        'kind_name': 'Al'
    }, {
        'specie': 'Al',
        'position': np.array([0., 2.02, 2.02]),
        'selective': [True, True, True],
        'velocities': None,
        'predictors': np.array([0., 0., 0.]),
        'direct': False,
        'symbol': 'Al',
        'kind_name': 'Al'
    }, {
        'specie': 'Al',
        'position': np.array([2.02, 0., 2.02]),
        'selective': [True, True, True],
        'velocities': None,
        'predictors': np.array([0., 0., 0.]),
        'direct': False,
        'symbol': 'Al',
        'kind_name': 'Al'
    }, {
        'specie': 'Al',
        'position': np.array([2.02, 2.02, 0.]),
        'selective': [True, True, True],
        'velocities': None,
        'predictors': np.array([0., 0., 0.]),
        'direct': False,
        'symbol': 'Al',
        'kind_name': 'Al'
    }]
    assert set(result.keys()) == set(keys)
    assert np.allclose(result['unitcell'], unitcell)
    for index, item in enumerate(species):
        assert item['specie'] == result['sites'][index]['specie']
        assert item['selective'] == result['sites'][index]['selective']
        assert item['direct'] == result['sites'][index]['direct']
        assert item['symbol'] == result['sites'][index]['symbol']
        assert item['kind_name'] == result['sites'][index]['kind_name']
        assert np.allclose(item['position'], result['sites'][index]['position'])
        if velocities:
            assert item['velocities'] == result['sites'][index]['velocities']
        if predictors:
            assert np.allclose(item['predictors'], result['sites'][index]['predictors'])
