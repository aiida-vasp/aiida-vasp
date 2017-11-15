"""Unittests for the PymatgenParser."""
# pylint: disable=unused-import,redefined-outer-name,unused-argument,unused-wildcard-import,wildcard-import

import numpy
import pytest
from aiida.common.exceptions import ParsingError

from aiida_vasp.parsers.pymatgen_vasp import PymatgenParser
from aiida_vasp.utils.fixtures.testdata import data_path
from aiida_vasp.utils.fixtures import *


@pytest.fixture(params=[-1])
def vasprun_path(request, tmpdir):
    """Truncate vasprun.xml at the given line number and parse."""
    original_path = data_path('phonondb', 'vasprun.xml')
    if request.param == -1:
        return original_path
    truncated_path = tmpdir.join('vasprun.xml')
    with open(original_path, 'r') as original_fo:
        truncated_content = '\n'.join(original_fo.readlines()[:request.param])
    truncated_path.write(truncated_content)
    return str(truncated_path)


@pytest.fixture(params=[None])
def parse_result(request, aiida_env, vasprun_path):
    """
    Give the result of parsing a retrieved calculation (emulated).

    Returns a function which does:

    1. create a calculation with parser settings
    2. update the parser settings with the extra_settings
    3. create a parser with the calculation
    4. populate a fake retrieved folder and pass it to the parser
    5. return the result
    """

    def parse(**extra_settings):
        """Run the parser using default settings updated with extra_settings."""
        from aiida.orm import CalculationFactory, DataFactory
        calc = CalculationFactory('vasp.vasp')()
        settings_dict = {'pymatgen_parser': {'parse_potcar_file': False, 'exception_on_bad_xml': request.param}}
        settings_dict.update(extra_settings)
        calc.use_settings(DataFactory('parameter')(dict=settings_dict))
        parser = PymatgenParser(calc=calc)
        retrieved = DataFactory('folder')()
        retrieved.add_path(vasprun_path, '')

        success, nodes = parser.parse_with_retrieved({'retrieved': retrieved})
        nodes = dict(nodes)
        return success, nodes

    return parse


@pytest.fixture()
def parse_nac(aiida_env):
    """Give the parsing result of a retrieved NAC calculation (emulated)."""
    from aiida.orm import CalculationFactory, DataFactory
    calc = CalculationFactory('vasp.vasp')()
    calc.use_settings(DataFactory('parameter')(dict={'pymatgen_parser': {'parse_potcar_file': False, 'exception_on_bad_xml': False}}))
    parser = PymatgenParser(calc=calc)
    retrieved = DataFactory('folder')()
    retrieved.add_path(data_path('born_effective_charge', 'vasprun.xml'), '')
    retrieved.add_path(data_path('born_effective_charge', 'OUTCAR'), '')

    def parse():
        success, nodes = parser.parse_with_retrieved({'retrieved': retrieved})
        nodes = dict(nodes)
        return success, nodes

    return parse


def test_kpoints_result(parse_result):
    """Test that the kpoints result node is a KpointsData instance."""
    from aiida.orm import DataFactory
    _, nodes = parse_result()
    assert isinstance(nodes['kpoints'], DataFactory('array.kpoints'))


def test_structure_result(parse_result):
    """Test that the structure result node is a StructureData instance."""
    from aiida.orm import DataFactory
    _, nodes = parse_result()
    assert isinstance(nodes['structure'], DataFactory('structure'))


def test_forces_result(parse_result):
    """Check the parsed forces result node."""
    from aiida.orm import DataFactory
    _, nodes = parse_result()
    assert isinstance(nodes['forces'], DataFactory('array'))
    assert numpy.all(nodes['forces'].get_array('forces')[0] == numpy.array([-0.23272115, -0.01115905, 0.03449686]))
    assert numpy.all(nodes['forces'].get_array('forces')[-1] == numpy.array([-0.00300438, 0.00453998, 0.00066599]))


def test_res(parse_result):
    """Check that the results manager can find scalar / low dim results."""
    _, nodes = parse_result()
    output_data = nodes['output_parameters'].get_dict()
    assert output_data['energy'] == -459.8761413
    assert output_data['efermi'] == 2.96801422
    assert 'stress' in output_data
    assert 'dielectric tensor' not in output_data


def test_no_born(parse_result):
    """Make sure no born charges output node exists if lepsilon is not True"""
    _, nodes = parse_result()
    assert 'born_charges' not in nodes


def test_bands(parse_result):
    """Check that bands are parsed and have the right shape."""
    _, nodes = parse_result()
    bands = nodes['bands']
    assert bands.get_bands().shape == (1, 2, 452)


def test_dos(parse_result):
    """Check that dos are parsed."""
    _, nodes = parse_result()
    dos = nodes['dos']
    name, array, units = dos.get_y()[0]
    assert name == 'dos_spin_up'
    assert array.shape == (301,)
    assert units == 'states/eV'


def test_suppress_options(parse_result):
    """Test that suppress options work."""
    _, nodes = parse_result(parser={'parse_dos': False, 'parse_bands': False})
    assert 'dos' not in nodes
    assert 'bands' not in nodes


@pytest.mark.parametrize(['vasprun_path', 'parse_result'], [(2331, False)], indirect=True)
def test_slightly_broken_vasprun(parse_result, recwarn):
    """Test that truncated vasprun (after one ionic step) can be read and warnings are emitted."""
    success, nodes = parse_result()
    assert success
    assert 'kpoints' in nodes
    assert 'structure' in nodes
    assert len(recwarn) >= 1


@pytest.mark.parametrize(['vasprun_path', 'parse_result'], [(2331, True)], indirect=True)
def test_broken_vasprun_exception(parse_result):
    """Test that the parser raises an error with exception_on_bad_xml=True."""
    with pytest.raises(ParsingError):
        _ = parse_result()  # noqa: F841


def test_born_charges(parse_nac, recwarn):
    """Test that the born effective charges get parsed."""
    success, nodes = parse_nac()
    assert success
    output_data = nodes['output_parameters'].get_dict()
    assert numpy.all(output_data['dielectric tensor'] == numpy.array([[5.894056, 0.0, 0.0], [0.0, 5.894056, 0.0], [0.0, 0.0, 6.171821]]))
    assert nodes['born_charges'].get_array('born_charges').shape == (26, 3, 3)
    assert len(recwarn) >= 1
