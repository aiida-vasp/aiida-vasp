"""Unittests for ParseVasp."""
# pylint:
# disable=unused-import,redefined-outer-name,unused-argument,unused-wildcard-import,wildcard-import

import pytest

from aiida_vasp.utils.fixtures.testdata import data_path
from aiida_vasp.utils.aiida_utils import get_data_class


def xml_path(folder):
    """Set the path to the XML file."""
    return data_path(folder, 'vasprun.xml')


def xml_truncate(index, original, tmp):
    """Truncate vasprun.xml at the given line number and parse."""
    with open(original, 'r') as xmlfile:
        content = xmlfile.read().splitlines()
        truncated_content = '\n'.join(content[:-index or None])
    with open(tmp, 'w') as xmlfile:
        xmlfile.write(str(truncated_content))


@pytest.fixture(params=[0, 1])
def _parse_me(request, tmpdir):
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
        settings_dict = {'parser_settings': {'add_bands': True, 'output_params': ['fermi_level']}}
        settings_dict.update(extra_settings)
        calc.use_settings(DataFactory('parameter')(dict=settings_dict))
        from aiida_vasp.parsers.vasp import VaspParser
        parser = VaspParser(calc=calc)
        retrieved = DataFactory('folder')()
        fldr = "basic"
        if "folder" in extra_settings:
            fldr = extra_settings["folder"]
        xml_file_path = xml_path(fldr)
        tmp_file_path = str(tmpdir.join('vasprun.xml'))
        #tmp_file_path = os.path.realpath(os.path.join(
        #    __file__, '../../../test_data/tmp/vasprun.xml'))
        xml_truncate(request.param, xml_file_path, tmp_file_path)
        retrieved.add_path(tmp_file_path, '')
        success, nodes = parser.parse_with_retrieved({'retrieved': retrieved})
        nodes = dict(nodes)
        return success, nodes

    return parse


def test_parameters_result(_parse_me):
    """Test that the parameters result node is a KpointsData instance."""

    _, nodes = _parse_me(folder='basic')
    parameters = nodes['output_parameters']
    bands = nodes['output_bands']
    kpoints = nodes['output_kpoints']
    assert isinstance(parameters, get_data_class('parameter'))
    assert isinstance(bands, get_data_class('array.bands'))
    assert isinstance(kpoints, get_data_class('array.kpoints'))
    assert parameters.get_dict()['fermi_level'] == 5.96764939
