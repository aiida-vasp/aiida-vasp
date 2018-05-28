"""Unittests for ParseVasp."""
# pylint:
# disable=unused-import,redefined-outer-name,unused-argument,unused-wildcard-import,wildcard-import

import numpy as np
import pytest

from aiida_vasp.parsers.vasp import VaspParser
from aiida_vasp.utils.fixtures.testdata import data_path
from aiida_vasp.utils.fixtures import *
#from aiida import load_dbenv
#load_dbenv()

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
def parse_result(request, aiida_env, tmpdir):
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
        settings_dict = {'parser_settings': {'parse_potcar_file': False,
                                             'exception_on_bad_xml': False,
                                             'should_parse_OUTCAR': False,
                                             'should_parse_CONTCAR': False}}
        settings_dict.update(extra_settings)
        calc.use_settings(DataFactory('parameter')(dict=settings_dict))
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


def test_kpoints_result(parse_result, vasp_xml):
    """Test that the kpoints result node is a KpointsData instance."""

    from aiida.orm import DataFactory
    settings = {'quantities_to_parse': ['kpoints']}
    quantity = vasp_xml.get_quantity('kpoints', settings = settings)
    data_obj = quantity['kpoints']
    assert isinstance(data_obj, DataFactory('array.kpoints'))
    assert np.all(data_obj.get_kpoints()[0] ==
                  np.array([0.0, 0.0, 0.0]))
    assert np.all(data_obj.get_kpoints()[-1] ==
                  np.array([0.42857143, -0.42857143, 0.42857143]))
