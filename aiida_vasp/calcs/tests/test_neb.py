"""
Unittests for VaspNEBCalculation
"""
# pylint: disable=unused-import,redefined-outer-name,unused-argument,unused-wildcard-import,wildcard-import, import-outside-toplevel
import contextlib
import os

import pytest

from aiida_vasp.parsers.content_parsers.potcar import MultiPotcarIo
from aiida_vasp.utils.fixtures import *
from aiida_vasp.utils.fixtures.calcs import ONLY_ONE_CALC
from aiida_vasp.utils.aiida_utils import get_data_class, get_data_node


# pylint: disable=protected-access
@ONLY_ONE_CALC
def test_prepare(fresh_aiida_env, aiida_instance, vasp_neb_calc, vasp_neb_inputs, localhost_dir):
    """Check that preparing creates all necessary files."""
    from aiida.common.folders import Folder
    from aiida_vasp.calcs.neb import VaspNEBCalculation

    inputs_dict = {'gga': 'PE', 'gga_compat': False, 'lorbit': 11, 'sigma': 0.5, 'magmom': '30 * 2*0.', 'images': 3}

    inputs = vasp_neb_inputs(parameters=inputs_dict)

    calc = vasp_neb_calc(inputs=inputs)
    temp_folder = Folder(str(localhost_dir))
    calcinfo = calc.prepare_for_submission(temp_folder)
    input_files = temp_folder.get_content_list()

    for file_name in ['INCAR', 'KPOINTS', 'POSCAR', 'POTCAR']:
        assert file_name in input_files

    assert ['01/OUTCAR', '.', 2] in calcinfo.retrieve_list
    assert ['02/CONTCAR', '.', 2] in calcinfo.retrieve_list
    assert 'vasprun.xml' in calcinfo.retrieve_list

    assert calcinfo.codes_info[0].stdout_name == VaspNEBCalculation._VASP_OUTPUT
    assert calcinfo.codes_info[0].join_files is True

    # Test retriving more files
    settings = get_data_node('core.dict', dict={
        'PER_IMAGE_ADDITIONAL_RETRIEVE_LIST': ['IBZKPT'],
    })
    inputs['settings'] = settings
    calc = vasp_neb_calc(inputs=inputs)
    temp_folder.erase(create_empty_folder=True)
    calcinfo = calc.prepare_for_submission(temp_folder)

    assert ['01/IBZKPT', '.', 2] in calcinfo.retrieve_list
