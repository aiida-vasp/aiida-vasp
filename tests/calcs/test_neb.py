"""
Unittests for VaspNEBCalculation
"""

import pytest
from aiida import orm

from aiida_vasp.calcs.neb import VaspNEBCalculation


@pytest.mark.parametrize(['vasp_structure', 'vasp_kpoints'], [('cif', 'mesh')], indirect=True)
def test_prepare(
    aiida_profile,
    vasp_neb_calc,
    vasp_neb_inputs,
    sandbox_folder,
):
    """Check that preparing creates all necessary files."""
    inputs_dict = {
        'gga': 'PE',
        'gga_compat': False,
        'lorbit': 11,
        'sigma': 0.5,
        'magmom': '30 * 2*0.',
        'images': 3,
    }

    inputs = vasp_neb_inputs(parameters=inputs_dict)

    calc = vasp_neb_calc(inputs=inputs)
    temp_folder = sandbox_folder
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
    settings = orm.Dict(
        dict={
            'PER_IMAGE_ADDITIONAL_RETRIEVE_LIST': ['IBZKPT'],
        },
    )
    inputs['settings'] = settings
    calc = vasp_neb_calc(inputs=inputs)
    temp_folder.erase(create_empty_folder=True)
    calcinfo = calc.prepare_for_submission(temp_folder)

    assert ['01/IBZKPT', '.', 2] in calcinfo.retrieve_list
