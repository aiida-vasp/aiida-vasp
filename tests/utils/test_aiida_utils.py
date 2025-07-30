"""Test aiida_utils functionss."""

import warnings

import pytest

from aiida_vasp.utils.aiida_utils import convert_dict_case, get_current_user

try:
    from aiida_vasp.utils.pmg import PymatgenAdapator, get_incar, get_kpoints, get_outcar, get_vasprun
except ImportError:
    PMG_INSTALLED = False
else:
    PMG_INSTALLED = True


def test_get_current_user(aiida_profile_clean):
    """Assert that get_current_user returns a user."""
    user = get_current_user()
    assert user.pk
    assert user.first_name == ''
    assert user.last_name == ''
    assert user.email


@pytest.mark.skipif(not PMG_INSTALLED, reason='pymatgen not installed')
@pytest.mark.parametrize(['vasp_structure', 'vasp_kpoints'], [('str', 'mesh')], indirect=True)
def test_pmg_adaptor(aiida_profile_clean, tmp_path, run_vasp_process):
    """
    Test export vasp calculation
    """

    _, node = run_vasp_process()

    adapt = PymatgenAdapator(node)
    assert adapt.vasprun

    assert 'vasprun' in adapt.pmg_objects

    with PymatgenAdapator(node) as adapt:
        adapt.vasprun

    assert 'vasprun_dict' in adapt.cache
    assert 'pmg_cache' in node.base.extras.all

    assert get_incar(node)
    assert get_kpoints(node)
    assert get_vasprun(node)
    assert get_outcar(node)


def test_convert_dict_case_lowercase():
    input_dict = {'KeyOne': 1, 'KEYTWO': 2}
    expected_output = {'keyone': 1, 'keytwo': 2}
    assert convert_dict_case(input_dict) == expected_output


def test_convert_dict_case_uppercase():
    input_dict = {'keyone': 1, 'keytwo': 2}
    expected_output = {'KEYONE': 1, 'KEYTWO': 2}
    assert convert_dict_case(input_dict, lower=False) == expected_output


def test_convert_dict_case_nested():
    input_dict = {'KeyOne': {'SubKeyOne': 1, 'SUBKEYTWO': 2}, 'KEYTWO': 2}
    expected_output = {'keyone': {'subkeyone': 1, 'subkeytwo': 2}, 'keytwo': 2}
    assert convert_dict_case(input_dict) == expected_output


def test_convert_dict_case_warning(caplog):
    input_dict = {'KeyOne': 1, 'KEYTWO': 2}
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter('always')
        convert_dict_case(input_dict, warn=True)
        assert len(w) == 2
        assert "Key 'KeyOne' converted to 'keyone' - please use lower case keys" in str(w[0].message)
        assert "Key 'KEYTWO' converted to 'keytwo' - please use lower case keys" in str(w[1].message)


def test_convert_dict_case_raise():
    input_dict = {'KeyOne': 1, 'KEYTWO': 2}
    with pytest.raises(ValueError, match="Key 'KeyOne' converted to 'keyone' - please use lower case keys"):
        convert_dict_case(input_dict, raise_convert=True)


def test_convert_dict_case_no_conversion():
    input_dict = {'keyone': 1, 'keytwo': 2}
    assert convert_dict_case(input_dict) == input_dict
