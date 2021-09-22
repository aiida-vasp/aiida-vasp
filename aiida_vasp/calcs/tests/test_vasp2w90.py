"""Unittests for Vasp2w90Calculation."""
# pylint: disable=unused-import,redefined-outer-name,unused-argument,unused-wildcard-import,wildcard-import, import-outside-toplevel
from __future__ import absolute_import
from __future__ import print_function
import contextlib
import os
import re

import pytest

from aiida_vasp.utils.fixtures import *
from aiida_vasp.utils.fixtures.calcs import ONLY_ONE_CALC


def normalize_contents(object_contents):
    """Remove trailing zeroes after floating point and normalize trailin newline to unix standard."""
    normalized = re.sub(r'(\d*.\d*?)0+(\s)', r'\g<1>0\2', object_contents)  # remove trailing zeroes
    if not re.match(r'\n', object_contents[-1]):  # add trailing newline if necessary
        normalized += '\n'
    return normalized


def assert_contents_equivalent(contents_a, contents_b):
    """Assert equivalence of objects with floating point numbers."""
    assert normalize_contents(contents_a) == normalize_contents(contents_b)


@ONLY_ONE_CALC
def test_wannier_parameters(vasp2w90_calc_and_ref):
    """Test that Wannier90 parameters and projections can be specified."""
    vasp_calc, _ = vasp2w90_calc_and_ref
    projections = vasp_calc.inputs.wannier_projections.get_list()
    parameters = vasp_calc.inputs.wannier_parameters.get_dict()
    assert projections == ['Ga : s; px; py; pz', 'As : px; py; pz']
    assert parameters == {'spinors': True, 'num_iter': 0, 'num_wann': 14, 'num_bands': 24, 'dis_num_iter': 1000}


@ONLY_ONE_CALC
def test_prepare_for_submission(vasp2w90_calc_and_ref, tmp_path):
    """Test that the lwannier90 flag is written at the prepare for submission step."""
    from aiida.common.folders import Folder
    vasp_calc, reference = vasp2w90_calc_and_ref
    folder = Folder(os.fspath(tmp_path))
    vasp_calc.prepare_for_submission(folder)
    with managed_temp_object() as temp_object:
        vasp_calc.write_incar(temp_object)
        with open(temp_object, 'r') as result_incar_fo:
            assert result_incar_fo.readlines() == reference['incar']


@ONLY_ONE_CALC
def test_write_win(vasp2w90_calc_and_ref):
    """Write wannier90.win input object and compare to reference."""
    vasp_calc, reference = vasp2w90_calc_and_ref
    with managed_temp_object() as temp_object:
        vasp_calc.write_win(temp_object)
        with open(temp_object, 'r') as result_incar_fo:
            assert result_incar_fo.read() == reference['win']


@contextlib.contextmanager
def managed_temp_object():
    """Create a temp file object for a with context, delete after use."""
    import tempfile
    _, temp_object = tempfile.mkstemp()
    try:
        yield temp_object
    finally:
        os.remove(temp_object)
