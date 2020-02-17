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


def normalize_contents(file_contents):
    """Remove trailing zeroes after floating point and normalize trailin newline to unix standard."""
    normalized = re.sub(r'(\d*.\d*?)0+(\s)', r'\g<1>0\2', file_contents)  # remove trailing zeroes
    if not re.match(r'\n', file_contents[-1]):  # add trailing newline if necessary
        normalized += '\n'
    return normalized


def assert_contents_equivalent(contents_a, contents_b):
    """Assert equivalence of files with floating point numbers."""
    assert normalize_contents(contents_a) == normalize_contents(contents_b)


@pytest.mark.skip(reason='aiida_wannier90 has not been migrated yet')
@pytest.mark.parametrize(['vasp_structure', 'vasp_kpoints'], [('cif', 'mesh'), ('str', 'list')], indirect=True)
def test_store(vasp2w90_calc_and_ref):
    vasp_calc, _ = vasp2w90_calc_and_ref
    vasp_calc.store_all()
    assert vasp_calc.pk is not None


@pytest.mark.skip(reason='aiida_wannier90 has not been migrated yet')
@ONLY_ONE_CALC
def test_write_incar(vasp2w90_calc_and_ref):
    """Write INCAR reference file and compare to reference."""
    vasp_calc, reference = vasp2w90_calc_and_ref
    with managed_temp_file() as temp_file:
        vasp_calc.write_incar(temp_file)
        with open(temp_file, 'r') as result_incar_fo:
            assert result_incar_fo.read() == reference['incar']


@pytest.mark.skip(reason='aiida_wannier90 has not been migrated yet')
@ONLY_ONE_CALC
def test_write_win(vasp2w90_calc_and_ref):
    """Write wannier90.win input file and compare to reference."""
    vasp_calc, reference = vasp2w90_calc_and_ref
    with managed_temp_file() as temp_file:
        vasp_calc.write_win(temp_file)
        with open(temp_file, 'r') as result_incar_fo:
            assert result_incar_fo.read() == reference['win']


@contextlib.contextmanager
def managed_temp_file():
    """Create a temp file for a with context, delete after use."""
    import tempfile
    _, temp_file = tempfile.mkstemp()
    try:
        yield temp_file
    finally:
        os.remove(temp_file)
