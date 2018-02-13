"""Parser unit (?) test"""
import os

import pytest

from aiida_vasp.io.win import WinParser


@pytest.fixture(params=['0.win', '1.win'])
def win_result(request, sample):
    filename = request.param
    win_parser = WinParser(sample(os.path.join('win_files', filename)))
    return win_parser.result


def test_keys(win_result):  # pylint: disable=redefined-outer-name
    assert set(win_result.keys()) == set([
        'projections',
        'spinors',
        'unit_cell_cart',
        'atoms_cart',
        'mp_grid',
        'kpoints',
        'num_wann',
    ])


def test_no_empty_key(win_result):  # pylint: disable=redefined-outer-name
    assert '' not in win_result.keys()


def test_projections(win_result):  # pylint: disable=redefined-outer-name
    assert all(isinstance(line, str) for line in win_result['projections'])
