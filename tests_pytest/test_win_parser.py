"""Parser unit (?) test"""
from aiida_vasp.utils.io.win import WinParser


def test_keys(sample):
    win_parser = WinParser(sample('wannier90.win'))
    result = win_parser.result
    assert set(result.keys()) == set([
        'projections', 'spinors', 'unit_cell_cart', 'atoms_cart', 'mp_grid',
        'kpoints'
    ])


def test_no_empty_key(sample):
    win_parser = WinParser(sample('wannier90.win'))
    result = win_parser.result
    assert '' not in result.keys()
