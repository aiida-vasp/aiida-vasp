from aiida_vasp.utils.io.win import WinParser

def test_keys(sample):
    wp = WinParser(sample('wannier90.win'))
    result = wp.result
    set(result.keys()) == set([
        'projections', 'spinors', 'unit_cell_cart', 'atoms_cart', 'mp_grid',
        'kpoints'
    ])

def test_no_empty_key(sample):
    wp = WinParser(sample('wannier90.win'))
    result = wp.result
    assert '' not in result.keys()
