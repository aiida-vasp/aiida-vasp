"""
Tests for the common module
"""

from aiida.orm import Dict, List, StructureData

from aiida_vasp.common import aiida_to_python, plain_python_args, site_magnetization_to_magmom
from aiida_vasp.common.transform import magnetic_structure_decorate, magnetic_structure_dedecorate

# pylint:disable=unused-argument


def test_conversion(aiida_profile):
    """Test type dispatch conversion"""
    node = Dict(dict={'x': 1})
    assert aiida_to_python(node) == {'x': 1}

    @plain_python_args
    def func(node, *args, **kwargs):
        return node

    assert func(node) == {'x': 1}

    node = List(list=[1, 2, 3])
    assert aiida_to_python(node) == [1, 2, 3]


def test_magmom_from_site(aiida_profile):
    """Test exacting magmom"""
    output = {
        'site_magnetization': {
            'sphere': {
                'x': {
                    'site_moment': {
                        '1': {
                            'd': 0.472,
                            'f': 0.0,
                            'p': 0.021,
                            's': 0.011,
                            'tot': 0.505,
                        },
                        '2': {
                            'd': 2.851,
                            'f': 0.0,
                            'p': 0.008,
                            's': 0.007,
                            'tot': 2.866,
                        },
                    },
                    'total_magnetization': {
                        'd': 13.307,
                        'f': -0.012,
                        'p': 2.148,
                        's': 0.247,
                        'tot': 15.69,
                    },
                },
                'y': {'site_moment': {}, 'total_magnetization': {}},
                'z': {'site_moment': {}, 'total_magnetization': {}},
            },
            'full_cell': [15.9999942],
        }
    }

    assert site_magnetization_to_magmom(output) == [0.505, 2.866]
    assert site_magnetization_to_magmom(Dict(dict=output)) == [0.505, 2.866]


def test_magmom_from_site_noncollinear(aiida_profile):
    """Test extracting non-collinear magmom (3-vectors)"""
    output = {
        'site_magnetization': {
            'sphere': {
                'x': {
                    'site_moment': {
                        '1': {'tot': 0.5},
                        '2': {'tot': 1.5},
                    },
                },
                'y': {
                    'site_moment': {
                        '1': {'tot': 0.0},
                        '2': {'tot': 0.0},
                    },
                },
                'z': {
                    'site_moment': {
                        '1': {'tot': -0.5},
                        '2': {'tot': 1.5},
                    },
                },
            },
        }
    }
    result = site_magnetization_to_magmom(output)
    assert result == [(0.5, 0.0, -0.5), (1.5, 0.0, 1.5)]
    assert site_magnetization_to_magmom(Dict(dict=output)) == [(0.5, 0.0, -0.5), (1.5, 0.0, 1.5)]


def _make_feo_structure():
    """Tiny FeO-like structure used in the magmom decoration tests."""
    s = StructureData()
    s.set_cell([[3.0, 0.0, 0.0], [0.0, 3.0, 0.0], [0.0, 0.0, 3.0]])
    s.append_atom(position=[0.0, 0.0, 0.0], symbols='Fe')
    s.append_atom(position=[1.5, 1.5, 1.5], symbols='Fe')
    s.append_atom(position=[0.0, 1.5, 1.5], symbols='O')
    s.append_atom(position=[1.5, 0.0, 0.0], symbols='O')
    return s


def test_magnetic_structure_decorate_collinear(aiida_profile):
    """Anti-ferromagnetic FeO collinear case: Fe atoms should split into Fe1/Fe2."""
    structure = _make_feo_structure()
    magmom = List(list=[2.0, -2.0, 0.0, 0.0])
    result = magnetic_structure_decorate(structure, magmom)
    assert [site.kind_name for site in result['structure'].sites] == ['Fe1', 'Fe2', 'O', 'O']
    mapping = result['mapping'].get_dict()
    assert mapping == {'Fe1': 2.0, 'Fe2': -2.0, 'O': 0.0}

    dedecorate = magnetic_structure_dedecorate(result['structure'], result['mapping'])
    assert dedecorate['magmom'].get_list() == [2.0, -2.0, 0.0, 0.0]


def test_magnetic_structure_decorate_noncollinear(aiida_profile):
    """Non-collinear magmom should split species per unique 3-vector."""
    structure = _make_feo_structure()
    magmom = List(
        list=[
            [2.0, 0.0, 0.0],
            [-2.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
        ]
    )
    result = magnetic_structure_decorate(structure, magmom)
    assert [site.kind_name for site in result['structure'].sites] == ['Fe1', 'Fe2', 'O', 'O']
    mapping = result['mapping'].get_dict()
    assert mapping == {'Fe1': [2.0, 0.0, 0.0], 'Fe2': [-2.0, 0.0, 0.0], 'O': [0.0, 0.0, 0.0]}

    dedecorate = magnetic_structure_dedecorate(result['structure'], result['mapping'])
    roundtrip = dedecorate['magmom'].get_list()
    assert roundtrip == [[2.0, 0.0, 0.0], [-2.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]]
