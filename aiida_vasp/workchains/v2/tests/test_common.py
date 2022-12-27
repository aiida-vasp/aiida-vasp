"""
Tests for the common module
"""
from aiida.orm import Dict, List

from aiida_vasp.workchains.v2.common import aiida_to_python, plain_python_args, site_magnetization_to_magmom

#pylint:disable=unused-argument


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
                'y': {
                    'site_moment': {},
                    'total_magnetization': {}
                },
                'z': {
                    'site_moment': {},
                    'total_magnetization': {}
                },
            },
            'full_cell': [15.9999942],
        }
    }

    assert site_magnetization_to_magmom(output) == [0.505, 2.866]
    assert site_magnetization_to_magmom(Dict(dict=output)) == [0.505, 2.866]
