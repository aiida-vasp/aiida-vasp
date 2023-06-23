"""Test the magnetization parser."""
# pylint: disable=unused-import,redefined-outer-name,unused-argument,unused-wildcard-import,wildcard-import

import pytest
import numpy as np

from aiida_vasp.utils.fixtures import *
from aiida_vasp.utils.aiida_utils import get_data_class
from aiida_vasp.parsers.node_composer import NodeComposer, get_node_composer_inputs_from_file_parser


@pytest.mark.parametrize('outcar_parser', ['magnetization'], indirect=True)
def test_magnetization_parser(fresh_aiida_env, outcar_parser):
    """
    Test that the magnetization node is a ParametersData instance.

    Should contain the site magnetization and the full cell magnetization

    """

    outcar_parser._settings._output_nodes_dict.update(  # pylint: disable=protected-access
        {'add_site_magnetization': {
            'link_name': 'site_magnetization',
            'type': 'dict',
            'quantities': ['site_magnetization']
        }})

    inputs = get_node_composer_inputs_from_file_parser(outcar_parser, quantity_keys=['site_magnetization'])
    data_obj = NodeComposer.compose('dict', inputs)
    ref_class = get_data_class('dict')
    assert isinstance(data_obj, ref_class)
    data_dict = data_obj.get_dict()
    # test symmetries
    test = {
        'sphere': {
            'x': {
                'site_moment': {
                    1: {
                        's': -0.014,
                        'p': -0.051,
                        'd': 1.687,
                        'tot': 1.621
                    },
                    2: {
                        's': -0.015,
                        'p': -0.052,
                        'd': 1.686,
                        'tot': 1.619
                    },
                    3: {
                        's': -0.014,
                        'p': -0.053,
                        'd': 1.708,
                        'tot': 1.64
                    },
                    4: {
                        's': -0.014,
                        'p': -0.053,
                        'd': 1.708,
                        'tot': 1.64
                    }
                },
                'total_magnetization': {
                    's': -0.057,
                    'p': -0.21,
                    'd': 6.788,
                    'tot': 6.521
                }
            },
            'y': {
                'site_moment': {},
                'total_magnetization': {}
            },
            'z': {
                'site_moment': {},
                'total_magnetization': {}
            }
        },
        'full_cell': np.asarray([6.4424922])
    }
    assert set(data_dict['site_magnetization']) == set(test)


@pytest.mark.parametrize('outcar_parser', ['magnetization_single'], indirect=True)
def test_magnetization_single_parser(fresh_aiida_env, outcar_parser):  # pylint: disable=invalid-name
    """
    Test that the magnetization node is a ParametersData instance.

    Should contain the site magnetization and the full cell magnetization when
    there is a single atom in the unit cell

    """

    outcar_parser._settings._output_nodes_dict.update(  # pylint: disable=protected-access
        {'add_site_magnetization': {
            'link_name': 'site_magnetization',
            'type': 'dict',
            'quantities': ['site_magnetization']
        }})

    inputs = get_node_composer_inputs_from_file_parser(outcar_parser, quantity_keys=['site_magnetization'])
    data_obj = NodeComposer.compose('dict', inputs)
    ref_class = get_data_class('dict')
    assert isinstance(data_obj, ref_class)
    data_dict = data_obj.get_dict()
    # test symmetries
    test = {
        'sphere': {
            'x': {
                'site_moment': {
                    1: {
                        's': -0.012,
                        'p': -0.043,
                        'd': 2.49,
                        'tot': 2.434
                    },
                },
                'total_magnetization': {
                    's': -0.012,
                    'p': -0.043,
                    'd': 2.49,
                    'tot': 2.434
                }
            },
            'y': {
                'site_moment': {},
                'total_magnetization': {}
            },
            'z': {
                'site_moment': {},
                'total_magnetization': {}
            }
        },
        'full_cell': np.asarray([2.4077611])
    }
    assert set(data_dict['site_magnetization']) == set(test)
