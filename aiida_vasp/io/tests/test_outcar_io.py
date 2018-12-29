"""Test the vasprun.xml io interface"""
# pylint: disable=unused-import,redefined-outer-name,unused-argument,unused-wildcard-import,wildcard-import

import pytest
import numpy as np

from aiida_vasp.utils.fixtures import *
from aiida_vasp.utils.fixtures.testdata import data_path
from aiida_vasp.utils.aiida_utils import get_data_class
from aiida_vasp.parsers.node_composer import NodeComposer


@pytest.mark.parametrize('outcar_parser', ['disp_details'], indirect=True)
def test_parameter_results(outcar_parser):
    """
    Test that the parameter node is a ParametersData instance.

    Should contain the symmetries and the elastic moduli.

    """

    outcar_parser.settings.nodes.update({
        'parameters': {
            'type': 'parameter',
            'quantities': ['symmetries', 'elastic_moduli'],
            'link_name': 'my_custom_node'
        }
    })

    composer = NodeComposer(file_parsers=[outcar_parser])
    data_obj = composer.compose('parameter', quantities=['symmetries', 'elastic_moduli'])
    ref_class = get_data_class('parameter')
    assert isinstance(data_obj, ref_class)
    data_dict = data_obj.get_dict()
    # test symmetries
    test = {
        'symmetrized_cell_type': {
            'static': [
                'face centered cubic supercell.', 'body centered tetragonal supercell.', 'body centered tetragonal supercell.',
                'body centered tetragonal supercell.', 'body centered tetragonal supercell.', 'body centered tetragonal supercell.',
                'body centered tetragonal supercell.', 'base centered monoclinic supercell.', 'base centered monoclinic supercell.',
                'base centered monoclinic supercell.', 'base centered monoclinic supercell.', 'base centered monoclinic supercell.',
                'base centered monoclinic supercell.', 'face centered cubic supercell.', 'face centered cubic supercell.',
                'face centered cubic supercell.'
            ],
            'dynamic': [
                'face centered cubic supercell.', 'body centered tetragonal supercell.', 'body centered tetragonal supercell.',
                'body centered tetragonal supercell.', 'body centered tetragonal supercell.', 'body centered tetragonal supercell.',
                'body centered tetragonal supercell.', 'base centered monoclinic supercell.', 'base centered monoclinic supercell.',
                'base centered monoclinic supercell.', 'base centered monoclinic supercell.', 'base centered monoclinic supercell.',
                'base centered monoclinic supercell.', 'face centered cubic supercell.', 'face centered cubic supercell.',
                'face centered cubic supercell.'
            ]
        },
        'num_point_group_operations': {
            'static': [24, 8, 8, 8, 8, 8, 8, 2, 2, 2, 2, 2, 2, 4, 4, 24],
            'dynamic': [24, 8, 8, 8, 8, 8, 8, 2, 2, 2, 2, 2, 2, 4, 4, 24]
        },
        'point_group': {
            'static': [
                'T_d', 'D_2d.', 'D_2d.', 'D_2d.', 'D_2d.', 'D_2d.', 'D_2d.', 'C_2', 'C_2', 'C_2', 'C_2', 'C_2', 'C_2', 'C_2v.', 'C_2v.',
                'T_d'
            ],
            'dynamic': [
                'T_d', 'D_2d.', 'D_2d.', 'D_2d.', 'D_2d.', 'D_2d.', 'D_2d.', 'C_2', 'C_2', 'C_2', 'C_2', 'C_2', 'C_2', 'C_2v.', 'C_2v.',
                'T_d'
            ]
        },
        'space_group': {
            'static': [
                'O_h', 'D_4h.', 'D_4h.', 'D_4h.', 'D_4h.', 'D_4h.', 'D_4h.', 'C_2h.', 'C_2h.', 'C_2h.', 'C_2h.', 'C_2h.', 'C_2h.', 'D_2h.',
                'D_2h.', 'O_h'
            ],
            'dynamic': [
                'O_h', 'D_4h.', 'D_4h.', 'D_4h.', 'D_4h.', 'D_4h.', 'D_4h.', 'C_2h.', 'C_2h.', 'C_2h.', 'C_2h.', 'C_2h.', 'C_2h.', 'D_2h.',
                'D_2h.', 'O_h'
            ]
        },
        'original_cell_type': {
            'static': [
                'primitive cell', 'primitive cell', 'primitive cell', 'primitive cell', 'primitive cell', 'primitive cell',
                'primitive cell', 'primitive cell', 'primitive cell', 'primitive cell', 'primitive cell', 'primitive cell',
                'primitive cell', 'primitive cell', 'primitive cell', 'primitive cell'
            ],
            'dynamic': [
                'primitive cell', 'primitive cell', 'primitive cell', 'primitive cell', 'primitive cell', 'primitive cell',
                'primitive cell', 'primitive cell', 'primitive cell', 'primitive cell', 'primitive cell', 'primitive cell',
                'primitive cell', 'primitive cell', 'primitive cell', 'primitive cell'
            ]
        },
        'num_space_group_operations': {
            'static': [48, 16, 16, 16, 16, 16, 16, 4, 4, 4, 4, 4, 4, 8, 8, 48],
            'dynamic': [48, 16, 16, 16, 16, 16, 16, 4, 4, 4, 4, 4, 4, 8, 8, 48]
        }
    }
    # then elastic moduli
    test = np.array([1674.5786, 704.739, 704.739, -0.0, 0.0, 0.0])
    np.testing.assert_allclose(data_dict['elastic_moduli']['symmetrized'][0], test)
    test = np.array([0.0, 0.0, 0.0, -0.0, -0.0, 1122.6622])
    np.testing.assert_allclose(data_dict['elastic_moduli']['symmetrized'][5], test)
    test = np.array([705.0238, 1674.8491, 705.0238, -0.0, -0.0, 0.0])
    np.testing.assert_allclose(data_dict['elastic_moduli']['non_symmetrized'][1], test)
    test = np.array([-0.0078, -0.0495, 0.0147, 0.0, 1123.0829, -0.0])
    np.testing.assert_allclose(data_dict['elastic_moduli']['non_symmetrized'][4], test)
    test = np.array([704.739, 704.739, 1674.5786, -0.0, -0.0, 0.0])
    np.testing.assert_allclose(data_dict['elastic_moduli']['total'][2], test)
    test = np.array([-0.0, -0.0, -0.0, 775.8054, 0.0, -0.0])
    np.testing.assert_allclose(data_dict['elastic_moduli']['total'][3], test)
