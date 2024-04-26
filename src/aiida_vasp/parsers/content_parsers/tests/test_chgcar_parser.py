"""Test the CHGCAR parser."""
# pylint: disable=unused-import,redefined-outer-name,unused-argument,unused-wildcard-import,wildcard-import

import numpy as np
import pytest
from aiida_vasp.utils.fixtures import *

compare_charge_density = np.array(
    [
        [
            [0.09329446, 0.18658892, 0.27988338],
            [0.37317784, 0.4664723, 0.55976676],
            [0.65306122, 0.74635569, 0.83965015],
            [0.93294461, 1.02623907, 1.11953353],
        ],
        [
            [1.21282799, 1.30612245, 1.39941691],
            [1.49271137, 1.58600583, 1.67930029],
            [1.77259475, 1.86588921, 1.95918367],
            [2.05247813, 2.14577259, 2.23906706],
        ],
        [
            [2.33236152, 2.42565598, 2.51895044],
            [2.6122449, 2.70553936, 2.79883382],
            [2.89212828, 2.98542274, 3.0787172],
            [3.17201166, 3.26530612, 3.35860058],
        ],
        [
            [3.45189504, 3.5451895, 3.63848397],
            [3.73177843, 3.82507289, 3.91836735],
            [4.01166181, 4.10495627, 4.19825073],
            [4.29154519, 4.38483965, 4.47813411],
        ],
        [
            [4.57142857, 4.66472303, 4.75801749],
            [4.85131195, 4.94460641, 5.03790087],
            [5.13119534, 5.2244898, 5.31778426],
            [5.41107872, 5.50437318, 5.59766764],
        ],
    ]
)


@pytest.mark.parametrize(['chgcar_parser'], [('chgcar',)], indirect=True)
def test_parse_charge(chgcar_parser):
    """Load a reference CHGCAR parser.

    We check that the it parses and provides the correct content for the charge density.

    """
    result = chgcar_parser.get_quantity('charge_density')
    assert result is not None
    assert 'charge_density' in result
    assert np.allclose(result['charge_density'], compare_charge_density)


@pytest.mark.parametrize(
    ['chgcar_parser'],
    [
        (
            [
                'chgcar',
                'CHGCAR.spin',
                {'quantities_to_parse': ['charge_density', 'magnetization_density']},
            ],
        )
    ],
    indirect=True,
)
def test_parse_magnetization(chgcar_parser):
    """Load a reference CHGCAR parser.

    We check that the it parses and provides the correct content for the magnetization density.

    """
    result = chgcar_parser.get_quantity('charge_density')
    assert result is not None
    assert 'charge_density' in result
    assert np.allclose(result['charge_density'], compare_charge_density)
    # Magnetization density if turned off by default, enable it
    result = chgcar_parser.get_quantity('magnetization_density')
    assert result is not None
    assert 'magnetization_density' in result
    assert np.allclose(result['magnetization_density'], compare_charge_density)


@pytest.mark.parametrize(
    ['chgcar_parser'],
    [
        (
            [
                'chgcar',
                'CHGCAR.ncl',
                {'quantities_to_parse': ['charge_density', 'magnetization_density']},
            ],
        )
    ],
    indirect=True,
)
def test_parse_magnetization_ncl(chgcar_parser):
    """Load a reference CHGCAR parser.

    We check that the it parses and provides the correct content for the magnetization density.

    """
    result = chgcar_parser.get_quantity('charge_density')
    assert result is not None
    assert 'charge_density' in result
    assert np.allclose(result['charge_density'], compare_charge_density)
    result = chgcar_parser.get_quantity('magnetization_density')
    assert result is not None
    assert 'magnetization_density' in result
    assert set(['x', 'y', 'z']) == set(result['magnetization_density'].keys())
    for item in result['magnetization_density'].values():
        assert np.allclose(item, compare_charge_density)
