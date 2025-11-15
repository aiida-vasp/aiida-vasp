import numpy as np
import pytest

from aiida_vasp.workchains.v2.relax import get_maximum_force


def test_get_maximum_forces():
    """Test the get maximum forces function"""

    forces = np.array([[0.1, 0.2, 0.3], [-0.1, -0.2, -0.3], [1.0, 1.0, 1.0], [-0.4, -0.5, -0.6], [0.0, 0.0, 0.0]])

    assert pytest.approx(get_maximum_force(forces)) == np.sqrt(3)
    # with mask
    mask = np.array(
        [
            [True, True, True],
            [True, True, True],
            [False, True, True],
            [True, True, True],
            [True, True, True],
        ]
    )
    assert get_maximum_force(forces, mask) == pytest.approx(np.sqrt(0.4**2 + 0.5**2 + 0.6**2), rel=1e-6, abs=1e-8)
