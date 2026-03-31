import numpy as np
import pytest
from aiida.common.extendeddicts import AttributeDict

from aiida_vasp.workchains.v2.relax import VaspRelaxWorkChain, get_maximum_force


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


def test_check_shape_convergence_uses_angle_threshold():
    """A pure angle change above the angle threshold should not be marked converged."""
    workchain = object.__new__(VaspRelaxWorkChain)
    object.__setattr__(
        workchain,
        '_context',
        AttributeDict(
            {
                'relax_settings': AttributeDict(
                    {
                        'convergence_shape_angles': 0.1,
                        'convergence_shape_lengths': 0.01,
                    }
                )
            }
        ),
    )
    workchain.report = lambda *args, **kwargs: None
    workchain.is_verbose = lambda: False

    delta = AttributeDict({'cell_lengths': np.array([0.0, 0.0, 0.0]), 'cell_angles': np.array([0.2, 0.0, 0.0])})

    assert workchain.check_shape_convergence(delta) is False
