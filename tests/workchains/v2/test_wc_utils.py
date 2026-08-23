import numpy as np
import pytest
from aiida.common.extendeddicts import AttributeDict

from aiida_vasp.utils.workchains import site_magnetization_to_magmom
from aiida_vasp.workchains.v2.bands import _magmom_list_to_incar as bands_magmom_list_to_incar
from aiida_vasp.workchains.v2.relax import VaspRelaxWorkChain, _magmom_list_to_incar, get_maximum_force


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


def _site_magnetization_output(site_moments):
    """Build a ``site_magnetization`` output dict from per-axis site moments."""
    sphere = {}
    for axis in 'xyz':
        sphere[axis] = {
            'site_moment': {str(i + 1): {'tot': value} for i, value in enumerate(site_moments.get(axis, []))}
        }
    return {'site_magnetization': {'sphere': sphere}}


def test_site_magnetization_to_magmom_collinear():
    """Single populated axis returns scalar magmoms per site."""
    output = _site_magnetization_output({'x': [0.505, 2.866]})
    assert site_magnetization_to_magmom(output) == [0.505, 2.866]


def test_site_magnetization_to_magmom_noncollinear():
    """Three populated axes return one 3-tuple per site."""
    output = _site_magnetization_output({'x': [0.5, 1.5], 'y': [0.0, 0.0], 'z': [-0.5, 1.5]})
    assert site_magnetization_to_magmom(output) == [(0.5, 0.0, -0.5), (1.5, 0.0, 1.5)]


def test_site_magnetization_to_magmom_two_axes():
    """A missing middle axis is filled with zeros."""
    output = _site_magnetization_output({'x': [0.5, 1.5], 'z': [-0.5, 1.5]})
    assert site_magnetization_to_magmom(output) == [(0.5, 0.0, -0.5), (1.5, 0.0, 1.5)]


def test_site_magnetization_to_magmom_no_axes_raises():
    """No populated axis raises a clear error."""
    with pytest.raises(ValueError):
        site_magnetization_to_magmom({'site_magnetization': {'sphere': {'x': {'site_moment': {}}}}})


def test_relax_magmom_list_to_incar():
    """The relax workchain flattens mixed scalar/vector magmoms into an INCAR string."""
    assert _magmom_list_to_incar([2.0, -2.0]) == '2.0 -2.0'
    assert _magmom_list_to_incar([[2.0, 0.0, 0.0], [-2.0, 0.0, 0.0]]) == '2.0 0.0 0.0 -2.0 0.0 0.0'
    assert _magmom_list_to_incar([2.0, (1.0, 0.0, 0.0)]) == '2.0 1.0 0.0 0.0'


def test_bands_magmom_list_to_incar():
    """The bands workchain flattens mixed scalar/vector magmoms into an INCAR string."""
    assert bands_magmom_list_to_incar([2.0, -2.0]) == '2.0 -2.0'
    assert bands_magmom_list_to_incar([(2.0, 0.0, 0.0), (-2.0, 0.0, 0.0)]) == '2.0 0.0 0.0 -2.0 0.0 0.0'
