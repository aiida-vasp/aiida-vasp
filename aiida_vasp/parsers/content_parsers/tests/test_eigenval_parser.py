"""Test the EIGENVAL parser."""
# pylint: disable=unused-import,redefined-outer-name,unused-argument,unused-wildcard-import,wildcard-import

import numpy as np
import pytest

from aiida_vasp.utils.fixtures import *
from aiida_vasp.utils.fixtures.testdata import data_path


@pytest.mark.parametrize(['eigenval_parser'], [('eigenval',)], indirect=True)
def test_parse_eigenval(eigenval_parser):
    """Parse a reference EIGENVAL with the EigenvalParser and compare the result to a reference."""
    result = eigenval_parser.get_quantity('eigenval-eigenvalues')
    eigenvalues = np.array([[[-1.439825, 2.964373, 2.964373, 2.964373, 7.254542, 7.254542, 7.254542, 11.451811, 11.670398, 11.670398]]])
    metadata = {
        0: [4, 4, 1, 1],
        1: [16.48482, 4.04e-10, 4.04e-10, 4.04e-10, 1e-16],
        2: 0.0001,
        'n_ions': 4,
        'n_atoms': 4,
        'p00': 1,
        'nspin': 1,
        'cartesian': True,
        'name': 'unknown system',
        'some_num': 12,
        'n_bands': 10,
        'n_kp': 1
    }
    assert result['metadata'] == metadata
    assert np.allclose(result['eigenvalues'], eigenvalues)
