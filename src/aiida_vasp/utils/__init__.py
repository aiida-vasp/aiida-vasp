"""
Utility functions for VASP calculations.
"""

import numpy as np


def get_maximum_force(forces):
    norm = np.linalg.norm(forces, axis=1)
    return np.amax(np.abs(norm))
