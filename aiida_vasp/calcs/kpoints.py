"""
Utilities for preparing VASP - KPOINTS files

provides write_kpoints(argdict), modified from ase.calculators.vasp to accomodate Line Mode
"""
import numpy as np


def write_kpoints(argdict):
    """Writes the KPOINTS file. Modified from ase version to allow using Line Mode for band structure calculation."""
    params = argdict
    kpoints = open('KPOINTS', 'w')
    kpoints.write('KPOINTS created by Atomic Simulation Environment\n')
    params['kpts'] = params['kpoints']
    shape = np.array(params['kpts']).shape
    if len(shape) == 1:
        kpoints.write('0\n')
        if params.get('gamma'):
            kpoints.write('Gamma\n')
        else:
            kpoints.write('Monkhorst-Pack\n')
        for kpt in params['kpts']:
            kpoints.write('%i ' % kpt)
        kpoints.write('\n0 0 0\n')
    elif len(shape) == 2:
        kpoints.write('%i \n' % (params.get('intersections') or len(params['kpts'])))
        kpoints.write('Line-mode \n')
        if params['reciprocal']:
            kpoints.write('Reciprocal\n')
        else:
            kpoints.write('Cartesian\n')
        for i in range(len(params['kpts'])):
            for kpt in params['kpts'][i]:
                kpoints.write('%f ' % kpt)
            if shape[1] == 4:
                kpoints.write('\n')
            elif shape[1] == 3:
                kpoints.write('1.0 \n')
    kpoints.close()
