import numpy as np

__doc__ = 'provides write_kpoints(argdict), modified from ase.calculators.vasp to accomodate Line Mode'


def write_kpoints(argdict):
    """Writes the KPOINTS file. Modified from ase version to allow using Line Mode for band structure calculation."""
    p = argdict
    kpoints = open('KPOINTS', 'w')
    kpoints.write('KPOINTS created by Atomic Simulation Environment\n')
    p['kpts'] = p['kpoints']
    shape = np.array(p['kpts']).shape
    if len(shape) == 1:
        kpoints.write('0\n')
        if p.get('gamma'):
            kpoints.write('Gamma\n')
        else:
            kpoints.write('Monkhorst-Pack\n')
        [kpoints.write('%i ' % kpt) for kpt in p['kpts']]
        kpoints.write('\n0 0 0\n')
    elif len(shape) == 2:
        kpoints.write('%i \n' % (p.get('intersections') or len(p['kpts'])))
        kpoints.write('Line-mode \n')
        if p['reciprocal']:
            kpoints.write('Reciprocal\n')
        else:
            kpoints.write('Cartesian\n')
        for n in range(len(p['kpts'])):
            [kpoints.write('%f ' % kpt) for kpt in p['kpts'][n]]
            if shape[1] == 4:
                kpoints.write('\n')
            elif shape[1] == 3:
                kpoints.write('1.0 \n')
    kpoints.close()
