__doc__ = '''Script to submit many vasp calculations
based on a single CIF file. Intended to test as many
types of VASP simulations and how well the plugin deals
with the variation in output files and formats.'''

# ~ from aiida_vasp.calcs.maker import VaspMaker
# ~ from aiida.orm import Code
# ~ from sys import argv

from_scratch = dict(istart=0)
continuing = dict(istart=1)

# check for out.structure, out.wavefunctions, out.charge_density
md = dict(ibrion=0, nsw=5)
relax = [dict(ibrion=i) for i in [1, 2, 3, 5, 6, 7, 8, 44]]

# check for out.wavefunctions and out.charge_density
static = dict(ibrion=-1, nsw=0)

# ~ modes = [md, relax[0], static]

# check for different PROCAR formats
procar = dict(lorbit=10)
lmprocar = dict(lorbit=11)
lmprocar_phase = dict(lorbit=12)
# ~ pcforms = [procar, lmprocar, lmprocar_phase]

# different combos should be different DOSCAR ouputs
noncol = dict(lsorbit=True, gga_compat=False)
spinpol = dict(ispin=2)

# should be continuations from standard dft, should not
# be able to do manual kpoint strings but only meshes
# should not influence output files
hy_pbe0 = dict(lhfcalc=True, gga='PE', precfock='Fast')
hy_hse03 = dict(lhfcalc=True, hfscreen=0.3, gga='PE', precfock='Fast')
hy_hse06 = dict(lhfcalc=True, hfscreen=0.2, gga='PE', precfock='Fast')
hy_b3lyp = dict(lhfcalc=True, hfscreen=0.2, gga='B3', aexx=0.2,
                aggax=0.72, aggac=0.81, aldac=0.19, precfock='Fast')
hy_hf = dict(lhfcalc=True, aexx=1.0, aggac=0.0, aldac=0.0, precfock='Fast')
# ~ funcs = [dict(), hy_hse06, hy_b3lyp, hy_hf]


def build_parameters(*args):
    parameters = {}
    for arg in args:
        parameters.update(arg)
    return parameters

# ~
# ~ if __name__ == '__main__':
    # ~ mkr = VaspMaker(structure=argv[0])
    # ~ mkr.kpoints.set_kpoints_mesh([2, 2, 2])
    # ~ mkr.code = Code.get_from_string('asevasp')
    # ~ mkr.queue = 'dphys_compute'
