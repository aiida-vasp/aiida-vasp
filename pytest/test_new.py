import sys

def test_v5(configure_with_daemon, sample):
    from aiida.orm import Code, Computer
    from aiida.orm.querytool import QueryTool
    from aiida_vasp.calcs.maker import VaspMaker

    mkcalc = VaspMaker(structure=sample('GaAs/POSCAR'))
    mkcalc.code = Code.get_from_string('vasp')
    mkcalc.kpoints.set_kpoints_mesh([8, 8, 8])
    mkcalc.add_settings(
        system=mkcalc.structure.filename,
        npar=8
    )
    mkcalc.recipe = 'test_sc'

    mkcalc.computer = Computer.get('monch')
    mkcalc.queue = 'dphys_compute'

    v5 = mkcalc.new()
    v5.set_resources({
        'num_machines': 8,
        'num_mpiprocs_per_machine': 2})
    v5.run()
