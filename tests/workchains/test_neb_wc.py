"""
Test for the NEB workchain
"""

from io import StringIO

import pytest
from aiida import orm
from aiida.engine import run_get_node
from aiida.plugins import DataFactory, WorkflowFactory

from aiida_vasp.parsers.content_parsers.poscar import PoscarParser
from aiida_vasp.parsers.vasp import get_structure_node
from aiida_vasp.utils.mock_code import VaspMockRegistry
from aiida_vasp.utils.neb import neb_interpolate


@pytest.fixture
def nh3_end_points(fresh_aiida_env):
    """Making NH3 structure for NEB example"""
    pos1 = StringIO(
        """#NO COMMENT
  1.00000000000000
    6.000000    0.000000    0.000000
    0.000000    7.000000    0.000000
    0.000000    0.000000    8.000000
    H     N
    3     1
Direct
 0.636428  0.567457  0.5491645
 0.500000  0.364985  0.5491330
 0.363572  0.567457  0.5491645
 0.500000  0.500000  0.5000000
"""
    )
    pos2 = StringIO(
        """ammonia flipping
  1.00000000000000
    6.000000    0.000000    0.000000
    0.000000    7.000000    0.000000
    0.000000    0.000000    8.000000
    H     N
    3     1
Direct
 0.636428  0.567457  0.4508355
 0.500000  0.364985  0.4508670
 0.363572  0.567457  0.4508355
 0.500000  0.500000  0.5000000
"""
    )
    init = PoscarParser(handler=pos1).structure
    final = PoscarParser(handler=pos2).structure

    return get_structure_node(init), get_structure_node(final)


@pytest.fixture
def neb_wc_input(fresh_aiida_env, upload_potcar, potcar_family_name, potcar_mapping, mock_vasp_strict, nh3_end_points):
    """Generate inputs for an NEB workchain"""
    # upload_real_pseudopotentials('/home/bonan/appdir/VASP/POTCARS/potpaw_PBE.54-2015_subset/')
    init, final = nh3_end_points
    neb_frames = neb_interpolate(init, final, orm.Int(3))
    parameters = {
        'images': 3,
        'spring': -5,
        'ibrion': 3,
        'nsw': 50,
        'algo': 'normal',
        'potim': 0.0,
        'iopt': 1,
        'ediffg': -0.02,
    }
    kpoints = orm.KpointsData()
    kpoints.set_kpoints_mesh((1, 1, 1))
    builder = WorkflowFactory('vasp.neb').get_builder()
    builder.parameters = orm.Dict(dict={'incar': parameters})
    builder.options = orm.Dict(dict={'resources': {'tot_num_mpiprocs': 1, 'num_machines': 1}, 'withmpi': False})

    builder.potential_family = orm.Str(potcar_family_name)
    builder.potential_mapping = orm.Dict(dict=potcar_mapping)

    builder.kpoints = kpoints
    builder.initial_structure = neb_frames['image_init']
    builder.final_structure = neb_frames['image_final']
    builder.neb_images = {f'image_{i:02d}': neb_frames[f'image_{i:02d}'] for i in (1, 2, 3)}
    builder.code = mock_vasp_strict
    #    create_authinfo(computer=mock_vasp_strict.computer, store=True)

    builder.metadata.label = 'NH3 NEB'

    return builder


def upload_real_pseudopotentials(path):
    """
    Upload real pseudopotentials for workchain test mock deposition


    This function should be called once before the REAL vasp calculation is launch to setup the
    correct POTCARs
    """
    global POTCAR_FAMILY_NAME  # noqa: PLW0603
    POTCAR_FAMILY_NAME = 'TEMP'
    potcar_data_cls = DataFactory('vasp.potcar')
    potcar_data_cls.upload_potcar_family(path, 'TEMP', 'TEMP-REALPOTCARS', stop_if_existing=False, dry_run=False)


def upload_real_workchain(node, name):
    """
    Upload the workchain to the repository to make it work with mocking

    This function should be called once after the REAL vasp calculation is run during the test
    """

    reg = VaspMockRegistry()
    print(reg.base_path)
    reg.upload_aiida_work(node, name)


def test_vasp_neb_wc(fresh_aiida_env, neb_wc_input):
    """Test the workchain"""

    _, out_node = run_get_node(neb_wc_input)
    assert out_node.exit_status == 0
    assert 'image_01' in out_node.outputs.structure
    assert 'image_02' in out_node.outputs.structure
    assert 'image_03' in out_node.outputs.structure
    assert 'total_energies' in out_node.outputs.misc
    assert 'forces' in out_node.outputs.misc
    # upload_real_workchain(out_node, "neb-workchain")
