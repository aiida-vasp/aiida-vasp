"""
Test for the NEB workchain
"""
# pylint: disable=unused-import,wildcard-import,unused-wildcard-import,unused-argument,redefined-outer-name, import-outside-toplevel

from io import StringIO
import pytest

from aiida import orm
from aiida.plugins import WorkflowFactory

from aiida_vasp.utils.fixtures import *

from aiida_vasp.parsers.content_parsers.poscar import PoscarParser
from aiida_vasp.utils.fixtures.data import POTCAR_FAMILY_NAME, POTCAR_MAP
from aiida_vasp.utils.aiida_utils import create_authinfo
from aiida_vasp.utils.neb import neb_interpolate
from aiida_vasp.parsers.node_composer import NodeComposer


@pytest.fixture
def nh3_end_points(fresh_aiida_env):
    """Making NH3 structure for NEB example"""
    pos1 = StringIO("""#NO COMMENT
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
""")
    pos2 = StringIO("""ammonia flipping
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
""")
    init = PoscarParser(handler=pos1).structure
    final = PoscarParser(handler=pos2).structure

    init_structure = NodeComposer.compose_structure('structure', {'structure': init})
    final_structure = NodeComposer.compose_structure('structure', {'structure': final})
    return init_structure, final_structure


@pytest.fixture
def neb_wc_input(fresh_aiida_env, potentials, mock_vasp_strict, nh3_end_points):
    """Generate inputs for an NEB workchain"""
    #upload_real_pseudopotentials('/home/bonan/appdir/VASP/POTCARS/potpaw_PBE.54-2015_subset/')
    init, final = nh3_end_points
    neb_frames = neb_interpolate(init, final, orm.Int(3))
    parameters = {
        'images': 3,
        'spring': -5,
        'ibrion': 3,
        'nsw': 50,
        'algo': 'normal',
        'potim': 0.,
        'iopt': 1,
        'ediffg': -0.02,
    }
    kpoints = orm.KpointsData()
    kpoints.set_kpoints_mesh((1, 1, 1))
    builder = WorkflowFactory('vasp.neb').get_builder()
    builder.parameters = orm.Dict(dict={'incar': parameters})
    builder.options = orm.Dict(dict={'resources': {'tot_num_mpiprocs': 1, 'num_machines': 1}, 'withmpi': False})

    builder.potential_family = orm.Str(POTCAR_FAMILY_NAME)
    builder.potential_mapping = orm.Dict(dict=POTCAR_MAP)

    builder.kpoints = kpoints
    builder.initial_structure = neb_frames['image_init']
    builder.final_structure = neb_frames['image_final']
    builder.neb_images = {f'image_{i:02d}': neb_frames[f'image_{i:02d}'] for i in (1, 2, 3)}
    builder.code = mock_vasp_strict
    create_authinfo(computer=mock_vasp_strict.computer, store=True)

    builder.metadata.label = 'NH3 NEB'

    return builder


def upload_real_pseudopotentials(path):
    """
    Upload real pseudopotentials for workchain test mock deposition


    This function should be called once before the REAL vasp calculation is launch to setup the
    correct POTCARs
    """
    from aiida.plugins import DataFactory
    global POTCAR_FAMILY_NAME  # pylint: disable=global-statement
    POTCAR_FAMILY_NAME = 'TEMP'
    potcar_data_cls = DataFactory('vasp.potcar')
    potcar_data_cls.upload_potcar_family(path, 'TEMP', 'TEMP-REALPOTCARS', stop_if_existing=False, dry_run=False)


def upload_real_workchain(node, name):
    """
    Upload the workchain to the repository to make it work with mocking

    This function should be called once after the REAL vasp calculation is run during the test
    """
    from aiida_vasp.utils.mock_code import VaspMockRegistry
    reg = VaspMockRegistry()
    print(reg.base_path)
    reg.upload_aiida_work(node, name)


def test_vasp_neb_wc(fresh_aiida_env, neb_wc_input):
    """Test the workchain"""

    from aiida.engine import run_get_node
    _, out_node = run_get_node(neb_wc_input)
    assert out_node.exit_status == 0
    #upload_real_workchain(out_node, "neb-workchain")
