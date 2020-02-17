"""
Test submitting a ConvergenceWorkChain.

Only `run` currently works.

"""
# pylint: disable=unused-import,wildcard-import,unused-wildcard-import,unused-argument,redefined-outer-name, too-many-statements, import-outside-toplevel
from __future__ import print_function

import numpy as np
import pytest
from aiida.common.extendeddicts import AttributeDict

from aiida_vasp.utils.fixtures import *
from aiida_vasp.utils.fixtures.data import POTCAR_FAMILY_NAME, POTCAR_MAP
from aiida_vasp.utils.fixtures.testdata import data_path
from aiida_vasp.utils.aiida_utils import get_data_node, get_current_user, aiida_version, cmp_version
from aiida_vasp.parsers.file_parsers.poscar import PoscarParser
from aiida_vasp.parsers.file_parsers.incar import IncarParser
from aiida_vasp.utils.aiida_utils import create_authinfo


@pytest.mark.skip(reason='Currently fails, need to find a better way to eject logs etc.')
@pytest.mark.wc
def test_bands_wc(fresh_aiida_env, potentials, mock_vasp):
    """Test with mocked vasp code."""
    from aiida.orm import Code, Log
    from aiida.plugins import WorkflowFactory
    from aiida.engine import run

    workchain = WorkflowFactory('vasp.bands')

    mock_vasp.store()
    create_authinfo(computer=mock_vasp.computer, store=True)

    structure = PoscarParser(file_path=data_path('test_bands_wc', 'inp', 'POSCAR')).structure
    parameters = IncarParser(file_path=data_path('test_bands_wc', 'inp', 'INCAR')).incar
    parameters['system'] = 'test-case:test_bands_wc'
    #chgcar = get_data_node('vasp.chargedensity', file=data_path('test_bands_wc', 'inp', 'CHGCAR'))

    inputs = AttributeDict()
    inputs.code = Code.get_from_string('mock-vasp@localhost')
    inputs.structure = structure
    inputs.parameters = get_data_node('dict', dict=parameters)
    #inputs.chgcar = chgcar
    inputs.potential_family = get_data_node('str', POTCAR_FAMILY_NAME)
    inputs.potential_mapping = get_data_node('dict', dict=POTCAR_MAP)
    inputs.options = get_data_node('dict',
                                   dict={
                                       'withmpi': False,
                                       'queue_name': 'None',
                                       'import_sys_environment': True,
                                       'resources': {
                                           'num_machines': 1,
                                           'num_mpiprocs_per_machine': 1
                                       },
                                       'max_wallclock_seconds': 3600
                                   })
    inputs.max_iterations = get_data_node('int', 1)
    inputs.clean_workdir = get_data_node('bool', False)
    inputs.verbose = get_data_node('bool', True)
    results, node = run.get_node(workchain, **inputs)
    print(Log.objects.get_logs_for(node))
    assert node.exit_status == 0
    assert 'bands' in results
    assert 'kpoints' in results
    kpoints = results['kpoints'].get_kpoints()
    np.set_printoptions(threshold=np.nan)
    np.testing.assert_allclose(kpoints[0, 0:3], np.array([0., 0., 0.]))
    np.testing.assert_allclose(kpoints[97, 0:3], np.array([0.5, 0., 0.5]))
    bands = results['bands'].get_bands()
    assert bands.shape == (1, 98, 20)
    np.testing.assert_allclose(bands[0, 0, 0:3], np.array([-6.0753, 6.0254, 6.0254]))
    np.testing.assert_allclose(bands[0, 2, 0:3], np.array([-6.0386, 5.7955, 5.8737]))
    np.testing.assert_allclose(bands[0, 97, 0:3], np.array([-1.867, -1.867, 3.1102]))
