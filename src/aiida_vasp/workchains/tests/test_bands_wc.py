"""
Test submitting a ConvergenceWorkChain.

Only `run` currently works.

"""

# pylint: disable=unused-import,wildcard-import,unused-wildcard-import,unused-argument,redefined-outer-name, too-many-statements, import-outside-toplevel
from __future__ import print_function

import numpy as np
from aiida.common.extendeddicts import AttributeDict
from aiida_vasp.parsers.content_parsers.incar import IncarParser
from aiida_vasp.parsers.content_parsers.poscar import PoscarParser
from aiida_vasp.parsers.node_composer import NodeComposer
from aiida_vasp.utils.aiida_utils import create_authinfo, get_data_node
from aiida_vasp.utils.fixtures import *
from aiida_vasp.utils.fixtures.data import POTCAR_FAMILY_NAME, POTCAR_MAP
from aiida_vasp.utils.fixtures.testdata import data_path


# @pytest.mark.skip(reason='Currently fails, need to find a better way to eject logs etc.')
# @pytest.mark.wc
def test_bands_wc(fresh_aiida_env, potentials, mock_vasp):
    """Test with mocked vasp code."""
    from aiida.engine import run
    from aiida.orm import (
        Code,
        RemoteData,  # pylint: disable=no-name-in-module
    )
    from aiida.plugins import WorkflowFactory

    workchain = WorkflowFactory('vasp.bands')

    mock_vasp.store()
    create_authinfo(computer=mock_vasp.computer, store=True)

    structure = None
    with open(data_path('test_bands_wc', 'inp', 'POSCAR'), 'r', encoding='utf8') as handler:
        structure_parser = PoscarParser(handler=handler)
        structure = structure_parser.get_quantity('poscar-structure')
        structure = NodeComposer.compose_core_structure('core.structure', {'structure': structure_parser.structure})
    parameters = None
    with open(data_path('test_bands_wc', 'inp', 'INCAR'), 'r', encoding='utf8') as handler:
        incar_parser = IncarParser(handler=handler)
        parameters = incar_parser.get_quantity('incar')

    parameters['system'] = 'test-case:test_bands_wc'
    # Make sure we replace encut with pwcutoff
    del parameters['encut']
    parameters = {'incar': parameters}
    parameters['electronic'] = {'pwcutoff': 200}

    inputs = AttributeDict()
    inputs.code = Code.get_from_string('mock-vasp@localhost')
    inputs.structure = structure
    inputs.parameters = get_data_node('core.dict', dict=parameters)
    inputs.potential_family = get_data_node('core.str', POTCAR_FAMILY_NAME)
    inputs.potential_mapping = get_data_node('core.dict', dict=POTCAR_MAP)
    inputs.options = get_data_node(
        'core.dict',
        dict={
            'withmpi': False,
            'queue_name': 'None',
            'resources': {'num_machines': 1, 'num_mpiprocs_per_machine': 1},
            'max_wallclock_seconds': 3600,
        },
    )
    inputs.max_iterations = get_data_node('core.int', 1)
    inputs.clean_workdir = get_data_node('core.bool', False)
    inputs.verbose = get_data_node('core.bool', True)
    # Also set the restart folder as we assume a bands data will start from
    # a previous calculation that is sitting in the restart folder
    inputs.restart_folder = RemoteData(computer=inputs.code.computer, remote_path=data_path('test_bands_wc', 'inp'))
    results, node = run.get_node(workchain, **inputs)
    assert node.exit_status == 0
    assert 'bands' in results
    kpoints = results['bands'].get_kpoints()
    test_array = np.array(
        [
            [0.0, 0.0, 0.0],
            [0.02272727, 0.0, 0.02272727],
            [0.04545454, 0.0, 0.04545454],
            [0.06818182, 0.0, 0.06818182],
            [0.09090909, 0.0, 0.09090909],
            [0.11363636, 0.0, 0.11363636],
            [0.13636364, 0.0, 0.13636364],
            [0.15909091, 0.0, 0.15909091],
            [0.18181818, 0.0, 0.18181818],
            [0.20454545, 0.0, 0.20454545],
            [0.22727273, 0.0, 0.22727273],
            [0.25, 0.0, 0.25],
            [0.27272727, 0.0, 0.27272727],
            [0.29545455, 0.0, 0.29545455],
            [0.31818182, 0.0, 0.31818182],
            [0.34090909, 0.0, 0.34090909],
            [0.36363636, 0.0, 0.36363636],
            [0.38636364, 0.0, 0.38636364],
            [0.40909091, 0.0, 0.40909091],
            [0.43181818, 0.0, 0.43181818],
            [0.45454545, 0.0, 0.45454545],
            [0.47727273, 0.0, 0.47727273],
            [0.5, 0.0, 0.5],
            [0.51785714, 0.03571429, 0.51785714],
            [0.53571429, 0.07142857, 0.53571429],
            [0.55357143, 0.10714286, 0.55357143],
            [0.57142857, 0.14285714, 0.57142857],
            [0.58928571, 0.17857143, 0.58928571],
            [0.60714286, 0.21428571, 0.60714286],
            [0.625, 0.25, 0.625],
            [0.375, 0.375, 0.75],
            [0.35869565, 0.35869565, 0.7173913],
            [0.3423913, 0.3423913, 0.68478261],
            [0.32608696, 0.32608696, 0.65217391],
            [0.30978261, 0.30978261, 0.61956522],
            [0.29347826, 0.29347826, 0.58695652],
            [0.27717391, 0.27717391, 0.55434783],
            [0.26086957, 0.26086957, 0.52173913],
            [0.24456522, 0.24456522, 0.48913043],
            [0.22826087, 0.22826087, 0.45652174],
            [0.21195652, 0.21195652, 0.42391304],
            [0.19565217, 0.19565217, 0.39130435],
            [0.17934783, 0.17934783, 0.35869565],
            [0.16304348, 0.16304348, 0.32608696],
            [0.14673913, 0.14673913, 0.29347826],
            [0.13043478, 0.13043478, 0.26086957],
            [0.11413044, 0.11413044, 0.22826087],
            [0.09782609, 0.09782609, 0.19565217],
            [0.08152174, 0.08152174, 0.16304348],
            [0.06521739, 0.06521739, 0.13043478],
            [0.04891304, 0.04891304, 0.09782609],
            [0.0326087, 0.0326087, 0.06521739],
            [0.01630435, 0.01630435, 0.0326087],
            [0.0, 0.0, 0.0],
            [0.02631579, 0.02631579, 0.02631579],
            [0.05263158, 0.05263158, 0.05263158],
            [0.07894737, 0.07894737, 0.07894737],
            [0.10526316, 0.10526316, 0.10526316],
            [0.13157895, 0.13157895, 0.13157895],
            [0.15789474, 0.15789474, 0.15789474],
            [0.18421053, 0.18421053, 0.18421053],
            [0.21052632, 0.21052632, 0.21052632],
            [0.2368421, 0.2368421, 0.2368421],
            [0.26315789, 0.26315789, 0.26315789],
            [0.28947368, 0.28947368, 0.28947368],
            [0.31578947, 0.31578947, 0.31578947],
            [0.34210526, 0.34210526, 0.34210526],
            [0.36842105, 0.36842105, 0.36842105],
            [0.39473684, 0.39473684, 0.39473684],
            [0.42105263, 0.42105263, 0.42105263],
            [0.44736842, 0.44736842, 0.44736842],
            [0.47368421, 0.47368421, 0.47368421],
            [0.5, 0.5, 0.5],
            [0.5, 0.48333333, 0.51666667],
            [0.5, 0.46666667, 0.53333333],
            [0.5, 0.45, 0.55],
            [0.5, 0.43333333, 0.56666667],
            [0.5, 0.41666667, 0.58333333],
            [0.5, 0.4, 0.6],
            [0.5, 0.38333333, 0.61666667],
            [0.5, 0.36666667, 0.63333333],
            [0.5, 0.35, 0.65],
            [0.5, 0.33333333, 0.66666667],
            [0.5, 0.31666667, 0.68333333],
            [0.5, 0.3, 0.7],
            [0.5, 0.28333333, 0.71666667],
            [0.5, 0.26666667, 0.73333333],
            [0.5, 0.25, 0.75],
            [0.5, 0.225, 0.725],
            [0.5, 0.2, 0.7],
            [0.5, 0.175, 0.675],
            [0.5, 0.15, 0.65],
            [0.5, 0.125, 0.625],
            [0.5, 0.1, 0.6],
            [0.5, 0.075, 0.575],
            [0.5, 0.05, 0.55],
            [0.5, 0.025, 0.525],
            [0.5, 0.0, 0.5],
        ]
    )
    np.testing.assert_allclose(kpoints, test_array)
    bands = results['bands'].get_bands()
    assert bands.shape == (1, 98, 20)
    np.testing.assert_allclose(bands[0, 0, 0:3], np.array([-6.0753, 6.0254, 6.0254]))
    np.testing.assert_allclose(bands[0, 2, 0:3], np.array([-6.0386, 5.7955, 5.8737]))
    np.testing.assert_allclose(bands[0, 97, 0:3], np.array([-1.867, -1.867, 3.1102]))

    assert results['bands'].labels
