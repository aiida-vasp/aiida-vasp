"""
Test submitting a RelaxWorkChain.

This does not seem to work, for `submit` the daemon will not pick up the workchain
and `run` just seems to get stuck after a while.
"""
# pylint: disable=unused-import,wildcard-import,unused-wildcard-import,unused-argument,redefined-outer-name,no-member, import-outside-toplevel
from __future__ import print_function

import pytest
from aiida.common.extendeddicts import AttributeDict

from aiida_vasp.utils.fixtures import *
from aiida_vasp.utils.fixtures.data import POTCAR_FAMILY_NAME, POTCAR_MAP
from aiida_vasp.utils.fixtures.testdata import data_path
from aiida_vasp.utils.aiida_utils import get_data_node, aiida_version, cmp_version, create_authinfo
from aiida_vasp.parsers.file_parsers.kpoints import KpointsParser
from aiida_vasp.parsers.file_parsers.poscar import PoscarParser
from aiida_vasp.parsers.file_parsers.incar import IncarParser


@pytest.mark.wc
def test_relax_wc(fresh_aiida_env, vasp_params, potentials, mock_vasp):
    # def test_relax_wc(fresh_aiida_env, vasp_params, potentials, mock_vasp, mock_relax_wc):
    """Test submitting only, not correctness, with mocked vasp code."""
    from aiida.orm import Code
    from aiida.plugins import WorkflowFactory
    from aiida.engine import run

    workchain = WorkflowFactory('vasp.relax')

    mock_vasp.store()
    create_authinfo(computer=mock_vasp.computer, store=True)

    structure = PoscarParser(file_path=data_path('test_relax_wc', 'inp', 'POSCAR')).structure
    kpoints = KpointsParser(file_path=data_path('test_relax_wc', 'inp', 'KPOINTS')).kpoints
    parameters = IncarParser(file_path=data_path('test_relax_wc', 'inp', 'INCAR')).incar
    parameters = {k: v for k, v in parameters.items() if k not in ['isif', 'ibrion', 'nsw', 'ediffg']}
    parameters['system'] = 'test-case:test_relax_wc'
    parameters['relax'] = {}
    parameters['relax']['perform'] = True
    parameters['relax']['algo'] = 'cg'
    parameters['relax']['force_cutoff'] = 0.01

    inputs = AttributeDict()
    inputs.code = Code.get_from_string('mock-vasp@localhost')
    inputs.structure = structure
    inputs.kpoints = kpoints
    inputs.parameters = get_data_node('dict', dict=parameters)
    inputs.potential_family = get_data_node('str', POTCAR_FAMILY_NAME)
    inputs.potential_mapping = get_data_node('dict', dict=POTCAR_MAP)
    inputs.options = get_data_node('dict',
                                   dict={
                                       'withmpi': False,
                                       'queue_name': 'None',
                                       'max_wallclock_seconds': 1,
                                       'import_sys_environment': True,
                                       'resources': {
                                           'num_machines': 1,
                                           'num_mpiprocs_per_machine': 1
                                       },
                                   })
    inputs.max_iterations = get_data_node('int', 1)
    inputs.clean_workdir = get_data_node('bool', False)
    inputs.verbose = get_data_node('bool', True)
    #raise ValueError(inputs.parameters.get_dict())
    results, node = run.get_node(workchain, **inputs)
    assert node.exit_status == 0
    assert 'relax' in results
    relax = results['relax']
    assert 'structure' in relax
    sites = relax['structure'].sites
    assert sites[0].kind_name == 'Si'
    assert sites[1].kind_name == 'Si'
    assert sites[0].position == (4.8125, 4.8125, 4.8125)
    assert sites[1].position == (0.6875, 0.6875, 0.715)


def test_mode_values():
    """Test that geometry relaxation modes either return a value or raise a ValueError"""
    from aiida_vasp.workchains.relax import RelaxWorkChain

    for pos in (True, False):
        for shape in (True, False):
            for volume in (True, False):
                try:
                    RelaxWorkChain.ModeEnum.get_from_dof(**{'positions': pos, 'shape': shape, 'volume': volume})
                except ValueError:
                    pass
                except Exception as exc:  # pylint: disable=broad-except
                    assert 'Get from DOF function has to either return the correct value or raise a ValueError ' \
                           'for invalid combinations, instead got {} exception'.format(type(exc))
