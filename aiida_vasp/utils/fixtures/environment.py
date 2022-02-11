"""
Fixtures related to environments.

---------------------------------
Here we set up pytest fixtures that yield the correct AiiDA environment
to run our tests. In particular, the fresh_aiida_env is crucial in order
to set up a dummy database and configuration that can be used during
testing. AiiDA contributes with its own fixture manager, which we use.
"""
# pylint: disable=unused-import,unused-argument,redefined-outer-name, protected-access
from __future__ import absolute_import
from __future__ import print_function
from pathlib import Path
import pytest

from aiida.orm import QueryBuilder, Code, CalculationNode
from aiida.manage.tests import TemporaryProfileManager
from aiida.manage.tests.pytest_fixtures import aiida_caplog
from aiida.cmdline.utils.common import get_calcjob_report, get_workchain_report
from aiida.cmdline.utils.ascii_vis import format_call_graph
from aiida.cmdline.utils.common import get_node_info
from aiida.tools.importexport.dbexport import export


@pytest.fixture()
def fresh_aiida_env(aiida_profile, aiida_caplog):
    """Reset the database before and after the test function."""
    if isinstance(aiida_profile._manager, TemporaryProfileManager):
        print('The root directory of the fixture manager is: {}'.format(aiida_profile._manager.root_dir))  # pylint: disable=protected-access
    else:
        print('Using existing profile: {}'.format(aiida_profile._manager._profile.name))
    yield aiida_profile
    print_and_export_failed_mock()
    aiida_profile.reset_db()


def print_and_export_failed_mock():
    """
    Print details about any failed mock
    """

    query = QueryBuilder()
    query.append(CalculationNode,
                 filters={'or': [
                     {
                         'attributes.process_state': 'excepted'
                     },
                     {
                         'attributes.exit_status': {
                             '!==': 0
                         }
                     },
                 ]},
                 project=['*'])
    query.append(Code, with_outgoing=CalculationNode, filters={'extras.is_mock_code': True})
    if query.count() == 0:
        return

    # Print information
    print('######## Information for FAILED mock code calculations ########')
    entities = []
    for (calcjob,) in query.all():
        root = get_call_root(calcjob)
        if root != calcjob:
            print('######## Information for the call root ########')
            print(format_call_graph(root))
            print(get_node_info(root))
            print(get_workchain_report(root, 'REPORT'))

        print('######## Information for the calcjob ########')
        print(get_node_info(calcjob))
        print(get_calcjob_report(calcjob))
        names = calcjob.outputs.retrieved.list_object_names()
        stdout = None
        for option in ['stdout', 'vasp_output']:
            if option in names:
                stdout = option
                break
        if stdout is not None:
            print('######## STDOUT from the calcjob ########')
            print(calcjob.outputs.retrieved.get_object_content(stdout))
        else:
            print('ERROR: No STDOUT found for the calculation')

    # Export the failed calculations
    output_file = Path(__file__).parent.parent.parent.parent / 'test_mock_error.aiida'
    export(entities, filename=str(output_file), file_format='zip', include_logs=True, overwrite=True)


def get_call_root(node):
    """Obtain the root of the caller"""
    caller = node
    while True:
        next_caller = caller.caller

        if next_caller is None:
            break
        caller = next_caller
    return caller
