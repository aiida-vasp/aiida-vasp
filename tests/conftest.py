from __future__ import annotations

import os
import pathlib
import subprocess as sp

import numpy as np
import pytest
from aiida import orm
from aiida.cmdline.utils.ascii_vis import format_call_graph
from aiida.cmdline.utils.common import get_calcjob_report, get_node_info, get_workchain_report
from aiida.common.exceptions import NotExistent
from aiida.common.extendeddicts import AttributeDict
from aiida.engine import run
from aiida.orm import CalculationNode, Code, Computer, Dict, InstalledCode, QueryBuilder, load_code
from aiida.plugins import CalculationFactory, DataFactory, WorkflowFactory
from aiida.tools.archive import create_archive

from aiida_vasp.data.potcar import OLD_POTCAR_FAMILY_TYPE, Group, PotcarData, PotcarFileData, PotcarGroup, temp_potcar
from aiida_vasp.utils.general import copytree
from aiida_vasp.workchains.v2.common.builder_updater import VaspBuilderUpdater

pytest_plugins = 'aiida.tools.pytest_fixtures'

cwd = pathlib.Path(__file__).parent


@pytest.fixture()
def fresh_aiida_env(aiida_profile_clean):
    """Reset the database before and after the test function."""
    try:
        yield aiida_profile_clean
    finally:
        print_and_export_failed_mock()


@pytest.fixture(scope='session')
def localhost_dir(tmp_path_factory):
    return tmp_path_factory.mktemp('localhost_work')


@pytest.fixture()
def localhost(aiida_profile, localhost_dir):
    """Fixture for a local computer called localhost. This is currently not in the AiiDA fixtures."""
    try:
        computer = Computer.collection.get(label='localhost')
    except NotExistent:
        computer = Computer(
            label='localhost',
            hostname='localhost',
            transport_type='core.local',
            scheduler_type='core.direct',
            workdir=str(localhost_dir),
        ).store()
        computer.set_minimum_job_poll_interval(0.0)
        computer.configure()
    return computer


@pytest.fixture()
def mock_vasp(aiida_profile, localhost):
    """
    Give an mock-up of the VASP executable

    This code will always create the output object even if no matching
    calculations from the registry is found. This makes it suitable for simple
    tests.
    """
    return _create_mock_vasp_code(aiida_profile, localhost, 'mock-vasp-loose')


@pytest.fixture()
def mock_vasp_strict(aiida_profile, localhost):
    """
    Give an mock-up of the VASP executable with strict input matching.

    This code will not create the output object unless matching calculations from the
    registry is found. It is suitable for testsing complex multi-step workchains.
    tests.
    """
    return _create_mock_vasp_code(aiida_profile, localhost, 'mock-vasp')


@pytest.fixture()
def data_path():
    """Give the path to test data."""

    def _data_path(*args):
        path = pathlib.Path(cwd / 'test_data', *args)
        path = path.resolve()
        assert path.exists()
        assert path.is_absolute()
        return str(path)

    return _data_path


@pytest.fixture()
def read_file(data_path):
    """Give the content (string) of a test data file."""

    def _read_file(*args, **kwargs):
        path = kwargs.get('path', None)
        mode = kwargs.pop('mode', None)
        encoding = kwargs.pop('encoding', None)
        if not mode:
            mode = 'r'
        if not path:
            path = data_path(*args)
        if encoding is not None and 'b' not in mode:
            encoding = 'utf8'
        with open(path, mode, encoding=encoding) as testdata_fo:
            testdata_content = testdata_fo.read()
        return testdata_content

    return _read_file


@pytest.fixture()
def vasp_code(localhost):
    """Fixture for a vasp code, the executable it points to does not exist."""

    if not localhost.pk:
        localhost.store()
    code = InstalledCode(localhost, '/usr/local/bin/vasp')
    code.label = 'vasp'
    code.description = 'VASP code'
    code.default_calc_job_plugin = 'vasp.vasp'
    return code


@pytest.fixture
def vasp_params(aiida_profile):
    incar_data = orm.Dict(dict={'gga': 'PE', 'gga_compat': False, 'lorbit': 11, 'sigma': 0.5, 'magmom': '30 * 2*0.'})
    return incar_data


@pytest.fixture(params=['cif', 'str'])
def vasp_structure(request, aiida_profile, data_path):
    """Fixture: StructureData or CifData."""
    if request.param == 'cif':
        cif_path = data_path('cif', 'EntryWithCollCode43360.cif')
        structure = DataFactory('core.cif').get_or_create(cif_path)[0]
    elif request.param == 'str':
        larray = np.array([[0, 0.5, 0.5], [0.5, 0, 0.5], [0.5, 0.5, 0]])
        alat = 6.058
        structure = DataFactory('core.structure')(cell=larray * alat)
        structure.append_atom(position=[0, 0, 0], symbols='In')
        structure.append_atom(position=[0.25, 0.25, 0.25], symbols='As')
        structure.append_atom(position=[0.25, 0.33, 0.34], symbols='As')
        structure.append_atom(position=[0.5, 0.5, 0.5], symbols='In', name='In_d')
        structure.append_atom(position=[0.7896, 0.6234, 0.5], symbols='In', name='In_d')
        structure.append_atom(position=[0.75, 0.75, 0.75], symbols='As')
    elif request.param == 'str-Al':
        larray = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
        alat = 4.04
        structure = DataFactory('core.structure')(cell=larray * alat)
        structure.append_atom(position=np.array([0, 0, 0]) * alat, symbols='Al')
        structure.append_atom(position=np.array([0, 0.5, 0.5]) * alat, symbols='Al')
        structure.append_atom(position=np.array([0.5, 0, 0.5]) * alat, symbols='Al')
        structure.append_atom(position=np.array([0.5, 0.5, 0]) * alat, symbols='Al')
    elif request.param == 'str-InAs':
        structure_cls = DataFactory('core.structure')
        structure = structure_cls(cell=np.array([[0, 0.5, 0.5], [0.5, 0, 0.5], [0.5, 0.5, 0]]) * 6.058)
        structure.append_atom(position=(0, 0, 0), symbols='In', name='Hamburger')
        structure.append_atom(position=(0.25, 0.25, 0.25), symbols='As', name='Pizza')
    return structure


@pytest.fixture(params=['mesh', 'list'])
def vasp_kpoints(request, aiida_profile, data_path):
    """Fixture: (kpoints object, resulting KPOINTS)."""

    def _ref_kp_list():
        with open(data_path('kpoints', 'KPOINTS_list'), 'r', encoding='utf8') as reference_kpoints_fo:
            ref_kp_str = reference_kpoints_fo.read()
        return ref_kp_str

    def _ref_kp_mesh():
        with open(data_path('kpoints', 'KPOINTS_mesh'), 'r', encoding='utf8') as reference_kpoints_fo:
            ref_kp_list = reference_kpoints_fo.read()
        return ref_kp_list

    if request.param == 'mesh':
        kpoints = DataFactory('core.array.kpoints')()
        kpoints.set_kpoints_mesh([2, 2, 2])
        ref_kpoints = _ref_kp_mesh()
    elif request.param == 'list':
        kpoints = DataFactory('core.array.kpoints')()
        kpoints.set_kpoints([[0.0, 0.0, 0.0], [0.0, 0.0, 0.5]], weights=[1.0, 1.0])
        ref_kpoints = _ref_kp_list()
    return kpoints, ref_kpoints


@pytest.fixture(scope='session')
def potcar_family_name() -> str:
    """Return the POTCAR family name."""
    _potcar_family_name = 'test_family'
    return _potcar_family_name


@pytest.fixture(scope='session')
def potcar_mapping() -> dict:
    """Return the POTCAR mapping."""
    _potcar_mapping = {
        'In': 'In_sv',
        'In_d': 'In_d',
        'As': 'As',
        'Ga': 'Ga',
        'Si': 'Si',
        'P': 'P',
        'S': 'S',
        'Zn': 'Zn',
        'N': 'N',
        'H': 'H',
    }
    return _potcar_mapping


@pytest.fixture
def temp_pot_folder(tmp_path, data_path):
    """A temporary copy of the potcar test data folder, to avoid extracting tar objects inside the repo."""
    potcar_ga = pathlib.Path(data_path('potcar')) / 'Ga'
    assert not potcar_ga.exists()
    pot_archive = pathlib.Path(data_path('potcar'))
    target = tmp_path / 'potentials'
    # Ensure that the target path exists
    pathlib.Path(target).mkdir(exist_ok=True)
    copytree(pot_archive, target)
    return target


@pytest.fixture
def upload_potcar(aiida_profile, temp_pot_folder, potcar_family_name, data_path):
    """Upload a POTCAR family to DB."""
    if PotcarData.get_potcar_group(potcar_family_name) is not None:
        # Already uploaded, so we do not need to do it again.
        return

    def duplicate_potcar_data(potcar_node):
        """Create and store (and return) a duplicate of a given PotcarData node."""
        file_node = PotcarFileData()
        with temp_potcar(potcar_node.get_content()) as potcar_file:
            file_node.add_file(potcar_file)
            file_node.base.attributes.set('sha512', 'abcd')
            file_node.base.attributes.set('full_name', potcar_node.full_name)
            file_node.store()
        data_node, _ = PotcarData.get_or_create(file_node)
        return data_node

    potcar_ga = pathlib.Path(data_path('potcar')) / 'Ga'
    family_name = potcar_family_name
    family_desc = 'A POTCAR family used as a test fixture. Contains only unusable POTCAR files.'
    potcar_cls = PotcarData
    potcar_cls.upload_potcar_family(str(temp_pot_folder), family_name, family_desc, stop_if_existing=False)
    if len(potcar_cls.find(full_name='In_d')) == 1:
        family_group = potcar_cls.get_potcar_group(potcar_family_name)
        in_d = potcar_cls.find(full_name='In_d')[0]
        in_d_double = duplicate_potcar_data(in_d)
        family_group.add_nodes(in_d_double)
        assert in_d.uuid == potcar_cls.find(full_name='In_d')[0].uuid
    assert 'As' in potcar_cls.get_full_names(potcar_family_name, 'As')
    assert 'Ga' in potcar_cls.get_full_names(potcar_family_name, 'Ga')
    assert 'In_d' in potcar_cls.get_full_names(potcar_family_name, 'In')
    assert not potcar_ga.exists()


@pytest.fixture
def potentials(upload_potcar, potcar_family_name, potcar_mapping):
    """Fixture for two incomplete POTPAW potentials."""
    potcar_cls = PotcarData
    potentials = potcar_cls.get_potcars_dict(
        ['In', 'In_d', 'As'], family_name=potcar_family_name, mapping=potcar_mapping
    )

    return potentials


@pytest.fixture
def legacy_potcar_family(upload_potcar, potcar_family_name):
    """
    Fixture from creating an legacy potcar group

    Returns a tuple of group label and the LegacyPotcarGroup with the old type_string
    """

    class LegacyPotcarGroup(Group):
        """Old style group with the old type string"""

    # Override the _type_string class property which is supposed to be loaded from the entrypoint.
    LegacyPotcarGroup._type_string = OLD_POTCAR_FAMILY_TYPE
    new_group = PotcarGroup.collection.get(label=potcar_family_name)
    old_group = LegacyPotcarGroup(label=potcar_family_name + '_migrate_test')
    old_group.store()

    # Add the nodes from the new group to the old group
    old_group.add_nodes(list(new_group.nodes))
    return old_group.label, LegacyPotcarGroup


@pytest.fixture
def phonondb_run(tmp_path, data_path):
    phonondb = data_path('phonondb')
    copytree(phonondb, tmp_path)
    yield tmp_path


@pytest.fixture()
def run_vasp_process(
    aiida_profile,
    upload_potcar,
    vasp_params,
    potcar_family_name,
    potcar_mapping,
    vasp_kpoints,
    vasp_structure,
    mock_vasp,
):
    """Setup a standard VaspCalculation or VaspWorkChain with the mock executable that accepts input overrides."""

    def inner(inputs=None, settings=None, test_case=None, process_type='calcjob'):
        """
        Run a VaspCalculation or VaspWorkChain with specified input and settings overrides.

        Specific outputs can be selected using the test_case parameter.

        The type of process is set with the process_type parameter.
        """
        inpts = AttributeDict()
        inpts.structure = vasp_structure
        parameters = vasp_params.get_dict()
        options = {
            'withmpi': False,
            'queue_name': 'None',
            'resources': {'num_machines': 1, 'num_mpiprocs_per_machine': 1},
            'max_wallclock_seconds': 3600,
        }
        if test_case is not None:
            # Allow to fetch special tests cases using the mock-vasp executable
            parameters['system'] = f'test-case:{test_case}'
        if process_type == 'calcjob':
            process = CalculationFactory('vasp.vasp')
            inpts.potential = PotcarData.get_potcars_from_structure(
                structure=inpts.structure,
                family_name=potcar_family_name,
                mapping=potcar_mapping,
            )
            inpts.parameters = orm.Dict(dict=parameters)
            inpts.metadata = {}
            inpts.metadata['options'] = options
        elif process_type == 'workchain':
            process = WorkflowFactory('vasp.vasp')
            inpts.potential_family = orm.Str(potcar_family_name)
            inpts.potential_mapping = orm.Dict(dict=potcar_mapping)
            inpts.parameters = orm.Dict(dict={'incar': parameters})
            inpts.options = orm.Dict(dict=options)
            inpts.max_iterations = orm.Int(1)
            inpts.clean_workdir = orm.Bool(False)
            inpts.verbose = orm.Bool(True)
        else:
            raise ValueError(
                f"The supplied process_type: {process_type} is not supported. Use either 'calcjob' or 'workchain.'"
            )

        mock_vasp.store()
        # create_authinfo(computer=mock_vasp.computer, store=True)
        inpts.code = load_code('mock-vasp-loose@localhost')
        kpoints, _ = vasp_kpoints
        inpts.kpoints = kpoints
        if inputs is not None:
            # Allow overrides of the input
            inpts.update(inputs)
        if settings is not None and isinstance(settings, dict):
            inpts.settings = Dict(dict=settings)
        results_and_node = run.get_node(process, **inpts)
        return results_and_node

    return inner


@pytest.fixture()
def mock_potcars(aiida_profile, temp_pot_folder):
    """Create family of potcars for mock-vasp"""
    path = os.environ.get('MOCK_VASP_POTCAR_PATH', None)
    if path:
        PotcarData.upload_potcar_family(path, 'PBE.54', 'Family for mock calculation', stop_if_existing=False)
        return 'PBE.54'
    else:
        PotcarData.upload_potcar_family(
            str(temp_pot_folder), 'PBE.54', 'Family for mock calculation', stop_if_existing=False
        )
    return None


@pytest.fixture()
def builder_updater(
    aiida_profile,
    mock_potcars,
    mock_vasp,
):
    """
    Return a Builder Updater object for mock-vasp
    """
    return VaspBuilderUpdater(code='mock-vasp@localhost')


def print_and_export_failed_mock():
    """
    Print details about any failed mock
    """

    query = QueryBuilder()
    query.append(
        CalculationNode,
        filters={
            'or': [
                {'attributes.process_state': 'excepted'},
                {'attributes.exit_status': {'!==': 0}},
            ],
        },
        project=['*'],
    )
    query.append(Code, with_outgoing=CalculationNode, filters={'extras.is_mock_code': True})
    if query.count() == 0:
        return

    # Print information
    print('######## Information for FAILED mock code calculations ########')
    entities = []
    for (calcjob,) in query.all():
        entities.append(calcjob)
        root = get_call_root(calcjob)
        if root != calcjob:
            print('######## Information for the call root ########')
            print(format_call_graph(root))
            print(get_node_info(root))
            print(get_workchain_report(root, 'REPORT'))

        print('######## Information for the calcjob ########')
        print(get_node_info(calcjob))
        print(get_calcjob_report(calcjob))
        names = calcjob.outputs.retrieved.base.repository.list_object_names()
        stdout = None
        for option in ['stdout', 'vasp_output']:
            if option in names:
                stdout = option
                break
        if stdout is not None:
            print('######## STDOUT from the calcjob ########')
            print(calcjob.outputs.retrieved.base.repository.get_object_content(stdout))
        else:
            print('ERROR: No STDOUT found for the calculation')

    # Export the failed calculations
    output_file = cwd / 'test_mock_error.aiida'
    create_archive(entities, filename=str(output_file), include_logs=True, overwrite=True)


def get_call_root(node):
    """Obtain the root of the caller"""
    caller = node
    while True:
        next_caller = caller.caller

        if next_caller is None:
            break
        caller = next_caller
    return caller


def _create_mock_vasp_code(aiida_profile, localhost, exec_name):
    """
    Points to a mock-up of a VASP executable.

    If environmental variable REAL_VASP_PATH is set, switch the code
    to point to the REAL VASP executable. This is used to generate the
    actual outputs for mock tests later
    """
    query_builder = QueryBuilder()
    query_builder.append(Code, tag='code')
    query_builder.add_filter('code', {'label': {'==': exec_name}})
    query_results = query_builder.all()
    if query_results:
        code = query_results[0][0]
    else:
        os_env = os.environ.copy()
        if not localhost.pk:
            localhost.store()
        # returns unicode
        mock_vasp_path = sp.check_output(['which', exec_name], env=os_env, universal_newlines=True).strip()

        # Allow overriding mock using REAL code, this is used for running the actual
        # calculation once and deposit the results in the registry
        if os.environ.get('REAL_VASP_PATH'):
            mock_vasp_path = os.environ['REAL_VASP_PATH']

        code = InstalledCode(localhost, mock_vasp_path)
        code.label = exec_name
        code.description = 'Mock VASP for tests'
        code.default_calc_job_plugin = 'vasp.vasp'
        code.store()
        code.base.extras.set('is_mock_code', True)

    return code
