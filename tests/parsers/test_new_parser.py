import pathlib

import numpy as np
import pytest
from aiida import orm

from aiida_vasp.parsers.vasp import NotificationComposer, ParserSettingsConfig, VaspParser, is_all_empty


def test_all_empty_check():
    """Test the is_all_empty function"""
    empty = {'A': [], 'B': {'C': {}}}
    assert is_all_empty(empty)
    empty = []
    assert is_all_empty(empty)
    empty = {'A': []}
    assert is_all_empty(empty)
    not_empty = {'A': {'B': [1]}}
    assert not is_all_empty(not_empty)
    not_empty = {'A': [1, 2], 'B': 2}
    assert not is_all_empty(not_empty)
    assert not is_all_empty([1])


def test_parser_bare(calc_with_retrieved, request):
    _relative_file_path = '../test_data/basic_run'
    file_path = str(pathlib.Path(request.fspath).parent / _relative_file_path)
    node = calc_with_retrieved(file_path, {})
    parser = VaspParser(node)
    parser.parse(retrieved_tempoary_folder=file_path)
    assert 'misc' in parser.outputs
    assert 'site_magnetization' not in parser.outputs.misc.get_dict()
    assert 'structure' in parser.outputs

    node = calc_with_retrieved(file_path, {'parser_settings': {'include_quantity': ['projectors']}})
    parser = VaspParser(node)
    parser.parse(retrieved_tempoary_folder=file_path)
    assert 'dielectrics`' not in parser.outputs
    assert 'born_charges' not in parser.outputs
    assert 'parameters' not in parser.outputs['misc']
    assert 'trajectory' not in parser.outputs
    assert 'energies' not in parser.outputs

    # Test the parameters outputs
    node = calc_with_retrieved(file_path, {'parser_settings': {'include_quantity': ['parameters']}})
    parser = VaspParser(node)
    parser.parse(retrieved_tempoary_folder=file_path)
    assert 'parameters' in parser.outputs['misc']
    assert 'system' in parser.outputs['misc']['parameters']


@pytest.fixture
def parser_with_retrieved(calc_with_retrieved, request):
    """Fixture to create a VaspParser instance with a given pre-computed data folder"""

    def wrapped(name, settings={}, parse=True):
        _relative_file_path = f'../test_data/{name}'
        file_path = str(pathlib.Path(request.fspath).parent / _relative_file_path)
        node = calc_with_retrieved(file_path, settings)
        parser = VaspParser(node)
        exit_code = None
        if parse:
            exit_code = parser.parse(retrieved_tempoary_folder=file_path)
        return parser, exit_code

    return wrapped


@pytest.fixture
def parser_with_vasprun(parser_with_retrieved):
    def wrapped(settings={}):
        return parser_with_retrieved('vasprun', settings)

    return wrapped


@pytest.fixture
def parser_incomplete_output(parser_with_retrieved):
    """Fixture for tests with incomplete output"""
    default_settings = {
        'parser_settings': {
            'check_completeness': False,
            'critical_objects': [],
            'required_quantity': [],
        }
    }

    def wrapped(case, settings={}):
        if 'parser_settings' in settings:
            default_settings['parser_settings'].update(settings['parser_settings'])
        return parser_with_retrieved(case, default_settings)

    return wrapped


@pytest.fixture
def parser_with_outcar(parser_with_retrieved):
    def wrapped(settings={}):
        return parser_with_retrieved('outcar', settings)

    return wrapped


def test_parser_born(parser_with_retrieved):
    """Test parsing born effective charges"""
    parser, _ = parser_with_retrieved(
        'born_effective_charge',
        {
            'parser_settings': {
                'include_quantity': ['born_charges'],
                'check_completeness': False,
                'critical_objects': [],
            }
        },
    )
    assert 'born_charges' in parser.outputs


def test_parser_dielectrics(parser_incomplete_output):
    """Test parsing dielectrics"""
    parser, _ = parser_incomplete_output('dielectric')
    assert 'dielectrics' in parser.outputs


def test_parser_magnetization(parser_incomplete_output):
    """Test parsing dielectrics"""
    parser, _ = parser_incomplete_output('magnetization')
    assert 'misc' in parser.outputs
    assert parser.outputs['misc']['magnetization'][0] == 6.4424922


def test_parser_lepsilon(parser_incomplete_output):
    """Test parsing dielectrics"""
    parser, _ = parser_incomplete_output('lepsilon')
    assert 'misc' in parser.outputs
    assert 'epsilon' in parser.outputs['dielectrics'].get_arraynames()
    assert 'born_charges' in parser.outputs['born_charges'].get_arraynames()


def test_parser_localfield(parser_incomplete_output):
    """Test parsing dielectrics"""
    parser, _ = parser_incomplete_output('localfield')
    assert 'misc' in parser.outputs
    assert 'epsilon' in parser.outputs['dielectrics'].get_arraynames()
    assert 'born_charges' in parser.outputs['born_charges'].get_arraynames()


def test_parser_disp_details(parser_with_retrieved):
    """Test parsing elastic moduli and symmetries from OUTCAR"""
    parser, _ = parser_with_retrieved(
        'disp_details',
        {
            'parser_settings': {
                'include_quantity': ['symmetries', 'elastic_moduli'],
                'check_completeness': False,
                'critical_objects': [],
            }
        },
    )
    assert 'symmetries' in parser.outputs['misc']
    assert 'elastic_moduli' in parser.outputs['misc']
    data_dict = parser.outputs['misc'].get_dict()
    # Test
    test = np.array([1674.5786, 704.739, 704.739, -0.0, 0.0, 0.0])
    np.testing.assert_allclose(data_dict['elastic_moduli']['symmetrized'][0], test)
    test = np.array([0.0, 0.0, 0.0, -0.0, -0.0, 1122.6622])
    np.testing.assert_allclose(data_dict['elastic_moduli']['symmetrized'][5], test)
    test = np.array([705.0238, 1674.8491, 705.0238, -0.0, -0.0, 0.0])
    np.testing.assert_allclose(data_dict['elastic_moduli']['non_symmetrized'][1], test)
    test = np.array([-0.0078, -0.0495, 0.0147, 0.0, 1123.0829, -0.0])
    np.testing.assert_allclose(data_dict['elastic_moduli']['non_symmetrized'][4], test)
    test = np.array([704.739, 704.739, 1674.5786, -0.0, -0.0, 0.0])
    np.testing.assert_allclose(data_dict['elastic_moduli']['total'][2], test)
    test = np.array([-0.0, -0.0, -0.0, 775.8054, 0.0, -0.0])
    np.testing.assert_allclose(data_dict['elastic_moduli']['total'][3], test)

    assert data_dict['run_stats']
    assert data_dict['run_stats']['total_cpu_time_used'] == pytest.approx(89.795)
    assert data_dict['run_stats']['average_memory_used'] == pytest.approx(0.0)

    assert data_dict['run_status']['last_iteration_index'] == [15, 5]
    assert data_dict['run_status']['finished']
    assert data_dict['run_status']['ionic_converged']
    assert data_dict['run_status']['electronic_converged']
    assert data_dict['run_status']['nelm'] == 60
    assert data_dict['run_status']['nsw'] == 61
    assert data_dict['fermi_level'] == pytest.approx(6.17267267)
    assert np.amax(np.linalg.norm(data_dict['stress'], axis=1)) == pytest.approx(42.96872956444064)
    assert np.amax(np.abs(np.linalg.norm(data_dict['forces'], axis=1))) == pytest.approx(0.21326679)
    assert data_dict['total_energies']['energy_extrapolated'] == pytest.approx(-10.823296)


def test_vasprun_parsing(parser_with_vasprun):
    parser, _ = parser_with_vasprun(
        {'parser_settings': {'check_completeness': False, 'required_quantity': [], 'critical_objects': []}}
    )
    misc = parser.outputs['misc'].get_dict()
    shoud_exists = ['fermi_level', 'forces', 'stress', 'band_properties']
    for name in shoud_exists:
        assert name in misc


def test_outcar_parsing(parser_with_outcar):
    parser, _ = parser_with_outcar(
        {'parser_settings': {'check_completeness': False, 'required_quantity': [], 'critical_objects': []}}
    )
    misc = parser.outputs['misc'].get_dict()
    shoud_exists = ['run_stats', 'run_status']
    for name in shoud_exists:
        assert name in misc


def test_basic_run(parser_with_retrieved):
    parser, _ = parser_with_retrieved(
        'basic_run', {'parser_settings': {'include_quantity': ['born_charges'], 'critical_objects': []}}
    )

    misc = parser.outputs['misc'].get_dict()
    assert misc['band_properties']['cbm'] == pytest.approx(5.075)
    assert misc['band_properties']['vbm'] == pytest.approx(4.2811)
    assert misc['band_properties']['band_gap'] == pytest.approx(0.793899999999999)
    assert not misc['band_properties']['is_direct_gap']
    assert misc['version'] == '5.3.5'
    assert misc['total_energies']['energy_extrapolated'] == pytest.approx(-36.09616894)
    assert np.amax(np.linalg.norm(misc['stress'], axis=1)) == pytest.approx(8.50955439)
    assert np.amax(np.linalg.norm(misc['forces'], axis=1)) == pytest.approx(0.0)
    assert misc['run_status']['nelm'] == 60
    assert misc['run_status']['last_iteration_index'] == [1, 2]
    assert misc['run_status']['nsw'] == 0
    assert misc['run_status']['finished']
    assert misc['run_status']['ionic_converged'] is None
    assert misc['run_status']['electronic_converged']
    assert not misc['run_status']['consistent_nelm_breach']
    assert not misc['run_status']['contains_nelm_breach']
    assert misc['run_stats']['mem_usage_base'] == pytest.approx(30000.0)
    assert misc['run_stats']['mem_usage_nonl-proj'] == pytest.approx(7219.0)
    assert misc['run_stats']['mem_usage_fftplans'] == pytest.approx(776.0)
    assert misc['run_stats']['mem_usage_grid'] == pytest.approx(1605.0)
    assert misc['run_stats']['mem_usage_one-center'] == pytest.approx(124.0)
    assert misc['run_stats']['mem_usage_wavefun'] == pytest.approx(814.0)
    assert misc['run_stats']['maximum_memory_used'] == pytest.approx(95344.0)
    assert misc['run_stats']['average_memory_used'] == pytest.approx(0.0)
    assert misc['run_stats']['total_cpu_time_used'] == pytest.approx(20.463)
    assert misc['run_stats']['user_time'] == pytest.approx(11.266)
    assert misc['run_stats']['system_time'] == pytest.approx(9.197)
    assert misc['run_stats']['elapsed_time'] == pytest.approx(22.518)
    assert misc.get('magnetization') is None
    assert misc.get('site_magnetization') is None


def test_relax_run(parser_with_retrieved):
    parser, _ = parser_with_retrieved(
        'relax',
        {
            'parser_settings': {
                'include_quantity': ['born_charges'],
                'include_node': ['energies', 'trajectory'],
                'check_completeness': False,
                'required_quantity': [],
                'electronic_step_energies': True,
                'critical_objects': [],
            }
        },
    )

    array = parser.outputs['energies']
    energies_ext = array.get_array('energy_extrapolated')
    energies_ext_elec = array.get_array('energy_extrapolated_electronic')
    energies_elec_steps = array.get_array('electronic_steps')
    test_array_energies = [
        np.array(
            [
                163.37398579,
                14.26925896,
                -23.05190509,
                -34.91615104,
                -40.20080347,
                -42.18390876,
                -42.97469852,
                -43.31556073,
                -43.60169068,
                -43.61723125,
                -43.61871511,
                -43.61879751,
                -43.12548175,
                -42.90647187,
                -42.91031846,
                -42.91099027,
                -42.91111107,
                -42.91113348,
            ]
        ),
        np.array([-43.34236449, -43.31102002, -43.27768275, -43.27791002, -43.27761357, -43.27757545]),
        np.array([-43.40320524, -43.38084022, -43.36835045, -43.36666248, -43.36666583, -43.36649036, -43.36648855]),
        np.array([-43.37749056, -43.37749102, -43.37734414, -43.37734069]),
        np.array([-43.38117265, -43.38082881, -43.38063293, -43.38062479]),
        np.array([-43.38337336, -43.38334165]),
        np.array([-43.38778922, -43.38766017, -43.38752953, -43.38753003]),
        np.array([-43.38714489, -43.38708193]),
        np.array([-43.38640951, -43.38641449]),
        np.array([-43.3874799, -43.3871553, -43.38701949, -43.38701639]),
        np.array([-43.38790942, -43.38727062, -43.38700335, -43.38699488]),
        np.array([-43.38774394, -43.38773717]),
        np.array([-43.38984942, -43.3899134, -43.38988315]),
        np.array([-43.38988117, -43.3898822]),
        np.array([-43.39032165, -43.39017866, -43.39011239]),
        np.array([-43.39021044, -43.39020751]),
        np.array([-43.39034135, -43.39034244]),
        np.array([-43.39044466, -43.39044584]),
        np.array([-43.39084354, -43.39088709, -43.39087657]),
    ]
    test_array_steps = np.array([18, 6, 7, 4, 4, 2, 4, 2, 2, 4, 4, 2, 3, 2, 3, 2, 2, 2, 3])
    # Build a flattened array (not using flatten from NumPy as the content is staggered) and
    # test number of electronic steps per ionic step
    test_array_energies_flattened = np.array([])
    for ionic_step in test_array_energies:
        test_array_energies_flattened = np.append(test_array_energies_flattened, ionic_step)
    assert energies_ext_elec.shape == test_array_energies_flattened.shape
    np.testing.assert_allclose(test_array_energies_flattened, energies_ext_elec, atol=0.0, rtol=1.0e-7)
    np.testing.assert_allclose(test_array_steps, energies_elec_steps, atol=0.0, rtol=1.0e-7)
    test_array_energies = np.array(
        [
            -0.00236637,
            -0.00048614,
            -0.00047201,
            -0.00043261,
            -0.00041668,
            -0.00042584,
            -0.00043637,
            -0.00042806,
            -0.00042762,
            -0.00043875,
            -0.00042731,
            -0.00042705,
            -0.00043064,
            -0.00043051,
            -0.00043161,
            -0.00043078,
            -0.00043053,
            -0.00043149,
            -0.00043417,
        ]
    )
    # Testing on VASP 5, where the extrapolated energy should be the following due to a bug
    with np.testing.assert_raises(AssertionError):
        np.testing.assert_allclose(test_array_energies, energies_ext, atol=0.0, rtol=1.0e-7)
    # Instead we correct and it should be
    test_array = np.array(
        [
            -42.911133,
            -43.277575,
            -43.366489,
            -43.377341,
            -43.380625,
            -43.383342,
            -43.38753,
            -43.387082,
            -43.386414,
            -43.387016,
            -43.386995,
            -43.387737,
            -43.389883,
            -43.389882,
            -43.390112,
            -43.390208,
            -43.390342,
            -43.390446,
            -43.390877,
        ]
    )
    np.testing.assert_allclose(test_array, energies_ext, atol=0.0, rtol=1.0e-7)

    traj = parser.outputs['trajectory']
    assert 'energy_extrapolated' in traj.get_arraynames()
    assert 'forces' in traj.get_arraynames()
    assert 'stress' in traj.get_arraynames()
    assert traj.get_array('forces').shape == (19, 8, 3)
    traj.get_stepids()
    traj.get_cells()
    traj.get_positions()
    traj.get_step_data(1)
    assert isinstance(traj.get_step_structure(1), orm.StructureData)


def test_basic(parser_with_retrieved):
    parser, _ = parser_with_retrieved(
        'basic',
        {
            'parser_settings': {
                'critical_objects': [],
                'check_completeness': False,
                'include_quantity': ['kpoints'],
                'required_quantity': [],
                'include_node': ['kpoints'],
                'check_errors': False,
            }
        },
    )
    assert 'kpoints' in parser.outputs
    kpoints = parser.outputs['kpoints']
    np.testing.assert_allclose(kpoints.get_kpoints()[0], np.array([0.0, 0.0, 0.0]), atol=0.0, rtol=1.0e-7)
    np.testing.assert_allclose(
        kpoints.get_kpoints()[-1], np.array([0.42857143, -0.42857143, 0.42857143]), atol=0.0, rtol=1.0e-7
    )


def test_stream(parser_with_retrieved):
    """Test the functionality of the stream parser."""
    parser, _ = parser_with_retrieved(
        'stdout/out',
        {
            'parser_settings': {
                'critical_objects': [],
                'check_completeness': False,
                'include_quantity': ['stream'],
                'required_quantity': [],
            }
        },
    )
    misc_dict = parser.outputs['misc'].get_dict()
    assert misc_dict['notifications'][0]['name'] == 'ibzkpt'
    assert misc_dict['notifications'][0]['kind'] == 'ERROR'
    assert misc_dict['notifications'][0]['regex'] == 'internal error in subroutine IBZKPT'
    assert misc_dict['notifications'][1]['name'] == 'nostart'
    assert misc_dict['notifications'][1]['kind'] == 'ERROR'
    assert misc_dict['notifications'][1]['regex'] == 'vasp.'


def test_parser_exception(request, calc_with_retrieved):
    """
    This calculation has a missing eigenvalues section in the vasprun.xml
    """
    # This should work as the parser does not output the band by default
    # But the diagonsis information is missing so the erorr code is not zero
    settings_dict = {
        'parser_settings': {'check_completeness': False, 'critical_objects': ['vasprun.xml', 'vasp_output', 'OUTCAR']}
    }
    file_path = str(pathlib.Path(request.fspath).parent / '../test_data/basic_run_ill_format')
    node = calc_with_retrieved(file_path, settings_dict)
    _, output = VaspParser.parse_from_node(node, store_provenance=False, retrieved_temporary_folder=file_path)

    assert output.is_finished
    assert output.exit_status == 0

    # 1004 - node cannot be created as the quantity is missing
    settings_dict = {'parser_settings': {'check_completeness': True, 'include_node': ['bands']}}
    file_path = str(pathlib.Path(request.fspath).parent / '../test_data/basic_run_ill_format')
    node = calc_with_retrieved(file_path, settings_dict)
    _, output = VaspParser.parse_from_node(node, store_provenance=False, retrieved_temporary_folder=file_path)

    assert output.is_finished
    assert output.exit_status == 1004

    # 1002 - we explicitly require the eigenvalues to be present, but it is not
    settings_dict = {
        'parser_settings': {
            'critical_objects': [],
            'check_completeness': True,
            'required_quantity': ['eigenvalues'],
            'include_node': ['band'],
        }
    }
    file_path = str(pathlib.Path(request.fspath).parent / '../test_data/basic_run_ill_format')
    node = calc_with_retrieved(file_path, settings_dict)
    _, output = VaspParser.parse_from_node(node, store_provenance=False, retrieved_temporary_folder=file_path)

    assert output.is_finished
    assert output.exit_status == 1002


def test_notification_composer(parser_with_retrieved):
    """Test the NotificationComposer class"""

    parser, exit_code = parser_with_retrieved(
        'basic',
        {
            'parser_settings': {
                'critical_objects': [],
                'check_completeness': False,
                'include_quantity': ['kpoints'],
                'required_quantity': [],
            }
        },
        parse=False,
    )

    notifications = [{'name': 'edwav', 'kind': 'ERROR', 'message': 'Error in EDWAV'}]
    config = ParserSettingsConfig()
    composer = NotificationComposer(
        notifications,
        {},
        {
            'parameters': orm.Dict(dict={'nelect': 10}),
        },
        parser.exit_codes,
        critical_notifications=config.critical_notification_errors,
    )
    exit_code = composer.compose()
    assert exit_code.status == 703

    # BRMIX error but has NELECT defined in the input
    notifications = [{'name': 'brmix', 'kind': 'ERROR', 'message': 'Error in BRMIX'}]
    composer = NotificationComposer(
        notifications,
        {},
        {
            'parameters': orm.Dict(dict={'nelect': 10}),
        },
        parser.exit_codes,
        critical_notifications=config.critical_notification_errors,
    )
    exit_code = composer.compose()
    assert exit_code is None

    # BRMIX error but no NELECT tag
    composer = NotificationComposer(
        notifications,
        {},
        {'parameters': orm.Dict(dict={})},
        parser.exit_codes,
        critical_notifications=config.critical_notification_errors,
    )
    exit_code = composer.compose()
    assert exit_code.status == 703
