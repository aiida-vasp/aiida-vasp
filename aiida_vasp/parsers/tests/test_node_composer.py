@pytest.mark.parametrize(['poscar_parser'], [('poscar',)], indirect=True)
@pytest.mark.parametrize(['vasp_structure'], [('str-Al',)], indirect=True)
def test_parse_poscar_data(fresh_aiida_env, vasp_structure, poscar_parser):
    """
    Parse a reference POSCAR.

    Check that the node composer works and returns the expected content.

    """

    # Compose the node using the defined content parser
    inputs = get_node_composer_inputs_from_object_parser(poscar_parser, quantity_keys=['poscar-structure'])
    result = NodeComposer.compose('structure', inputs)

    # Compare
    assert result.cell == vasp_structure.cell
    assert result.get_site_kindnames() == vasp_structure.get_site_kindnames()
    assert result.sites[2].position == vasp_structure.sites[2].position


@pytest.mark.parametrize(['vasp_structure'], [('str-Al',)], indirect=True)
def test_parse_poscar_write(fresh_aiida_env, vasp_structure, tmpdir):
    """
    Parse (write) a reference POSCAR.

    Using the PoscarParser, write, read, compose node and compare to reference
    structure.

    """

    # Initialize content parser
    parser = PoscarParser(data=vasp_structure)

    # Write the loaded structure to file
    temp_path = str(tmpdir.join('POSCAR'))
    parser.write(temp_path)

    # Load the written structure using a new content parser instance
    parser = None
    with open(temp_path, 'r') as handler:
        parser = PoscarParser(handler=handler)

    # Compose the node
    inputs = get_node_composer_inputs_from_object_parser(parser, quantity_keys=['poscar-structure'])
    result = NodeComposer.compose('structure', inputs)

    # Then compare
    assert result.cell == vasp_structure.cell
    assert result.get_site_kindnames() == \
        vasp_structure.get_site_kindnames()
    assert result.sites[2].position == \
        vasp_structure.sites[2].position


@pytest.mark.xfail(aiida_version() < cmp_version('1.0.0a1'), reason='Element X only present in Aiida >= 1.x')
def test_parse_poscar_silly_read(fresh_aiida_env):
    """
    Parse (read) a reference POSCAR with silly elemental names.

    Using the PoscarParser and compare the result to a reference
    structure.

    """

    # Setup parser
    path = data_path('poscar', 'POSCARSILLY')
    parser = None
    with open(path, 'r') as handler:
        parser = PoscarParser(handler=handler)

    # Compose the node
    inputs = get_node_composer_inputs_from_object_parser(parser, quantity_keys=['poscar-structure'])
    result = NodeComposer.compose('structure', inputs)
    names = result.get_site_kindnames()
    assert names == ['Hamburger', 'Pizza']
    symbols = result.get_symbols_set()
    assert symbols == set(['X', 'X'])


@pytest.mark.parametrize(['vasp_structure'], [('str-InAs',)], indirect=True)
def test_parse_poscar_silly_write(fresh_aiida_env, vasp_structure, tmpdir):
    """
    Parse (read, write) a reference POSCAR with silly elemental names.

    Using the PoscarParser and compare the result to a reference structure.

    """

    # Initialize the content parser
    parser = PoscarParser(data=vasp_structure)
    # Compose the node
    inputs = get_node_composer_inputs_from_object_parser(parser, quantity_keys=['poscar-structure'])
    print(inputs)
    result = NodeComposer.compose('structure', inputs)

    # Compare
    names = result.get_site_kindnames()
    assert names == ['Hamburger', 'Pizza']
    symbols = result.get_symbols_set()
    assert symbols == set(['As', 'In'])

    temp_path = str(tmpdir.join('POSCAR'))
    parser.write(temp_path)

    parser = None
    with open(temp_path, 'r') as handler:
        parser = PoscarParser(handler=handler)

    # Compose a new node
    inputs = get_node_composer_inputs_from_object_parser(parser, quantity_keys=['poscar-structure'])
    result = NodeComposer.compose('structure', inputs)

    # Compare
    names = result.get_site_kindnames()
    assert names == ['Hamburger', 'Pizza']
    symbols = result.get_symbols_set()
    assert symbols == set(['X', 'X'])


@pytest.mark.parametrize(['vasp_structure'], [('str',)], indirect=True)
def test_parse_poscar_undercase(fresh_aiida_env, vasp_structure, tmpdir):
    """
    Parse a reference POSCAR.

    With potential elemental names using the PoscarParser and compare
    the result to a reference structure.

    """

    parser = PoscarParser(data=vasp_structure)
    result = parser.get_quantity('poscar-structure')
    names = result.get_site_kindnames()
    assert names == ['In', 'As', 'As', 'In_d', 'In_d', 'As']
    symbols = result.get_symbols_set()
    assert symbols == set(['As', 'In'])
    temp_path = str(tmpdir.join('POSCAR'))
    parser.write(temp_path)
    parser = None
    with open(temp_path, 'r') as handler:
        parser = PoscarParser(handler=handler)
    result_reparse = parser.structure
    names = result_reparse.get_site_kindnames()
    assert names == ['In', 'As', 'As', 'In_d', 'In_d', 'As']
    symbols = result_reparse.get_symbols_set()
    assert symbols == set(['As', 'In'])


def test_parse_kpoints(vasp_kpoints):
    """
    Parse a reference KPOINTS.

    Using the KpointsParser and compare the result to a reference
    kpoints-node.

    """

    kpoints, _ = vasp_kpoints

    try:
        _ = kpoints.get_attribute('mesh')
        path = data_path('kpoints', 'KPOINTS_mesh')
        method = 'get_kpoints_mesh'
        param = 'mesh'
    except AttributeError:
        pass

    try:
        _ = kpoints.get_attribute('array|kpoints')
        path = data_path('kpoints', 'KPOINTS_list')
        method = 'get_kpoints'
        param = 'list'
    except AttributeError:
        pass

    parser = None
    with open(path, 'r') as handler:
        parser = KpointsParser(file_handler=handler)
    result = parser.kpoints
    if param == 'list':
        assert getattr(result, method)().all() == getattr(kpoints, method)().all()
    if param == 'mesh':
        assert getattr(result, method)() == getattr(kpoints, method)()


@pytest.mark.parametrize('outcar_parser', ['disp_details'], indirect=True)
def test_parameter_results(fresh_aiida_env, outcar_parser):
    """
    Test that the parameter node is a ParametersData instance.

    Should contain the symmetries and the elastic moduli.

    """

    outcar_parser._settings._output_nodes_dict.update(  # pylint: disable=protected-access
        {'misc': {
            'type': 'dict',
            'quantities': ['symmetries_extended', 'elastic_moduli', 'run_stats', 'run_status'],
            'link_name': 'my_custom_node'
        }})

    inputs = get_node_composer_inputs_from_object_parser(outcar_parser,
                                                         quantity_keys=['symmetries_extended', 'elastic_moduli', 'run_stats', 'run_status'])
    data_obj = NodeComposer.compose('dict', inputs)
    ref_class = get_data_class('dict')
    assert isinstance(data_obj, ref_class)
    data_dict = data_obj.get_dict()
    # test symmetries
    test = {
        'symmetrized_cell_type': {
            'static': [
                'face centered cubic supercell.', 'body centered tetragonal supercell.', 'body centered tetragonal supercell.',
                'body centered tetragonal supercell.', 'body centered tetragonal supercell.', 'body centered tetragonal supercell.',
                'body centered tetragonal supercell.', 'base centered monoclinic supercell.', 'base centered monoclinic supercell.',
                'base centered monoclinic supercell.', 'base centered monoclinic supercell.', 'base centered monoclinic supercell.',
                'base centered monoclinic supercell.', 'face centered cubic supercell.', 'face centered cubic supercell.',
                'face centered cubic supercell.'
            ],
            'dynamic': [
                'face centered cubic supercell.', 'body centered tetragonal supercell.', 'body centered tetragonal supercell.',
                'body centered tetragonal supercell.', 'body centered tetragonal supercell.', 'body centered tetragonal supercell.',
                'body centered tetragonal supercell.', 'base centered monoclinic supercell.', 'base centered monoclinic supercell.',
                'base centered monoclinic supercell.', 'base centered monoclinic supercell.', 'base centered monoclinic supercell.',
                'base centered monoclinic supercell.', 'face centered cubic supercell.', 'face centered cubic supercell.',
                'face centered cubic supercell.'
            ]
        },
        'point_group': {
            'static': [
                'O_h', 'D_4h.', 'D_4h.', 'D_4h.', 'D_4h.', 'D_4h.', 'D_4h.', 'C_2h.', 'C_2h.', 'C_2h.', 'C_2h.', 'C_2h.', 'C_2h.', 'D_2h.',
                'D_2h.', 'O_h'
            ],
            'dynamic': [
                'O_h', 'D_4h.', 'D_4h.', 'D_4h.', 'D_4h.', 'D_4h.', 'D_4h.', 'C_2h.', 'C_2h.', 'C_2h.', 'C_2h.', 'C_2h.', 'C_2h.', 'D_2h.',
                'D_2h.', 'O_h'
            ]
        },
        'original_cell_type': {
            'static': [
                'primitive cell', 'primitive cell', 'primitive cell', 'primitive cell', 'primitive cell', 'primitive cell',
                'primitive cell', 'primitive cell', 'primitive cell', 'primitive cell', 'primitive cell', 'primitive cell',
                'primitive cell', 'primitive cell', 'primitive cell', 'primitive cell'
            ],
            'dynamic': [
                'primitive cell', 'primitive cell', 'primitive cell', 'primitive cell', 'primitive cell', 'primitive cell',
                'primitive cell', 'primitive cell', 'primitive cell', 'primitive cell', 'primitive cell', 'primitive cell',
                'primitive cell', 'primitive cell', 'primitive cell', 'primitive cell'
            ]
        },
        'num_space_group_operations': {
            'static': [48, 16, 16, 16, 16, 16, 16, 4, 4, 4, 4, 4, 4, 8, 8, 48],
            'dynamic': [48, 16, 16, 16, 16, 16, 16, 4, 4, 4, 4, 4, 4, 8, 8, 48]
        },
        'primitive_translations': [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
    }
    assert set(data_dict['symmetries']) == set(test)

    # then elastic moduli
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


_TEST_DATA = [
    (['outcar_extras', 'OUTCAR.converged'], [True, True, True, False, False]),
    (['outcar_extras', 'OUTCAR.nelm-breach-consistent'], [True, False, False, True, True]),
    (['outcar_extras', 'OUTCAR.nelm-breach-partial'], [True, False, True, False, True]),
    (['outcar_extras', 'OUTCAR.unfinished'], [False, False, False, False, False]),
    (['outcar_extras', 'OUTCAR.not-converged'], [True, False, True, False, False]),
]


@pytest.mark.parametrize('outcar_parser,expected', _TEST_DATA, indirect=['outcar_parser'])
def test_run_status(fresh_aiida_env, outcar_parser, expected):
    """
    Test the run_status obtained by checking the convergence problems of a calculation,
    finished or unfinished.
    """

    outcar_parser._settings._output_nodes_dict.update(  # pylint: disable=protected-access
        {'misc': {
            'type': 'dict',
            'quantities': ['run_stats', 'run_status'],
            'link_name': 'my_custom_node'
        }})

    inputs = get_node_composer_inputs_from_object_parser(outcar_parser, quantity_keys=['run_stats', 'run_status'])
    data_obj = NodeComposer.compose('dict', inputs)
    ref_class = get_data_class('dict')
    assert isinstance(data_obj, ref_class)
    data_dict = data_obj.get_dict()

    assert data_dict['run_status']['finished'] is expected[0]
    assert data_dict['run_status']['ionic_converged'] is expected[1]
    assert data_dict['run_status']['electronic_converged'] is expected[2]
    assert data_dict['run_status']['consistent_nelm_breach'] is expected[3]
    assert data_dict['run_status']['contains_nelm_breach'] is expected[4]
