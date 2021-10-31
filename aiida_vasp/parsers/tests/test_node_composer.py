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
