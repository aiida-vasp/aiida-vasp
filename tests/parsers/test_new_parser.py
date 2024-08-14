import pathlib


def test_parser_bare(calc_with_retrieved, request):
    from aiida_vasp.parsers.parser_new import VaspParser

    _relative_file_path = '../test_data/basic_run'
    file_path = str(pathlib.Path(request.fspath).parent / _relative_file_path)
    node = calc_with_retrieved(file_path, {})
    parser = VaspParser(node)
    parser.parse(retrieved_tempoary_folder=file_path)
    assert 'misc' in parser.outputs
    assert 'structure' in parser.outputs

    node = calc_with_retrieved(file_path, {'parser_settings': {'include_quantity': ['projectors']}})
    parser = VaspParser(node)
    parser.parse(retrieved_tempoary_folder=file_path)
    assert 'arrays' in parser.outputs
