import pathlib

import pytest
from aiida_vasp.parsers.parser_new import VaspParser


def test_parser_bare(calc_with_retrieved, request):
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


@pytest.fixture
def parser_with_retrieved(calc_with_retrieved, request):
    def wrapped(name, settings={}):
        _relative_file_path = f'../test_data/{name}'
        file_path = str(pathlib.Path(request.fspath).parent / _relative_file_path)
        node = calc_with_retrieved(file_path, settings)
        parser = VaspParser(node)
        parser.parse(retrieved_tempoary_folder=file_path)
        return parser

    return wrapped


def test_parser_born(parser_with_retrieved):
    parser = parser_with_retrieved('born_effective_charge', {'parser_settings': {'include_quantity': ['born_charges']}})
    assert 'born_charges' in parser.outputs['arrays'].get_arraynames()
