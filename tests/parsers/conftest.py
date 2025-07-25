import pytest
from aiida import orm
from aiida.common.links import LinkType
from aiida.orm import CalcJobNode
from aiida.plugins import DataFactory

from aiida_vasp.parsers.content_parsers.chgcar import ChgcarParser
from aiida_vasp.parsers.content_parsers.doscar import DoscarParser
from aiida_vasp.parsers.content_parsers.eigenval import EigenvalParser
from aiida_vasp.parsers.content_parsers.incar import IncarParser
from aiida_vasp.parsers.content_parsers.kpoints import KpointsParser
from aiida_vasp.parsers.content_parsers.outcar import OutcarParser, VtstNebOutcarParser
from aiida_vasp.parsers.content_parsers.poscar import PoscarParser
from aiida_vasp.parsers.content_parsers.stream import StreamParser
from aiida_vasp.parsers.content_parsers.vasprun import VasprunParser


@pytest.fixture()
def calc_with_retrieved(localhost):
    """A rigged CalcJobNode for testing the parser and that the calculation retrieve what is expected."""

    def _inner(file_path, input_settings=None):
        # Create a test computer
        computer = localhost

        process_type = 'aiida.calculations:vasp.vasp'

        node = CalcJobNode(computer=computer, process_type=process_type)
        node.base.attributes.set('input_filename', 'INCAR')
        node.base.attributes.set('output_filename', 'OUTCAR')
        # node.set_attribute('error_filename', 'aiida.err')
        node.base.attributes.set('scheduler_stderr', '_scheduler-stderr.txt')
        node.base.attributes.set('scheduler_stdout', '_scheduler-stdout.txt')
        node.set_option('resources', {'num_machines': 1, 'num_mpiprocs_per_machine': 1})
        node.set_option('max_wallclock_seconds', 1800)

        if input_settings is None:
            input_settings = {}

        settings = DataFactory('core.dict')(dict=input_settings)
        node.base.links.add_incoming(settings, link_type=LinkType.INPUT_CALC, link_label='settings')
        settings.store()
        node.store()

        # Create a `FolderData` that will represent the `retrieved` folder. Store the test
        # output fixture in there and link it.
        retrieved = DataFactory('core.folder')()
        retrieved.base.repository.put_object_from_tree(file_path)
        retrieved.base.links.add_incoming(node, link_type=LinkType.CREATE, link_label='retrieved')
        retrieved.store()

        return node

    return _inner


@pytest.fixture()
def vasp_structure_poscar(vasp_structure):
    """Fixture: Well formed POSCAR contents."""
    aiida_structure = vasp_structure
    if isinstance(vasp_structure, orm.CifData):
        ase_structure = vasp_structure.get_ase()
        aiida_structure = orm.StructureData(ase=ase_structure)
    writer = PoscarParser(data=aiida_structure)
    return writer


@pytest.fixture()
def vasprun_parser(request, data_path):
    """Return an instance of VasprunParser for a reference vasprun.xml. Remember rb mode."""
    path, settings = path_file_and_settings('vasprun.xml', request.param, data_path)
    with open(path, 'rb') as handler:
        parser = VasprunParser(handler=handler, settings=settings)
    return parser


@pytest.fixture()
def vasprun_parser_v621(request, data_path):
    """Return an instance of VasprunParser for a reference vasprun.xml of VASP6."""
    path, settings = path_file_and_settings('vasprun621.xml', request.param, data_path)
    with open(path, 'r', encoding='utf8') as handler:
        parser = VasprunParser(handler=handler, settings=settings)
    return parser


@pytest.fixture()
def neb_outcar_parser(request, data_path):
    """Return an instance of OutcarParser for a reference OUTCAR."""
    path, settings = path_file_and_settings('OUTCAR', request.param, data_path)
    with open(path, 'r', encoding='utf8') as handler:
        parser = VtstNebOutcarParser(handler=handler, settings=settings)
    return parser


@pytest.fixture()
def outcar_parser(request, data_path):
    """Return an instance of OutcarParser for a reference OUTCAR."""
    path, settings = path_file_and_settings('OUTCAR', request.param, data_path)
    with open(path, 'r', encoding='utf8') as handler:
        parser = OutcarParser(handler=handler, settings=settings)
    return parser


@pytest.fixture()
def poscar_parser(request, data_path):
    """Return an instance of PoscarParser for a reference POSCAR."""
    path, _ = path_file_and_settings('POSCAR', request.param, data_path)
    with open(path, 'r', encoding='utf8') as handler:
        parser = PoscarParser(handler=handler)
    return parser


@pytest.fixture()
def incar_parser(request, data_path):
    """Return an instance of IncarParser for a reference INCAR."""
    path, settings = path_file_and_settings('INCAR', request.param, data_path)
    with open(path, 'r', encoding='utf8') as handler:
        parser = IncarParser(handler=handler, settings=settings)
    return parser


@pytest.fixture()
def doscar_parser(request, data_path):
    """Return an instance of DoscarParser for a reference DOSCAR."""
    path, settings = path_file_and_settings('DOSCAR', request.param, data_path)
    with open(path, 'r', encoding='utf8') as handler:
        parser = DoscarParser(handler=handler, settings=settings)
    return parser


@pytest.fixture()
def chgcar_parser(request, data_path):
    """Return an instance of ChgcarParser for a reference CHGCAR."""
    path, settings = path_file_and_settings('CHGCAR', request.param, data_path)
    with open(path, 'r', encoding='utf8') as handler:
        parser = ChgcarParser(handler=handler, settings=settings)
    return parser


@pytest.fixture()
def eigenval_parser(request, data_path):
    """Return an instance of EigenvalParser for a reference EIGENVAL."""
    path, settings = path_file_and_settings('EIGENVAL', request.param, data_path)
    with open(path, 'r', encoding='utf8') as handler:
        parser = EigenvalParser(handler=handler, settings=settings)
    return parser


@pytest.fixture()
def kpoints_parser(request, data_path):
    """Return an instance of KpointsParser for a reference KPOINTS."""
    path, settings = path_file_and_settings('KPOINTS', request.param, data_path)
    with open(path, 'r', encoding='utf8') as handler:
        parser = KpointsParser(handler=handler, settings=settings)
    return parser


@pytest.fixture()
def stream_parser(request, data_path):
    """Return an instance of StreamParser for a reference stream capture."""
    path, settings = path_file_and_settings('vasp_output', request.param, data_path)
    with open(path, 'r', encoding='utf8') as handler:
        parser = StreamParser(handler=handler, settings=settings)
    return parser


def path_file_and_settings(name, param, data_path):
    """Locate folder, filename and settings from param. Return the path and settings."""
    settings = {}
    if isinstance(param, list):
        if len(param) == 3:
            folder, name, settings = param
        elif len(param) == 2:
            folder, name = param
        else:
            raise IndexError(
                'Please supply either folder and name, or folder, name and settings to the parser fixtures'
            )
    else:
        folder = param
    path = data_path(folder, name)

    return path, settings


@pytest.fixture
def compare_symmetries():
    return {
        'symmetrized_cell_type': {
            'static': [
                'face centered cubic supercell.',
                'body centered tetragonal supercell.',
                'body centered tetragonal supercell.',
                'body centered tetragonal supercell.',
                'body centered tetragonal supercell.',
                'body centered tetragonal supercell.',
                'body centered tetragonal supercell.',
                'base centered monoclinic supercell.',
                'base centered monoclinic supercell.',
                'base centered monoclinic supercell.',
                'base centered monoclinic supercell.',
                'base centered monoclinic supercell.',
                'base centered monoclinic supercell.',
                'face centered cubic supercell.',
                'face centered cubic supercell.',
                'face centered cubic supercell.',
            ],
            'dynamic': [
                'face centered cubic supercell.',
                'body centered tetragonal supercell.',
                'body centered tetragonal supercell.',
                'body centered tetragonal supercell.',
                'body centered tetragonal supercell.',
                'body centered tetragonal supercell.',
                'body centered tetragonal supercell.',
                'base centered monoclinic supercell.',
                'base centered monoclinic supercell.',
                'base centered monoclinic supercell.',
                'base centered monoclinic supercell.',
                'base centered monoclinic supercell.',
                'base centered monoclinic supercell.',
                'face centered cubic supercell.',
                'face centered cubic supercell.',
                'face centered cubic supercell.',
            ],
        },
        'original_cell_type': {
            'static': [
                'primitive cell',
                'primitive cell',
                'primitive cell',
                'primitive cell',
                'primitive cell',
                'primitive cell',
                'primitive cell',
                'primitive cell',
                'primitive cell',
                'primitive cell',
                'primitive cell',
                'primitive cell',
                'primitive cell',
                'primitive cell',
                'primitive cell',
                'primitive cell',
            ],
            'dynamic': [
                'primitive cell',
                'primitive cell',
                'primitive cell',
                'primitive cell',
                'primitive cell',
                'primitive cell',
                'primitive cell',
                'primitive cell',
                'primitive cell',
                'primitive cell',
                'primitive cell',
                'primitive cell',
                'primitive cell',
                'primitive cell',
                'primitive cell',
                'primitive cell',
            ],
        },
        'num_space_group_operations': {
            'static': [48, 16, 16, 16, 16, 16, 16, 4, 4, 4, 4, 4, 4, 8, 8, 48],
            'dynamic': [48, 16, 16, 16, 16, 16, 16, 4, 4, 4, 4, 4, 4, 8, 8, 48],
        },
    }


@pytest.fixture()
def neb_calc_with_retrieved(localhost):
    """A rigged CalcJobNode for testing the parser and that the calculation retrieve what is expected."""

    def _inner(file_path, input_settings=None, nimgs=3):
        # Create a test computer
        computer = localhost

        process_type = 'aiida.calculations:vasp.vasp'

        node = CalcJobNode(computer=computer, process_type=process_type)
        node.base.attributes.set('input_filename', 'INCAR')
        node.base.attributes.set('output_filename', 'OUTCAR')
        # node.set_attribute('error_filename', 'aiida.err')
        node.base.attributes.set('scheduler_stderr', '_scheduler-stderr.txt')
        node.base.attributes.set('scheduler_stdout', '_scheduler-stdout.txt')
        node.set_option('resources', {'num_machines': 1, 'num_mpiprocs_per_machine': 1})
        node.set_option('max_wallclock_seconds', 1800)

        if input_settings is None:
            input_settings = {}

        settings = DataFactory('core.dict')(dict=input_settings)
        node.base.links.add_incoming(settings, link_type=LinkType.INPUT_CALC, link_label='settings')
        settings.store()

        # Add inputs with the number of images
        param = orm.Dict(dict={'images': nimgs})
        node.base.links.add_incoming(param, link_type=LinkType.INPUT_CALC, link_label='parameters')
        param.store()

        node.store()

        # Create a `FolderData` that will represent the `retrieved` folder. Store the test
        # output fixture in there and link it.
        retrieved = DataFactory('core.folder')()
        retrieved.base.repository.put_object_from_tree(file_path)
        retrieved.base.links.add_incoming(node, link_type=LinkType.CREATE, link_label='retrieved')
        retrieved.store()

        return node

    return _inner
