"""Module defining sets of FileParsers to be used by the VaspParser"""

from aiida_vasp.io.doscar import DosParser
from aiida_vasp.io.eigenval import EigParser
from aiida_vasp.io.kpoints import KpParser
from aiida_vasp.io.outcar import OutcarParser
from aiida_vasp.io.vasprun import VasprunParser
from aiida_vasp.io.chgcar import ChgcarParser
from aiida_vasp.io.wavecar import WavecarParser
from aiida_vasp.io.poscar import PoscarParser

DEFAULT_PRIORITY = 100

FILE_PARSER_SETS = {
    'default': {
        'DOSCAR': {
            'parser_class': DosParser,
            'is_critical': False,
            'status': 'Unknown',
            'priority': DEFAULT_PRIORITY
        },
        'EIGENVAL': {
            'parser_class': EigParser,
            'is_critical': False,
            'status': 'Unknown',
            'priority': DEFAULT_PRIORITY
        },
        'IBZKPT': {
            'parser_class': KpParser,
            'is_critical': False,
            'status': 'Unknown',
            'priority': DEFAULT_PRIORITY
        },
        'OUTCAR': {
            'parser_class': OutcarParser,
            'is_critical': False,
            'status': 'Unknown',
            'priority': DEFAULT_PRIORITY
        },
        'vasprun.xml': {
            'parser_class': VasprunParser,
            'is_critical': True,
            'status': 'Unknown',
            'priority': DEFAULT_PRIORITY + 1
        },
        'CHGCAR': {
            'parser_class': ChgcarParser,
            'is_critical': False,
            'status': 'Unknown',
            'priority': DEFAULT_PRIORITY
        },
        'WAVECAR': {
            'parser_class': WavecarParser,
            'is_critical': False,
            'status': 'Unknown',
            'priority': DEFAULT_PRIORITY
        },
        'CONTCAR': {
            'parser_class': PoscarParser,
            'is_critical': False,
            'status': 'Unknown',
            'priority': DEFAULT_PRIORITY
        },
    },
}


def get_file_parser_set(parser_set='default'):
    if parser_set not in FILE_PARSER_SETS:
        return None
    return FILE_PARSER_SETS[parser_set]
