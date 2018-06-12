"""Module defining sets of FileParsers to be used by the VaspParser"""

from aiida_vasp.io.doscar import DosParser
from aiida_vasp.io.eigenval import EigParser
from aiida_vasp.io.kpoints import KpParser
from aiida_vasp.io.outcar import OutcarParser
from aiida_vasp.io.vasprun import VasprunParser
from aiida_vasp.io.chgcar import ChgcarParser
from aiida_vasp.io.wavecar import WavecarParser
from aiida_vasp.io.poscar import PoscarParser

DEFAULT_PARSABLE_ITEMS = {
    'structure': {
        'inputs': [],
        'parsers': ['vasprun.xml', 'CONTCAR'],
        'nodeName': 'structure',
        'prerequisites': []
    },
    'bands': {
        'inputs': [],
        'parsers': ['vasprun.xml', 'EIGENVAL'],
        'nodeName': 'bands',
        'prerequisites': []
    },
    'dos': {
        'inputs': [],
        'parsers': ['vasprun.xml', 'DOSCAR'],
        'nodeName': 'dos',
        'prerequisites': []
    },
    'kpoints': {
        'inputs': [],
        'parsers': ['vasprun.xml', 'IBZKPT'],
        'nodeName': 'kpoints',
        'prerequisites': []
    },
    'occupations': {
        'inputs': [],
        'parsers': ['vasprun.xml'],
        'nodeName': 'occupations',
        'prerequisites': []
    },
    'trajectory': {
        'inputs': [],
        'parsers': ['vasprun.xml'],
        'nodeName': 'trajectory',
        'prerequisites': []
    },
    'energies': {
        'inputs': [],
        'parsers': ['vasprun.xml', 'OUTCAR'],
        'nodeName': 'energies',
        'prerequisites': []
    },
    'projectors': {
        'inputs': [],
        'parsers': ['vasprun.xml'],
        'nodeName': 'projectors',
        'prerequisites': []
    },
    'dielectrics': {
        'inputs': [],
        'parsers': ['vasprun.xml'],
        'nodeName': 'dielectrics',
        'prerequisites': []
    },
    'final_stress': {
        'inputs': [],
        'parsers': ['vasprun.xml'],
        'nodeName': '',
        'prerequisites': []
    },
    'final_forces': {
        'inputs': [],
        'parsers': ['vasprun.xml'],
        'nodeName': '',
        'prerequisites': []
    },
    'final_structure': {
        'inputs': [],
        'parsers': ['vasprun.xml'],
        'nodeName': '',
        'prerequisites': []
    },
    'born_charges': {
        'inputs': [],
        'parsers': ['vasprun.xml'],
        'nodeName': 'born_charges',
        'prerequisites': []
    },
    'hessian': {
        'inputs': [],
        'parsers': ['vasprun.xml'],
        'nodeName': 'hessian',
        'prerequisites': []
    },
    'dynmat': {
        'inputs': [],
        'parsers': ['vasprun.xml'],
        'nodeName': 'dynmat',
        'prerequisites': []
    },
    'parameters': {
        'inputs': [],
        'parsers': ['vasprun.xml', 'OUTCAR'],
        'nodeName': 'parameters',
        'prerequisites': []
    }
    # eFL: Aiida can already calculate volume, right?
    # no reason to parse (or store)
    
    #'volume': {
    #    'inputs': ['parameters'],
    #    'parsers': ['vasprun.xml', 'OUTCAR'],
    #    'nodeName': 'parameters',
    #    'prerequisites': []
    #}

    # eFL: symmetries are in fact not presently available
    # in vasprun.xml, but will be included in the future
    #        'symmetries': {
    #            'inputs': [],
    #            'parsers': ['OUTCAR'],
    #            'nodeName': 'parameters',
    #            'prerequisites': []
    #        }
}


FILE_PARSER_SETS = {
    'default': {
        'DOSCAR': {
            'parser_class': DosParser,
            'is_critical': False,
            'status': 'Unknown'
        },
        'EIGENVAL': {
            'parser_class': EigParser,
            'is_critical': False,
            'status': 'Unknown'
        },
        'IBZKPT': {
            'parser_class': KpParser,
            'is_critical': False,
            'status': 'Unknown'
        },
        'OUTCAR': {
            'parser_class': OutcarParser,
            'is_critical': False,
            'status': 'Unknown'
        },
        'vasprun.xml': {
            'parser_class': VasprunParser,
            'is_critical': True,
            'status': 'Unknown'
        },
        'CHGCAR': {
            'parser_class': ChgcarParser,
            'is_critical': False,
            'status': 'Unknown'
        },
        'WAVECAR': {
            'parser_class': WavecarParser,
            'is_critical': False,
            'status': 'Unknown'
        },
        'CONTCAR': {
            'parser_class': PoscarParser,
            'is_critical': False,
            'status': 'Unknown'
        },
    },
}


def get_file_parser_set(parser_set='default'):
    if parser_set not in FILE_PARSER_SETS:
        return None
    return FILE_PARSER_SETS[parser_set]
