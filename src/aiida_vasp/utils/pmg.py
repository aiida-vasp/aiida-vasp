"""
Tools for using the Pymatgen library with aiida-vasp.
"""

import shutil
import tempfile
from pathlib import Path

from .aiida_utils import ensure_node_first_arg, ensure_node_kwargs
from .export import export_vasp_calc


@ensure_node_first_arg
@ensure_node_kwargs
def get_pymatgen_objects(node, parse_xml=True, parse_potcar_file=False, parse_outcar=True, **kwargs):
    """
    Return pymatgen objects for a VaspCalculation or a VaspWorkChain

    When given a VaspWorkChain, the final, completed calculation is used.
    """
    try:
        from pymatgen.io.vasp import Outcar, Vasprun
    except ImportError:
        raise ImportError('You need to install pymatgen to use this feature.')

    tmpf = Path(tempfile.mkdtemp())
    export_vasp_calc(node, tmpf)

    def gz_if_necessary(fname):
        if not fname.is_file():
            fname = str(fname) + '.gz'
        else:
            fname = str(fname)
        return fname

    if parse_xml:
        vrun = Vasprun(
            gz_if_necessary(tmpf / 'vasprun.xml'),
            parse_potcar_file=parse_potcar_file,
            **kwargs,
        )
    else:
        vrun = None
    if parse_outcar:
        outcar = Outcar(gz_if_necessary(tmpf / 'OUTCAR'))
    else:
        outcar = None

    # Clean up the temporary directory
    shutil.rmtree(tmpf)

    return vrun, outcar
