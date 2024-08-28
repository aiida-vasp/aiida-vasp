"""
Tools for using the Pymatgen library with aiida-vasp.
"""

import shutil
import tempfile
from contextlib import contextmanager
from pathlib import Path
from typing import List, Optional

import pymatgen.io.vasp as pvasp

from aiida_vasp.utils.export import export_vasp

from .aiida_utils import ensure_node_first_arg, ensure_node_kwargs
from .export import export_vasp_calc


@contextmanager
def temporary_folder():
    """Get a temporary folder and delete it after use."""
    tmpf = Path(tempfile.mkdtemp())
    yield tmpf
    shutil.rmtree(tmpf)


class PymatgenAdapator:
    """Adaptor for getting pymatgen objects from a VASP calculation/workflow"""

    FILES = {
        'potcar': ('Potcar', 'POTCAR'),
        'vasprun': ('Vasprun', 'vasprun.xml'),
        'kpoints': ('Kpoints', 'KPOINTS'),
        'ibzkpt': ('Kpoints', 'IBZKPT'),
        'incar': ('Incar', 'INCAR'),
        'outcar': ('Outcar', 'OUTCAR'),
        'contcar': ('Poscar', 'CONTCAR'),
        'poscar': ('Poscar', 'POSCAR'),
        'chgcar': ('Chgcar', 'CHGCAR'),
    }
    # Classes where from_dict is not implemented but still MSONable
    NO_FROM_DICT = ['vasprun', 'outcar', 'chgcar']

    def __init__(self, node):
        """Adaptor for getting pymatgen objects from a VASP calculation/workflow"""
        self.node = node
        self.pmg_objects = {}
        self.cache = {}

    def _parse_full(self, names: Optional[List[str]] = None):
        """
        Parse all files and save to the pmg_objects attribute
        The assumption is that exporting the calculation folder is the slowest part of the process.
        """
        if names is None:
            names = self.FILES
        else:
            names = {key: self.FILES[key] for key in names}

        with temporary_folder() as tmpf:
            export_vasp(self.node, tmpf)
            for name, (cls_name, file) in names.items():
                # Instantiate the pymatgen object
                cls = self.pmg_objects[name] = getattr(pvasp, cls_name)
                if not Path(tmpf / file).is_file():
                    continue
                if hasattr(cls, 'from_file'):
                    obj = cls.from_file(str(tmpf / file))
                else:
                    obj = cls(str(tmpf / file))
                self.pmg_objects[name] = obj

    def _get_pmg_object(self, name):
        """
        Get a pymatgen object

        1. If we can find the object in the cache, then just return it.
        """
        # We already parsed the calculation, so we can just return the object
        # Since we have access it - save it to the cache
        if name in self.pmg_objects:
            if name not in self.cache:
                self.cache[name + '_dict'] = self.pmg_objects[name].as_dict()
            return self.pmg_objects[name]

        # Check if the object was accessed before and previously save to the cache
        if not self.cache:
            self.cache = self.node.base.extras.get('pmg_cache', {})

        if name + '_dict' in self.cache and name not in self.NO_FROM_DICT:
            # Already in the cache - return the object
            return getattr(pvasp, self.FILES[name][0]).from_dict(self.cache[name])
        else:
            # Not in the cache
            if name not in self.pmg_objects:
                # We have not paresed the calculation yet - do it now
                self._parse_full()
            # Get the object
            obj = self.pmg_objects[name]
            # Save the parsed object to the cache
            self.cache[name + '_dict'] = obj.as_dict()
            return obj

    def _get_pmg_dict(self, name):
        """
        Get a pymatgen object as a dictionary
        """
        if not self.cache:
            self.cache = self.node.base.extras.get('pmg_cache', {})
        if name + '_dict' in self.cache:
            return self.cache[name + '_dict']
        return self._get_pmg_object(name).as_dict()

    def _flush_cache(self):
        """Close the adaptor and save the cache"""
        self.node.base.extras.set('pmg_cache', self.cache)

    def __enter__(self):
        """Enter the adaptor"""
        return self

    def __exit__(self, *args, **kwargs):
        self._flush_cache()

    @property
    def vasprun(self):
        """Return the Vasprun object"""
        return self._get_pmg_object('vasprun')

    @property
    def vasprun_dict(self):
        """Return the Vasprun object as dictionary (will trigger caching)"""
        return self._get_pmg_dict('vasprun')

    @property
    def outcar(self):
        """Return the Outcar object"""
        return self._get_pmg_object('outcar')

    @property
    def outcar_dict(self):
        """Return the Outcar object as dictionary (will trigger caching)"""
        return self._get_pmg_dict('outcar')

    @property
    def poscar(self):
        """Return the Poscar object"""
        return self._get_pmg_object('poscar')

    @property
    def poscar_dict(self):
        """Return the Poscar object as dictionary (will trigger caching)"""
        return self._get_pmg_dict('poscar')

    @property
    def incar(self):
        """Return the Incar object"""
        return self._get_pmg_object('incar')

    @property
    def incar_dict(self):
        """Return the Incar object as dictionary"""
        return self._get_pmg_dict('incar')

    @property
    def kpoints(self):
        """Return the Kpoints object"""
        return self._get_pmg_object('kpoints')

    @property
    def kpoints_dict(self):
        """Return the Kpoints object as dictionary"""
        return self._get_pmg_dict('kpoints')

    @property
    def ibzkpt(self):
        """Return the IBZKPT object"""
        return self._get_pmg_object('ibzkpt')

    @property
    def ibzkpt_dict(self):
        """Return the IBZKPT object as dictionary"""
        return self._get_pmg_dict('ibzkpt')

    def save_msonable(self, name, obj):
        """Save msonable object to the node extras"""
        dobj = obj.as_dict()
        assert '@module' in dobj
        self.node.base.extras.set(f'pmg_cache_{name}', dobj)

    def load_msonable(self, name):
        """Load msonable object from the node extras"""
        from monty.json import MontyDecoder

        return MontyDecoder().process_decoded(self.node.base.extras.get(f'pmg_cache_{name}'))


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

    with temporary_folder() as tmpf:
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

    return vrun, outcar
