"""
Tools for using the Pymatgen library with aiida-vasp.
"""

from __future__ import annotations

import shutil
import tempfile
from contextlib import contextmanager
from pathlib import Path
from typing import Any, Dict, Generator, List, Optional

try:
    import pymatgen.io.vasp as pvasp
except ImportError:
    raise ImportError('You need to install pymatgen to use this module.')


import gzip

from aiida import orm

from .aiida_utils import ensure_node_first_arg
from .export import export_vasp


@contextmanager
def temporary_folder() -> Generator[Path, None, None]:
    """Get a temporary folder and delete it after use."""
    tmpf = Path(tempfile.mkdtemp())
    yield tmpf
    shutil.rmtree(tmpf)


class PymatgenAdapator:
    """
    Adaptor for getting pymatgen objects from a VASP calculation/workflow
    This work by first exporting the calculation to a temporary folder and then parsing the files using pymmatgen.

    Some of the pymatgen objects does not have the from_dict method implemented as required by MSONable.
    Hence, they can only be reconstructed as a dictionary.
    """

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
    NO_RECONSTRUCT = ['vasprun', 'outcar', 'chgcar']

    def __init__(self, node: orm.CalcJobNode, store_cache: bool = True) -> None:
        """Adaptor for getting pymatgen objects from a VASP calculation/workflow"""
        self.node = node
        self.pmg_objects = {}
        self.cache = {}
        self.store_cache = store_cache

    def _parse_full(self, file_names: Optional[List[str]] = None) -> None:
        """
        Parse all files and save to the pmg_objects attribute
        The assumption is that exporting the calculation folder is the slowest part of the process.
        """
        if file_names is None:
            file_names = self.FILES
        else:
            file_names = {key: self.FILES[key] for key in file_names}

        with temporary_folder() as tmpf:
            export_vasp(self.node, tmpf)
            for name, (cls_name, file) in file_names.items():
                # Instantiate the pymatgen object
                cls = getattr(pvasp, cls_name)
                if not Path(tmpf / file).is_file():
                    continue
                # Try use the from_file method if it exists
                fname = str(tmpf / file)
                if hasattr(cls, 'from_file'):
                    try:
                        obj = cls.from_file(fname)
                    # Skip if the file is not found or the parsing fails
                    except Exception:
                        continue
                # If using Vasprun, try to parse the potcar file but fall back when needed
                elif cls == pvasp.Vasprun:
                    try:
                        obj = cls(fname)
                    except ValueError:
                        obj = cls(fname, parse_potcar=False)
                else:
                    obj = cls(fname)
                self.pmg_objects[name] = obj

    def export_files(self, dst: str | Path) -> None:
        """Export the VASP calculation files to a destination folder"""
        export_vasp(self.node, dst)

    def _get_pmg_object(self, name: str) -> Any:
        """
        Get a pymatgen object

        1. If we can find the object in parsed object , then just return it.
        2. If it is not already parsed, try to load the cache (stored in the extras)
        3. Otherwise, try to export and parse from the files explicitly. (slow)

        :param name: Name of the object to get (e.g. 'vasprun', 'outcar', 'poscar', 'incar', 'kpoints', 'ibzkpt')
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

        if name + '_dict' in self.cache and name not in self.NO_RECONSTRUCT:
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

    def _get_pmg_dict(self, name: str) -> Dict[str, Any]:
        """
        Get a pymatgen object as a dictionary
        """
        if not self.cache:
            self.cache = self.node.base.extras.get('pmg_cache', {})
        if name + '_dict' in self.cache:
            return self.cache[name + '_dict']
        return self._get_pmg_object(name).as_dict()

    def _flush_cache(self) -> None:
        """Close the adaptor and save the cache"""
        self.node.base.extras.set('pmg_cache', self.cache)

    def __enter__(self) -> 'PymatgenAdapator':
        """Enter the adaptor"""
        return self

    def __exit__(self, *args, **kwargs) -> None:
        if self.store_cache:
            self._flush_cache()

    @property
    def vasprun(self) -> pvasp.Vasprun:
        """Return the Vasprun object"""
        return self._get_pmg_object('vasprun')

    @property
    def vasprun_dict(self) -> Dict[str, Any]:
        """Return the Vasprun object as dictionary (will trigger caching)"""
        return self._get_pmg_dict('vasprun')

    @property
    def outcar(self) -> pvasp.Outcar:
        """Return the Outcar object"""
        return self._get_pmg_object('outcar')

    @property
    def outcar_dict(self) -> Dict[str, Any]:
        """Return the Outcar object as dictionary (will trigger caching)"""
        return self._get_pmg_dict('outcar')

    @property
    def poscar(self) -> pvasp.Poscar:
        """Return the Poscar object"""
        return self._get_pmg_object('poscar')

    @property
    def poscar_dict(self) -> Dict[str, Any]:
        """Return the Poscar object as dictionary (will trigger caching)"""
        return self._get_pmg_dict('poscar')

    @property
    def incar(self) -> pvasp.Incar:
        """Return the Incar object"""
        return self._get_pmg_object('incar')

    @property
    def incar_dict(self) -> Dict[str, Any]:
        """Return the Incar object as dictionary"""
        return self._get_pmg_dict('incar')

    @property
    def kpoints(self) -> pvasp.Kpoints:
        """Return the Kpoints object"""
        return self._get_pmg_object('kpoints')

    @property
    def kpoints_dict(self) -> Dict[str, Any]:
        """Return the Kpoints object as dictionary"""
        return self._get_pmg_dict('kpoints')

    @property
    def ibzkpt(self) -> pvasp.Kpoints:
        """Return the IBZKPT object"""
        return self._get_pmg_object('ibzkpt')

    @property
    def ibzkpt_dict(self) -> Dict[str, Any]:
        """Return the IBZKPT object as dictionary"""
        return self._get_pmg_dict('ibzkpt')

    def save_msonable(self, name: str, obj: Any) -> None:
        """Save msonable object to the node extras"""
        dobj = obj.as_dict()
        assert '@module' in dobj
        self.node.base.extras.set(f'pmg_cache_{name}', dobj)

    def load_msonable(self, name: str) -> Any:
        """Load msonable object from the node extras"""
        from monty.json import MontyDecoder  # noqa: PLC0415

        return MontyDecoder().process_decoded(self.node.base.extras.get(f'pmg_cache_{name}'))


@ensure_node_first_arg
def get_vasprun(node: Any, store_cache: bool = True) -> pvasp.Vasprun:
    """Return the Vasprun object"""
    return PymatgenAdapator(node, store_cache=store_cache).vasprun


@ensure_node_first_arg
def get_outcar(node: Any, store_cache: bool = True) -> pvasp.Outcar:
    """Return the OUTCAR object"""
    return PymatgenAdapator(node, store_cache=store_cache).outcar


@ensure_node_first_arg
def get_incar(node: Any, store_cache: bool = True) -> pvasp.Incar:
    """Return the INCAR object"""
    return PymatgenAdapator(node, store_cache=store_cache).incar


@ensure_node_first_arg
def get_kpoints(node: Any, store_cache: bool = True) -> pvasp.Kpoints:
    """Return the Kpoints object"""
    return PymatgenAdapator(node, store_cache=store_cache).kpoints


@ensure_node_first_arg
def get_ibzkpt(node: Any, store_cache: bool = True) -> pvasp.Kpoints:
    """Return the Kpoints object using the IBZKPT file"""
    return PymatgenAdapator(node, store_cache=store_cache).ibzkpt


def convert_pymatgen_potcar_folder(src: Path | str, dst: Path | str) -> None:
    """
    Convert pymatgen potcar folder to a structure used by aiida-vasp

    :param src: Path to the pymatgen potcar folder
    :param dst: Path to the aiida-vasp potcar folder

    :return: None
    """

    src = Path(src)
    dst = Path(dst)
    for fpath in Path(src).glob('POTCAR.*.gz'):
        symbol = fpath.name.split('.')[1]
        folder = dst / symbol
        folder.mkdir(exist_ok=True, parents=True)
        # unzip the file
        with gzip.open(fpath, 'rb') as f_in:
            with (folder / 'POTCAR').open('wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
