"""
Mock vasp code.

A more advanced way of mocking. The input of a calculation can be
hash and match from a repository of calculation that has been run.

This way we can perform tests for workchain without the need for
injecting test code into the workchain logic itself.
"""
import hashlib
import shutil
import os
from typing import Union
from pathlib import Path
import logging

import numpy as np

from aiida.repository import FileType

from parsevasp.kpoints import Kpoints
from parsevasp.incar import Incar
from parsevasp.poscar import Poscar

from aiida_vasp.utils.fixtures.testdata import data_path

# pylint: disable=logging-format-interpolation, import-outside-toplevel

INPUT_OBJECTS = ('POSCAR', 'INCAR', 'KPOINTS')
EXCLUDED = ('POTCAR', '.aiida')


def get_hash(dict_obj):
    """
    Return the hash for a dictionary of arbitrary items.

    This is not meant to be robust for all cases, but should be OK for
    matching hashes of calculation inputs.

    The point here is to make the value invariant to the permutation of key orders.
    """

    # If a list is passed - convert it to a dictionary with keys being the indices
    if isinstance(dict_obj, list):
        dict_obj = dict(enumerate(dict_obj))

    rec = []
    for key, value in dict_obj.items():
        key = repr(key)
        # For numpy/list with floating point zero (0.0) we have to converge -0.0 to 0.0
        # as they should be equivalent
        if isinstance(value, np.ndarray):
            value[value == 0] = 0
        elif isinstance(value, list) and 0 in value:
            value = [type(tmp)(0) if tmp == 0 else tmp for tmp in value]

        # Handle if value itself is float zero
        if isinstance(value, float) and value == 0:
            value = 0.

        if isinstance(value, (dict, list)):
            rec.append(key + ':' + get_hash(value)[0])
        else:
            # Use the string representation
            rec.append(key + ':' + repr(value) + ':' + repr(type(value)))

    # Update, use sorted so the original order does not matter, in force case so
    # sting keys with upper/lower cases are treated as the same
    base = [record.encode().lower() for record in sorted(rec)]
    # Compute the hash
    md5 = hashlib.md5()
    for item in base:
        md5.update(item)

    return md5.hexdigest(), base


class MockRegistry:
    """
    A class to create and manage a registry of completed calculations.

    Calculations are identified using the hash of the parsed inputs.
    """

    def __init__(self, base_path=data_path('.')):
        """
        Instantiate and Registry
        """
        self.base_path = Path(base_path)
        self.reg_hash = {}
        self.reg_name = {}
        self.logger = logging.getLogger('aiida_vasp.utils.mock_code.MockRegistry')
        self._setup_logger()
        self.scan()

    def _setup_logger(self, level=logging.INFO):
        """Setup the logger"""
        formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
        handler = logging.StreamHandler()
        handler.setLevel(level)
        handler.setFormatter(formatter)
        self.logger.setLevel(level)

    def scan(self):
        """
        Scan the base folder and locate input/output folders
        """
        for output_folder in Path(self.base_path).glob('**/out'):
            calc_base_folder = output_folder.parent.absolute()
            self._register_folder(calc_base_folder)

    def get_path_by_hash(self, hash_val):
        """
        Return the output folder for a given hash
        """
        return Path(self.reg_hash[hash_val])

    def get_path_by_name(self, name):
        """
        Return the output folder for a given hash
        """
        return Path(self.reg_hash[self.reg_name[name]])

    @staticmethod
    def compute_hash(input_folder: Path):
        """
        Compute the hash of a input folder
        """
        items = {}
        kpt_path = input_folder / 'KPOINTS'
        if kpt_path.is_file():
            kpoints = Kpoints(file_path=str(kpt_path))
            items['kpoints'] = kpoints.get_dict()
            items['kpoints'].pop('comment', None)

        incar_path = input_folder / 'INCAR'
        if incar_path.is_file():
            incar = Incar(file_path=str(incar_path), validate_tags=False)
            items['incar'] = incar.get_dict()

        poscar_path = input_folder / 'POSCAR'
        if not poscar_path.is_file():
            poscar_path = input_folder / '00/POSCAR'
        if poscar_path.is_file():
            poscar = Poscar(file_path=str(poscar_path))
            items['poscar'] = poscar.get_dict()
            items['poscar'].pop('comment', None)

        return get_hash(items)[0]

    def extract_calc_by_path(self, rel_path: Path, dst_path: Path, include_inputs: bool = True):
        """
        Copy the content of a give hash to a destination.

        :param rel_path: The relative path of the calculation folder to be
          extracted.
        :param dst: The destination path to be extracted to - must already exists.
        """
        rel_path = Path(rel_path)
        dst_path = Path(dst_path)

        base_out = self.base_path / rel_path / 'out'
        base_in = self.base_path / rel_path / 'inp'

        if not base_out.exists() or not base_in.exists():
            raise ValueError(f'Relative path: {rel_path} is invalid')

        # Copy the content of input and then the output folder
        paths = [base_in, base_out] if include_inputs else [base_out]
        for folder in paths:
            for fpath in folder.glob('*'):
                if fpath.is_file():
                    shutil.copy2(fpath, dst_path)
                elif fpath.is_dir():
                    shutil.copytree(fpath, dst_path / fpath.name, dirs_exist_ok=True)

    def extract_calc_by_hash(self, hash_val, dst, include_inputs=False):
        """
        Extract an registerred calculation using hash.
        """
        self.extract_calc_by_path(self.get_path_by_hash(hash_val), dst, include_inputs)

    def upload_calc(self, folder: Path, rel_path: Union[Path, str], excluded_object=None):
        """
        Register a calculation folder to the repository
        """
        inp = list(INPUT_OBJECTS)
        excluded = list(EXCLUDED)
        if excluded_object:
            excluded.extend(excluded_object)

        # Check if the repository folder already exists
        repo_calc_base = self.base_path / rel_path
        if repo_calc_base.exists():
            raise FileExistsError(f'There is already a direcotry at {repo_calc_base.resolve()}.')

        # Deposit the objects
        repo_calc_base.mkdir(parents=True)
        repo_in = repo_calc_base / 'inp'
        repo_out = repo_calc_base / 'out'
        repo_in.mkdir(parents=True)
        repo_out.mkdir(parents=True)

        for obj in folder.glob('*'):
            if obj.name in inp:
                shutil.copy2(obj, repo_in)
            elif obj.name not in excluded:
                if obj.is_file():
                    shutil.copy2(obj, repo_out)
                elif obj.is_dir():
                    shutil.copytree(obj, repo_out / obj.name)

        # Update the hash table
        self._register_folder(repo_calc_base)

    def _register_folder(self, calc_base: Path):
        """
        Register a folder inside the repository
        """
        # Get the relative path to the base
        rel = calc_base.relative_to(self.base_path)
        # Compute the hash
        hash_val = self.compute_hash(calc_base / 'inp')
        # Link absolute path to hash, and hash to relative path (used as name)
        self.reg_hash[hash_val] = calc_base.absolute()
        self.reg_name[str(rel)] = hash_val

    def upload_aiida_calc(self, calc_node, rel_path: Union[str, Path], excluded_names=None):
        """
        Register an aiida calc_class
        """
        # Check if the repository folder already exists
        repo_calc_base = self.base_path / rel_path
        if repo_calc_base.exists():
            raise FileExistsError(f'There is already a directory at {repo_calc_base.resolve()}.')

        # Deposit the objects
        repo_calc_base.mkdir(parents=True)
        repo_in = repo_calc_base / 'inp'
        repo_out = repo_calc_base / 'out'
        repo_in.mkdir(parents=True)
        repo_out.mkdir(parents=True)

        exclude = list(EXCLUDED)
        if excluded_names:
            exclude.extend(excluded_names)

        # Copy the input objects
        for obj in calc_node.list_objects():
            if obj.name in exclude:
                continue
            copy_from_aiida(obj.name, calc_node, repo_in)

        # Copy the retrieved objects
        for obj in calc_node.outputs.retrieved.list_objects():
            if obj.name in exclude:
                continue
            copy_from_aiida(obj.name, calc_node.outputs.retrieved, repo_out)

        self.logger.info('Calculation %s has been registered', calc_node)
        self._register_folder(repo_calc_base)

    def upload_aiida_work(self, worknode, rel_path: Union[str, Path]):
        """
        Upload all calculations in a workchain node
        """
        from aiida.orm import CalcJobNode
        from aiida.plugins import CalculationFactory
        calc_class = CalculationFactory('vasp.vasp')
        neb_class = CalculationFactory('vasp.neb')
        to_upload = []
        for node in worknode.called_descendants:
            if isinstance(node, CalcJobNode) and (node.process_class in [calc_class, neb_class]):
                to_upload.append(node)
        to_upload.sort(key=lambda x: x.ctime)
        self.logger.info('Collected %s nodes to upload under name %s.', to_upload, rel_path)

        for idx, node in enumerate(to_upload):
            rel = Path(rel_path) / f'calc-{idx:03d}'
            self.upload_aiida_calc(node, rel)
        self.logger.info('WorkChain %s has been uploaded.', worknode)


class MockVasp:
    """
    Mock VaspExecutable
    """

    def __init__(self, workdir: Union[str, Path], registry: MockRegistry):
        """
        Mock VASP executable that copies over outputs from existing calculations.
        Inputs are hash and looked for.
        """
        self.workdir = workdir
        self.registry = registry

    def run(self, debug=True):
        """
        Run the mock vasp
        """
        hash_val = self.registry.compute_hash(self.workdir)
        with open('/tmp/testing', 'a') as handler:
            handler.write(str(self.workdir))
        if debug:
            print(f'Target hash value: {hash_val}')
        if hash_val in self.registry.reg_hash:
            self.registry.extract_calc_by_hash(hash_val, self.workdir)
        else:
            if debug:
                print(f'Registered hashes: {self.registry.reg_hash}')
            raise ValueError('The calculation is not registered!!')

    @property
    def is_runnable(self) -> bool:
        """Return wether the mock code can be run"""
        hash_val = self.registry.compute_hash(self.workdir)
        return hash_val in self.registry.reg_hash


def copy_from_aiida(name: str, node, dst: Path):
    """
    Copy objects from aiida repository.

    :param name: The full name (including the parent path) of the object.
    :param node: Node object for which the objects in the repo to be copied.
    :param dst: Path of the destination folder.

    This is a recursive function so directory copying also works.
    """
    obj = node.get_object(name)

    # If it is a directory, copy the contents one by one
    if obj.file_type == FileType.DIRECTORY:
        for sub_obj in node.list_objects(name):
            copy_from_aiida(os.path.join(name, sub_obj.name), node, dst)
    else:
        # Anything else
        with node.open(name) as fsource:
            # Make parent directory if needed
            frepo_path = dst / name
            Path(frepo_path.parent).mkdir(exist_ok=True, parents=True)
            # Write the object
            with open(frepo_path, 'w') as fdst:
                shutil.copyfileobj(fsource, fdst)
