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

from aiida.repository import FileType

from parsevasp.kpoints import Kpoints
from parsevasp.incar import Incar
from parsevasp.poscar import Poscar

from aiida_vasp.utils.fixtures.testdata import data_path


def get_hash(dict_obj):
    """
    Return the hash for a dictionary of arbitrary items.

    This is not meant to be robust for all cases, but should be OK for
    matching hashes of calculation inputs.

    The point here is to make the value invariant to the permutation of key orders.
    """

    rec = []
    for key, value in dict_obj.items():
        key = str(key)
        if isinstance(value, dict):
            rec.append(key + ':' + get_hash(value)[0])
        else:
            rec.append(key + ':' + str(value) + ':' + str(type(value)))
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
        self.base_path = base_path
        self.reg_hash = {}
        self.reg_name = {}

        self.scan()

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
        kpt_file = input_folder / 'KPOINTS'
        if kpt_file.is_file():
            kpoints = Kpoints(file_path=str(kpt_file))
            items['kpoints'] = kpoints.get_dict()

        incar_file = input_folder / 'INCAR'
        if incar_file.is_file():
            incar = Incar(file_path=str(incar_file))
            items['incar'] = incar.get_dict()

        poscar_file = input_folder / 'POSCAR'
        if poscar_file.is_file():
            poscar = Poscar(file_path=str(poscar_file))
            items['poscar'] = poscar.get_dict()

        return get_hash(items)[0]

    def extract_calc_by_path(self, rel_path: Path, dst_path: Path):
        """
        Copy the content of a give hash to a destination.

        :param rel_path: The relative path of the calcualtion folder to be
          extracted.
        :param dst: The destination path to be extracted to - must already exists.
        """
        rel_path = Path(rel_path)
        dst_path = Path(dst_path)

        base_out = self.base_path / rel_path / 'out'
        base_in = self.base_path / rel_path / 'inp'

        # Copy the content of input and then the output folder
        for folder in [base_in, base_out]:
            for fpath in folder.glob('*'):
                if fpath.is_file():
                    shutil.copy2(fpath, dst_path)
                elif fpath.is_dir():
                    shutil.copytree(fpath, dst_path / fpath.name)

    def extract_calc_by_hash(self, hash_val, dst):
        """
        Extract an registerred calculation using hash.
        """
        self.extract_calc_by_path(self.get_path_by_hash(hash_val), dst)

    def upload_calc(self, folder: Path, rel_path: Union[Path, str], excluded_file=None):
        """
        Register a calculation folder to the repository
        """
        inp = ['POSCAR', 'INCAR', 'KPOINTS']
        excluded = ['POTCAR']
        if excluded_file:
            excluded.extend(excluded_file)

        # Check if the repository folder already exists
        repo_calc_base = self.base_path / rel_path
        if repo_calc_base.exists():
            raise FileExistsError(f'There is already a direcotry at {repo_calc_base.resolve()}.')

        # Deposit the files
        repo_calc_base.mkdir()
        repo_in = repo_calc_base / 'inp'
        repo_out = repo_calc_base / 'out'
        repo_in.mkdir()
        repo_out.mkdir()

        for file in folder.glob('*'):
            if file.name in inp:
                shutil.copy2(file, repo_in)
            elif file.name not in excluded:
                if file.is_file():
                    shutil.copy2(file, repo_out)
                elif file.is_dir():
                    shutil.copytree(file, repo_out / file.name)

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

    def upload_aiida_calc(self, calc_node, rel_path, excluded_names=None):
        """
        Register an aiida VaspCalculation
        """
        # Check if the repository folder already exists
        repo_calc_base = self.base_path / rel_path
        if repo_calc_base.exists():
            raise FileExistsError(f'There is already a directory at {repo_calc_base.resolve()}.')

        # Deposit the files
        repo_calc_base.mkdir()
        repo_in = repo_calc_base / 'inp'
        repo_out = repo_calc_base / 'out'
        repo_in.mkdir()
        repo_out.mkdir()

        exclude = ['.aiida']
        if excluded_names:
            exclude.extend(excluded_names)

        # Copy input files
        for obj in calc_node.list_objects():
            if obj.name in exclude:
                continue
            copy_from_aiida(obj.name, calc_node, repo_in)

        # Copy retrieved files
        for obj in calc_node.outputs.retrieved.list_objects():
            if obj.name in exclude:
                continue
            copy_from_aiida(obj.name, calc_node, repo_out)

        print(f'Calculation {calc_node} has been registered')
        self._register_folder(repo_calc_base)


class MockVasp:
    """
    Mock VaspExecutable
    """

    def __init__(self, workdir: Union[str, Path], repository: MockRegistry):
        """
        Mock VASP executable that copies over outputs from existing calculations.
        Inputs are hash and looked for.
        """
        self.workdir = workdir
        self.repository = repository

    def run(self):
        """
        Run the mock vasp
        """
        hash_val = self.repository.compute_hash(self.workdir)
        if self.repository.has_hash(hash_val):
            self.repository.extract_calc_by_hash(hash_val)
        else:
            raise ValueError('The calculation is not registered!!')


def copy_from_aiida(name: str, node, dst: Path):
    """
    Copy objects from aiida repository.

    :param name: The full name (including the parent path) of the object.
    :param node: Node object for which the files in the repo to be copied.
    :param dst: Path of the destination folder.

    This is a recursive function so directory copying also works.
    """
    obj = node.get_object(name)

    # If it is a directory, copy the contents one by one
    if obj.file_type == FileType.DIRECTORY:
        for sub_obj in node.list_objects(name):
            copy_from_aiida(os.path.join(name, sub_obj.name), node, dst)
    else:
        # It is a file
        with node.open(name) as fsource:
            # Make parent directory if needed
            frepo_path = dst / name
            Path(frepo_path.parent).mkdir(exist_ok=True, parents=True)
            # Write the file
            with open(frepo_path, 'w') as fdst:
                shutil.copyfileobj(fsource, fdst)
