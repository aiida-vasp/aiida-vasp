"""
Mock vasp code.

A more advanced way of mocking. The input of a calculation can be
hash and match from a repository of calculation that has been run.

This way we can perform tests for workchain without the need for
injecting test code into the workchain logic itself.
"""

from __future__ import annotations

import hashlib
import logging
import os
import pathlib
import shutil
import warnings
from subprocess import run

import numpy as np
from aiida import orm
from aiida.plugins import CalculationFactory
from aiida.repository import FileType
from parsevasp.incar import Incar
from parsevasp.kpoints import Kpoints
from parsevasp.poscar import Poscar

# pylint: disable=logging-format-interpolation, import-outside-toplevel

INPUT_OBJECTS = ('POSCAR', 'INCAR', 'KPOINTS')
# Objects (files) that should be excluded from storing in the repository
# These objects are mostly not used for parsing and can be large in size
DEFAULT_EXCLUDED = (
    'POTCAR',
    '.aiida',
    'LOCPOT',
    'DOSCAR',
    'PCDAT',
    'EIGENVAL',
    'OSZICAR',
    'PCDAT',
    'XDATCAR',
    'REPORT',
    'WAVECAR',
    'CHG',
)


def data_path(*args: str) -> str:
    """Return a path to a file in the test data directory."""
    path = pathlib.Path(__file__).parent.parent.parent.parent / 'tests' / 'test_data' / pathlib.Path(*args)
    path = path.resolve()
    if not path.exists():
        warnings.warn(f'Path {path} for test data does not exist.')
    return str(path)


def get_hash(dict_obj: dict | list) -> tuple[str, list[bytes]]:
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
    for key_, value_ in dict_obj.items():
        key = repr(key_)
        value = value_
        # For numpy/list with floating point zero (0.0) we have to converge -0.0 to 0.0
        # as they should be equivalent
        if isinstance(value, np.ndarray):
            value[value == 0] = 0

        # Handle if value itself is float zero
        if isinstance(value, float) and value == 0:
            value = 0.0

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

    The class uses environmental variables to control its behaviour:

    - MOCK_{CODE}_REG_BASE: Prefix to the upload relative path
    - MOCK_{CODE}_UPLOAD_PREFIX: Prefix to the upload relative path

    The `{CODE}` is replaced with the `CODE_NAME` class attribute of the subclass.

    """

    CODE_NAME = 'ABSTRACT'

    def __init__(self, base_path: str | pathlib.Path | list[str | pathlib.Path] | None = None) -> None:
        """
        Instantiate and Registry
        """
        if base_path is None:
            base_path = os.environ.get(f'MOCK_{self.CODE_NAME}_REG_BASE', data_path('.'))
        base_path = pathlib.Path(base_path).absolute()

        if isinstance(base_path, (pathlib.Path, str)):
            self._search_paths = [pathlib.Path(base_path)]
        else:
            self._search_paths = [pathlib.Path(path) for path in base_path]

        self.reg_hash = {}
        self.reg_name = {}
        self.logger = logging.getLogger('aiida_vasp.utils.mock_code.MockRegistry')
        self._setup_logger()
        self.scan()

    def append_search_path(self, path: str | pathlib.Path) -> None:
        """Add a path to the list of search paths"""
        self.search_paths.append(pathlib.Path(path))

    @property
    def base_path(self) -> pathlib.Path:
        """Return the base repository path of the registry"""
        return self._search_paths[0]

    @property
    def search_paths(self) -> list[pathlib.Path]:
        """Return a list of all search paths"""
        return self._search_paths

    def _setup_logger(self, level: int = logging.INFO) -> None:
        """Setup the logger"""
        formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
        handler = logging.StreamHandler()
        handler.setLevel(level)
        handler.setFormatter(formatter)
        self.logger.setLevel(level)

    def scan(self) -> None:
        """
        Scan the base folder and locate input/output folders
        """
        for registry_path in self.search_paths:
            for output_folder in registry_path.glob('**/out'):
                calc_base_folder = output_folder.parent.absolute()
                self._register_folder(calc_base_folder)

    def get_path_by_hash(self, hash_val: str) -> pathlib.Path:
        """
        Return the output folder for a given hash
        """
        return pathlib.Path(self.reg_hash[hash_val])

    def get_path_by_name(self, name: str) -> pathlib.Path:
        """
        Return the output folder for a given hash
        """
        return pathlib.Path(self.reg_hash[self.reg_name[name]])

    def extract_calc_by_path(
        self, rel_path: pathlib.Path | str, dst_path: pathlib.Path | str, include_inputs: bool = True
    ) -> None:
        """
        Copy the content of a give hash to a destination.

        :param rel_path: The relative path of the calculation folder to be
          extracted.
        :param dst: The destination path to be extracted to - must already exists.
        """
        rel_path = pathlib.Path(rel_path)
        dst_path = pathlib.Path(dst_path)
        found = False
        for reg_path in self.search_paths:
            base_out = reg_path / rel_path / 'out'
            base_in = reg_path / rel_path / 'inp'

            # Not a valid folder - skip this
            if not (base_out.exists() and base_in.exists()):
                continue

            found = True
            # Copy the content of input and then the output folder
            paths = [base_in, base_out] if include_inputs else [base_out]
            for folder in paths:
                for fpath in folder.glob('*'):
                    if fpath.is_file():
                        shutil.copy2(fpath, dst_path)
                    # Directory - then copy the sub files - this only handles one level down
                    elif fpath.is_dir():
                        for subfile in fpath.glob('*'):
                            shutil.copy2(subfile, dst_path / fpath.name / subfile.name)
            break
        if not found:
            raise ValueError(f'The path give: {rel_path}, is not found in any search paths.')

    def extract_calc_by_hash(self, hash_val: str, dst: pathlib.Path | str, include_inputs: bool = False) -> None:
        """
        Extract an registerred calculation using hash.
        """
        self.extract_calc_by_path(self.get_path_by_hash(hash_val), dst, include_inputs)

    def upload_calc(
        self,
        folder: pathlib.Path,
        rel_path: pathlib.Path | str,
        excluded_object: list[str] | None = None,
        included_object: list[str] | None = None,
    ) -> None:
        """
        Register a calculation folder to primary search path of the registry
        """
        inp = list(INPUT_OBJECTS)
        excluded = list(DEFAULT_EXCLUDED)
        # Exclude certain objects
        if excluded_object:
            excluded.extend(excluded_object)
        # Include certain objects
        if included_object is not None:
            for name in included_object:
                if name in excluded:
                    del excluded[excluded.index(name)]

        # Check if the repository folder already exists
        repo_calc_base = self.base_path / pathlib.Path(self.get_upload_prefix() + str(rel_path))
        if repo_calc_base.exists():
            raise FileExistsError(f'There is already a directory at {repo_calc_base.resolve()}.')

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

    def _register_folder(self, calc_base: pathlib.Path) -> None:
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

    @classmethod
    def from_env(cls) -> MockRegistry:
        """Instantiate from environmental variable"""
        path = os.environ.get(f'{cls.CODE_NAME}_MOCK_CODE_BASE')
        if path is None:
            raise ValueError(f'The {cls.CODE_NAME}_MOCK_CODE_BASE environmental variable is not set!')

        paths = path.split(':')
        return cls(paths)

    @staticmethod
    def compute_hash(folder: pathlib.Path) -> str:
        """Compute the hash for a target folder"""
        raise NotImplementedError

    def upload_aiida_calc(
        self, calc_node: orm.CalcJobNode, rel_path: str | pathlib.Path, excluded_names: list[str] | None = None
    ) -> None:
        """Update a calculation into the registry"""
        raise NotImplementedError

    def upload_aiida_work(self, work_node: orm.WorkChainNode, rel_path: str | pathlib.Path) -> None:
        """Update all calculations run by an workflow into the registry"""
        raise NotImplementedError

    def get_upload_prefix(self) -> str:
        """Prefix of the name of the calculation folder"""
        prefix = os.environ.get(f'MOCK_{self.CODE_NAME}_UPLOAD_PREFIX')
        if prefix:
            prefix = prefix + '-'
        else:
            prefix = ''
        return prefix


class VaspMockRegistry(MockRegistry):
    """
    Registry of mock code for VASP

    The registry's base folder defaults to the `test_data` folder, but can be modified with the
    `MOCK_VASP_REG_BASE` environmental variable.
    """

    CODE_NAME = 'VASP'

    def upload_aiida_calc(
        self, calc_node: orm.CalcJobNode, rel_path: str | pathlib.Path, excluded_names: list[str] | None = None
    ) -> None:
        """
        Register an aiida calc_class
        """

        assert isinstance(calc_node, orm.CalcJobNode), f'{calc_node} is not an CalcJobNode!'

        # Check if the repository folder already exists
        repo_calc_base = self.base_path / pathlib.Path(self.get_upload_prefix() + str(rel_path))
        if repo_calc_base.exists():
            raise FileExistsError(f'There is already a directory at {repo_calc_base.resolve()}.')

        # Deposit the objects
        repo_calc_base.mkdir(parents=True)
        repo_in = repo_calc_base / 'inp'
        repo_out = repo_calc_base / 'out'
        repo_in.mkdir(parents=True)
        repo_out.mkdir(parents=True)

        exclude = list(DEFAULT_EXCLUDED)
        if excluded_names:
            exclude.extend(excluded_names)

        # Copy the input objects
        for obj in calc_node.base.repository.list_objects():
            if obj.name in exclude:
                continue
            copy_from_aiida(obj.name, calc_node, repo_in)

        # Copy the retrieved objects
        for obj in calc_node.outputs.retrieved.base.repository.list_objects():
            if obj.name in exclude:
                continue
            copy_from_aiida(obj.name, calc_node.outputs.retrieved, repo_out)

        self.logger.info('Calculation %s has been registered', calc_node)
        self._register_folder(repo_calc_base)

    def upload_aiida_work(self, work_node: orm.WorkChainNode, rel_path: str | pathlib.Path) -> None:
        """
        Upload all calculations in a workchain node
        """

        calc_class = CalculationFactory('vasp.vasp')
        neb_class = CalculationFactory('vasp.neb')
        to_upload = []
        for node in work_node.called_descendants:
            # Only upload VASP calculations
            if isinstance(node, orm.CalcJobNode) and (node.process_class in [calc_class, neb_class]):
                to_upload.append(node)
        to_upload.sort(key=lambda x: x.ctime)
        self.logger.info('Collected %s nodes to upload under name %s.', to_upload, rel_path)

        for idx, node in enumerate(to_upload):
            rel = pathlib.Path(rel_path) / f'calc-{idx:03d}'
            self.upload_aiida_calc(node, rel)
        self.logger.info('WorkChain %s has been uploaded.', work_node)

    @staticmethod
    def compute_hash(folder: pathlib.Path) -> str:
        """
        Compute the hash of a input folder
        """
        items = {}
        kpt_path = folder / 'KPOINTS'
        if kpt_path.is_file():
            kpoints = Kpoints(file_path=str(kpt_path))
            items['kpoints'] = kpoints.get_dict()
            items['kpoints'].pop('comment', None)

        incar_path = folder / 'INCAR'
        if incar_path.is_file():
            incar = Incar(file_path=str(incar_path), validate_tags=False)
            items['incar'] = incar.get_dict()

        poscar_path = folder / 'POSCAR'
        if not poscar_path.is_file():
            poscar_path = folder / '00/POSCAR'
        if poscar_path.is_file():
            poscar = Poscar(file_path=str(poscar_path))
            items['poscar'] = poscar.get_dict()
            items['poscar'].pop('comment', None)
        return get_hash(items)[0]


class MockVasp:
    """
    Mock VaspExecutable
    """

    def __init__(
        self,
        workdir: str | pathlib.Path,
        registry: VaspMockRegistry,
        vasp_cmd: str | list[str] | None = None,
        stdout_fname: str = 'vasp_output',
    ) -> None:
        """
        Mock VASP executable that copies over outputs from existing calculations.
        Inputs are hash and looked for.

        Notice that we do not set the hash value at init of workdir as we allow
        the unit of the MockVasp at any point, typically, you are preparing for
        a VASP calculation. Only when you execute VASP is the files checked, in this
        case when executing run. Thus, we calculate the hash of the workdir only then.

        If the `vasp_cmd` is provided the mock vasp will run the command if needed and
        upload the results to the registry. This can be useful for generating test/demo
        data.
        """
        self.workdir = workdir
        self.registry = registry
        if isinstance(vasp_cmd, str):
            vasp_cmd = [vasp_cmd]
        self.vasp_cmd = vasp_cmd
        self.stdout_fname = stdout_fname

    def run(self, debug: bool = True) -> None:
        """
        Run the mock vasp
        """

        if not os.listdir(self.workdir):
            # Directory is empty, no point of trying to find matching calcs
            raise ValueError('No input files given, so we can not find the associated test data.')

        hash_val = self.registry.compute_hash(self.workdir)
        if debug:
            print(f'Target hash value: {hash_val}')
        if hash_val in self.registry.reg_hash:
            self.registry.extract_calc_by_hash(hash_val, self.workdir)
        else:
            if debug:
                print(f'Registered hashes: {self.registry.reg_hash}')
            if self.vasp_cmd is not None:
                with open(pathlib.Path(self.workdir) / self.stdout_fname, 'w') as stdout_handle:
                    out = run(self.vasp_cmd, cwd=self.workdir, stdout=stdout_handle, check=False)
                    if out.returncode != 0:
                        raise ValueError(f'The command {self.vasp_cmd} failed with return code {out.returncode}')
                    self.registry.upload_calc(self.workdir, hash_val)
                if debug:
                    print(f'Uploaded current calculation with hash: {hash_val}')
            else:
                raise ValueError('The calculation is not registered.')

    @property
    def is_runnable(self) -> bool:
        """Check if the mock code can be executed."""
        hash_val = self.registry.compute_hash(self.workdir)
        if hash_val in self.registry.reg_hash:
            return True
        # Can we run vasp it self?
        if self.vasp_cmd is not None:
            out = run(['which', self.vasp_cmd[0]], check=False)
            if out.returncode == 0:
                return True
        return False


def copy_from_aiida(name: str, node: orm.Node, dst: pathlib.Path) -> None:
    """
    Copy objects from aiida repository.

    :param name: The full name (including the parent path) of the object.
    :param node: Node object for which the objects in the repo to be copied.
    :param dst: Path of the destination folder.

    This is a recursive function so directory copying also works.
    """
    obj = node.base.repository.get_object(name)

    # If it is a directory, copy the contents one by one
    if obj.file_type == FileType.DIRECTORY:
        for sub_obj in node.base.repository.list_objects(name):
            copy_from_aiida(os.path.join(name, sub_obj.name), node, dst)
    else:
        # Anything else
        with node.base.repository.open(name) as fsource:
            # Make parent directory if needed
            frepo_path = dst / name
            pathlib.Path(frepo_path.parent).mkdir(exist_ok=True, parents=True)
            # Write the object
            with open(frepo_path, 'w', encoding='utf8') as fdst:
                shutil.copyfileobj(fsource, fdst)
