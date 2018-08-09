# pylint: disable=abstract-method
# explanation: pylint wrongly complains about (aiida) Node not implementing query
"""VaspImmigrant calculation: Immigrate externally run VASP calculations into AiiDA."""
import os

from py import path as py_path  # pylint: disable=no-name-in-module,no-member
from aiida.common.folders import SandboxFolder
from aiida.common.exceptions import InputValidationError

from aiida_vasp.calcs.vasp import VaspCalculation
from aiida_vasp.data.potcar import PotcarData
from aiida_vasp.io.incar import IncarParser
from aiida_vasp.io.kpoints import KpParser
from aiida_vasp.io.poscar import PoscarParser
from aiida_vasp.io.potcar import MultiPotcarIo
from aiida_vasp.utils.aiida_utils import get_data_node


def get_quantity_node(parser, quantity):
    return parser.get_quantity(quantity, None)[quantity]


class VaspImmigrant(VaspCalculation):
    """
    Takes a VASP run directory as input and creates inputs and outputs like for a VaspCalculation.

    Current limitations:
        * KPOINTS must be given as a separate file and not in INCAR only
        * settings has to be given manually
    """

    def _process_remote_workdir(self, remote_workdir):
        if remote_workdir:
            self.set_remote_workdir(remote_workdir)
        elif not self.get_attr('remote_workdir', None):
            raise InputValidationError(
                'expected keyword parameter ``remote_workdir`` as it has not been set previously with ``.set_remote_workdir()``')

    def set_remote_workdir(self, remote_workdir):
        self._set_attr('remote_workdir', remote_workdir)

    def _create_incar_input(self, open_transport, sandbox_path):
        """Copy the INCAR file, create a parameter node from it and use that."""
        remote_path = os.path.join(self._get_remote_workdir(), 'INCAR')
        open_transport.get(remote_path, sandbox_path.strpath)
        local_incar = sandbox_path.join('INCAR')
        incar_in = IncarParser(file_path=local_incar.strpath)
        incar = get_quantity_node(incar_in, 'incar')
        self.use_parameters(incar)
        return incar

    def _create_poscar_input(self, open_transport, sandbox_path):
        """Copy the POSCAR file, create a structure node from it and use that."""
        remote_path = os.path.join(self._get_remote_workdir(), 'POSCAR')
        open_transport.get(remote_path, sandbox_path.strpath)
        local_poscar = sandbox_path.join('POSCAR')
        poscar_in = PoscarParser(file_path=local_poscar.strpath)
        structure = get_quantity_node(poscar_in, 'poscar-structure')
        self.use_structure(structure)
        return structure

    def _create_potcar_input(self, open_transport, sandbox_path, structure, potcar_spec=None):
        """Copy the POTCAR files, retrieve the potcar nodes for them and use those."""
        remote_path = os.path.join(self._get_remote_workdir(), 'POTCAR')
        open_transport.get(remote_path, sandbox_path.strpath, ignore_nonexisting=True)
        local_potcar = sandbox_path.join('POTCAR')
        multi_potcar_io = None
        if local_potcar.exists():
            multi_potcar_io = MultiPotcarIo.read(local_potcar.strpath)
        elif potcar_spec:
            potentials_dict = PotcarData.get_potcars_dict(structure.get_kind_names(), potcar_spec['family'], potcar_spec['map'])
            multi_potcar_io = MultiPotcarIo.from_structure(structure, potentials_dict)
        else:
            raise InputValidationError('no POTCAR found in remote folder and potcar_spec was not passed')
        for kind_name, potcar in multi_potcar_io.get_potentials_dict(structure).items():
            self.use_potential(potcar, kind=kind_name)

        return multi_potcar_io

    def _create_kpoints_input(self, open_transport, sandbox_path, structure):
        """Copy the KPOINTS file, make it into a kpoints node and use that."""
        remote_path = os.path.join(self._get_remote_workdir(), 'KPOINTS')
        open_transport.get(remote_path, sandbox_path.strpath)
        local_kpoints = sandbox_path.join('KPOINTS')
        kpoints_in = KpParser(file_path=local_kpoints.strpath)
        kpoints = get_quantity_node(kpoints_in, 'kpoints-kpoints')
        kpoints.set_cell_from_structure(structure)
        self.use_kpoints(kpoints)
        return kpoints

    def create_input_nodes(self, open_transport, remote_workdir=None, settings_dict=None, potcar_spec=None):
        """
        Create a Calculation based on a folder in which VASP was run.

        :param open_transport: An open instance of the transport used by the calculations's computer.
        :param remote_workdir: optional only if ``remote_workdir`` has already been set previously using ``set_remote_workdir()``
                Absolute path to the directory where the calculation was run. All input and output files are expected to be found
                at the given path by the given transport instance.
        :param settings_dict: An optional dictionary with a dictionary for the ``settings`` input.
        """
        self._process_remote_workdir(remote_workdir)
        # Copy the input file and psuedo files to a temp folder for parsing.
        with SandboxFolder() as sandbox:

            sandbox_path = py_path.local(sandbox.abspath)

            _ = self._create_incar_input(open_transport, sandbox_path)
            structure = self._create_poscar_input(open_transport, sandbox_path)
            _ = self._create_potcar_input(open_transport, sandbox_path, structure, potcar_spec=potcar_spec)
            _ = self._create_kpoints_input(open_transport, sandbox_path, structure)

        if settings_dict:
            self.use_settings(get_data_node('parameter', dict=settings_dict))

        self._set_attr('input_nodes_created', True)
