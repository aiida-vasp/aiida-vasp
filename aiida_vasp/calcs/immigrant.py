# pylint: disable=abstract-method
# explanation: pylint wrongly complains about (aiida) Node not implementing query
"""VaspImmigrant calculation: Immigrate externally run VASP calculations into AiiDA."""
import os

from py import path as py_path  # pylint: disable=no-name-in-module,no-member
from aiida.common.folders import SandboxFolder

from aiida_vasp.calcs.vasp import VaspCalculation
from aiida_vasp.io.incar import IncarIo
from aiida_vasp.io.kpoints import KpParser
from aiida_vasp.io.poscar import PoscarParser


## TODO: add a test
class VaspImmigrant(VaspCalculation):
    """
    Takes a VASP run directory as input and creates inputs and outputs like for a VaspCalculation.

    Current limitations:
        * KPOINTS must be given as a separate file and not in INCAR only
        * settings has to be given manually
    """

    def _process_remote_workdir(remote_workdir):
        if remote_workdir:
            self.set_remote_workdir(remote_workdir)
        elif not self._get_attr('remote_workdir', None):
            raise InputValidationEror(
                'expeced keyword parameter ``remote_workdir`` as it has not been set previously with ``.set_remote_workdir()``')

    def create_input_nodes(self, open_transport, remote_workdir=None):
        """
        Create a Calculation based on a folder in which VASP was run.

        :param open_transport: An open instance of the transport used by the calculations's computer.
        :param remote_workdir: optional only if ``remote_workdir`` has already been set previously using ``set_remote_workdir()``
                Absolute path to the directory where the calculation was run. All input and output files are expected to be found
                at the given path by the given transport instance.
        """
        self._process_remote_workdir(remote_workdir)
        # Copy the input file and psuedo files to a temp folder for parsing.
        with SandboxFolder() as sandbox:

            sandbox_path = py_path.local(sandbox.abspath)
            ## Copy the input file to the temp folder.
            remote_path = os.path.join(self._get_remote_workdir(), self._INPUT_FILE_NAME)
            open_transport.get(remote_path, folder.abspath)

            local_incar = sandbox_path.join(self._INPUT_FILE_NAME)
            incar_io = IncarIo(file_path=local_incar.strpath)
            self.use_parameters(incar_io.get_param_node())

            ## Create the input structure from POSCAR
            ## TODO: allow user input for symbols / kind names
            remote_path = os.path.join(self._get_remote_workdir(), 'POSCAR')
            local_poscar = sandbox_path.join('POSCAR')
            poscar_in = PoscarParser(file_path=local_poscar.strpath)
            structure = poscar_in.get_quantity('poscar', None)['structure']
            self.use_structure(structure)

            ## Do the same with KPOINTS
            remote_path = os.path.join(self._get_remote_workdir(), 'KPOINTS')
            open_transport.get(remote_path, folder.abspath)
            local_kpoints = sandbox_path.join('KPOINTS')
            kpoints_in = KpParser(file_path=local_kpoints.strpath)
            kpoints = kpoints_in.get_quantity('kpoints', None)['kpoints']
            kpoints.set_cell_from_structure(structure)
            self.use_kpoints(kpoints)

            ## Now get the PotcarData associated with the used POTCARs
            remote_path = os.path.join(self._get_remote_workdir(), 'POTCAR')
            open_transport.get(remote_path, folder.abspath)
            local_potcar = sandbox_path.join('POTCAR')
            multi_potcar_io = MultiPotcarIo.read(local_potcar.strpath)
            for kind_name, potcar in multi_potcar_io.get_potentials_dict().items():
                self.use_potential(potcar, kind=kind_name)
