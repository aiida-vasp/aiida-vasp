"""
Module for importing existing VASP calculations by reading input files and running
dummy calculations through aiida-core.
"""

from __future__ import annotations

from logging import getLogger
from pathlib import Path
from typing import Any

import plumpy
from aiida import orm
from aiida.common import InputValidationError
from aiida.common.extendeddicts import AttributeDict
from aiida.common.folders import SandboxFolder
from aiida.engine import run_get_node, submit

from aiida_vasp.calcs.vasp import VaspCalculation
from aiida_vasp.data.chargedensity import ChargedensityData
from aiida_vasp.data.potcar import PotcarData
from aiida_vasp.data.wavefun import WavefunData
from aiida_vasp.parsers.content_parsers.incar import IncarParser
from aiida_vasp.parsers.content_parsers.kpoints import KpointsParser
from aiida_vasp.parsers.content_parsers.poscar import PoscarParser
from aiida_vasp.parsers.content_parsers.potcar import MultiPotcarIo
from aiida_vasp.parsers.vasp import get_kpoints_node, get_structure_node

logger = getLogger(__name__)


class VaspCalcImporter:
    """
    Importer for VASP calculations.

    The importer is responsible for parsing an existing folder to recreate the inputs nodes,
    then the calculation can be submitted to the AiiDA engine with a `remote_folder` input node
    so parsing is triggered immediately without actually running the calculation.
    """

    @classmethod
    def get_builder_from_folder(
        cls,
        code: orm.AbstractCode,
        remote_path: str | None = None,
        remote_folder: orm.RemoteData | None = None,
        options: None | dict[str, Any] = None,
        settings: None | dict[str, Any] = None,
        potential_family: None | str = None,
        potential_mapping: None | dict[str, str] = None,
        include_wavecar: bool = False,
        include_chgcar: bool = False,
        stdout_file_name: str = 'vasp_output',
        dummy_resources: dict[str, str] | None = None,
        **kwargs,
    ) -> AttributeDict:
        """
        Create inputs to launch a dummy calculation from a code and a remote path on the associated computer.

        If POTCAR does not exist, the provided ``potential_family`` and
        ``potential_mapping`` are used to link potential to inputs. In this
        case, at least ``potential_family`` has to be provided. Unless
        ``potential_mapping``, this mapping is generated from structure, i.e.,

        ::

            potential_mapping = {element: element for element in structure.get_kind_names()}

        :param code: a Code instance for the code originally used.
        :param remote_path: Directory where VASP inputs and outputs are stored on the remote machine.
        :param settings: dict. This is used as the input port of VaspCalculation.
        :param potential_family: str potential family to load the POTCAR from, if the POTCAR is missing.
        :param potential_mapping: dict mapping of elements to pseudopotentials, if the POTCAR is missing.
        :param include_wavecar: bool. Try to read WAVECAR.
        :param include_chgcar: bool. Try to read CHGCAR.
        :param stdout_file_name: str. Name of the stdout file to look for.
        :param dummy_resources: dict[str, str]. Dummy resources to use for the calculation.

        """

        computer = code.computer
        # Create the remote folder
        if remote_folder is None:
            if computer is None or remote_path is None:
                raise ValueError('If remote_folder is not provided, both computer and remote_path have to be set.')
            remote_folder = orm.RemoteData(remote_path, computer=computer)

        builder = VaspCalculation.get_builder()
        if options:
            builder.metadata.options = options
        builder.remote_folder = remote_folder
        if settings:
            builder.settings = settings
        _remote_workdir = Path(remote_folder.get_remote_path())
        with remote_folder.computer.get_transport() as transport:
            with SandboxFolder() as sandbox:
                sandbox_path = Path(sandbox.abspath)
                # TODO: use get_async instead?
                transport.get(str(_remote_workdir / 'INCAR'), str(sandbox_path))
                transport.get(str(_remote_workdir / 'POSCAR'), str(sandbox_path))
                transport.get(str(_remote_workdir / 'POTCAR'), str(sandbox_path), ignore_nonexisting=True)
                transport.get(str(_remote_workdir / 'KPOINTS'), str(sandbox_path))
                if include_wavecar:
                    transport.get(str(_remote_workdir / 'WAVECAR'), str(sandbox_path), ignore_nonexisting=True)
                    builder.wavefunctions = get_wavecar_input(sandbox_path)
                if include_chgcar:
                    transport.get(str(_remote_workdir / 'CHGCAR'), str(sandbox_path), ignore_nonexisting=True)
                    builder.charge_density = get_chgcar_input(sandbox_path)
                # Check if there is a stdout file
                if not transport.isfile(str(_remote_workdir / stdout_file_name)):
                    logger.warning(
                        f'No stdout file `{stdout_file_name}` found in remote folder.'
                        f' Creating a empty stdout file {VaspCalculation._VASP_OUTPUT} '
                        'to be parsed by aiida-vasp.'
                    )
                    transport.exec_command_wait(f'touch {_remote_workdir / VaspCalculation._VASP_OUTPUT}')
                # The stdout file has a different name - we copy it over
                elif stdout_file_name != VaspCalculation._VASP_OUTPUT:
                    transport.exec_command_wait(
                        f'cp {_remote_workdir / stdout_file_name} {_remote_workdir / VaspCalculation._VASP_OUTPUT}'
                    )
                    logger.info(
                        f'Copied `{stdout_file_name}` found in remote folder to {VaspCalculation._VASP_OUTPUT}.'
                    )

                builder.parameters = get_incar_input(sandbox_path)
                builder.structure = get_poscar_input(sandbox_path)
                builder.kpoints = get_kpoints_input(sandbox_path, structure=builder.structure)
                builder.potential = get_potcar_input(
                    sandbox_path,
                    structure=builder.structure,
                    potential_family=potential_family,
                    potential_mapping=potential_mapping,
                )
        dummy_resources = dummy_resources or {'num_machine': 1}
        builder.metadata.options.resources = dummy_resources
        builder.code = code
        return builder

    @classmethod
    def run_import(
        cls,
        code: orm.AbstractCode,
        remote_path: str | None = None,
        remote_folder: orm.RemoteData | None = None,
        options: None | dict[str, Any] = None,
        settings: None | dict[str, Any] = None,
        potential_family: None | str = None,
        potential_mapping: None | dict[str, str] = None,
        include_wavecar: bool = False,
        include_chgcar: bool = False,
        stdout_file_name: str = 'vasp_output',
        dummy_resources: dict[str, str] | None = None,
        **kwargs,
    ) -> VaspCalculation:
        """
        Run the import process for a VASP calculation.

        This is blocking function - the import process is done within the current python session.

        :param code: a Code instance for the code originally used.
        :param remote_path: Directory where VASP inputs and outputs are stored on the remote machine.
        :param remote_folder: The remote folder containing the VASP calculation files.
        :param options: Additional options for the calculation.
        :param settings: dict. This is used as the input port of VaspCalculation.
        :param potential_family: str potential family to load the POTCAR from, if the POTCAR is missing.
        :param potential_mapping: dict mapping of elements to pseudopotentials, if the POTCAR is missing.
        :param include_wavecar: bool. Try to read WAVECAR.
        :param include_chgcar: bool. Try to read CHGCAR.
        :param stdout_file_name: str. Name of the stdout file to look for.
        :param dummy_resources: dict[str, str]. Dummy resources to use for the calculation.
        :return: The VASP calculation instance.
        """
        builder = cls.get_builder_from_folder(
            code=code,
            remote_path=remote_path,
            remote_folder=remote_folder,
            options=options,
            settings=settings,
            potential_family=potential_family,
            potential_mapping=potential_mapping,
            include_wavecar=include_wavecar,
            include_chgcar=include_chgcar,
            stdout_file_name=stdout_file_name,
            dummy_resources=dummy_resources,
            **kwargs,
        )
        return run_get_node(builder)[1]

    @classmethod
    def run_import_daemon(
        cls,
        code: orm.AbstractCode,
        remote_path: str | None = None,
        remote_folder: orm.RemoteData | None = None,
        options: None | dict[str, Any] = None,
        settings: None | dict[str, Any] = None,
        potential_family: None | str = None,
        potential_mapping: None | dict[str, str] = None,
        include_wavecar: bool = False,
        include_chgcar: bool = False,
        stdout_file_name: str = 'vasp_output',
        dummy_resources: dict[str, str] | None = None,
        **kwargs,
    ) -> plumpy.Process:
        """
        Submit the import process for a VASP calculation to the daemon.

        This is a none-blocking function - the import process is carried out by the daemon.

        :param code: a Code instance for the code originally used.
        :param remote_path: Directory where VASP inputs and outputs are stored on the remote machine.
        :param remote_folder: The remote folder containing the VASP calculation files.
        :param options: Additional options for the calculation.
        :param settings: dict. This is used as the input port of VaspCalculation.
        :param potential_family: str potential family to load the POTCAR from, if the POTCAR is missing.
        :param potential_mapping: dict mapping of elements to pseudopotentials, if the POTCAR is missing.
        :param include_wavecar: bool. Try to read WAVECAR.
        :param include_chgcar: bool. Try to read CHGCAR.
        :param stdout_file_name: str. Name of the stdout file to look for.
        :param dummy_resources: dict[str, str]. Dummy resources to use for the calculation.
        :return: The process instance.
        """
        builder = cls.get_builder_from_folder(
            code=code,
            remote_path=remote_path,
            remote_folder=remote_folder,
            options=options,
            settings=settings,
            potential_family=potential_family,
            potential_mapping=potential_mapping,
            include_wavecar=include_wavecar,
            include_chgcar=include_chgcar,
            stdout_file_name=stdout_file_name,
            dummy_resources=dummy_resources,
            **kwargs,
        )
        return submit(builder)


def get_incar_input(dir_path: Path) -> orm.Dict:
    """Create a node that contains the INCAR content."""

    with open(str(dir_path / 'INCAR'), 'r', encoding='utf8') as handler:
        incar_parser = IncarParser(handler=handler, raise_errors=True)
    node = orm.Dict(incar_parser.incar)

    return node


def get_poscar_input(dir_path: Path) -> orm.StructureData:
    """Create a node that contains the POSCAR content."""

    with open(str(dir_path / 'POSCAR'), 'r', encoding='utf8') as handler:
        poscar_parser = PoscarParser(handler=handler, raise_errors=True)
    return get_structure_node(poscar_parser.structure)


def get_potcar_input(
    dir_path: Path,
    structure: orm.StructureData | None = None,
    potential_family: str | None = None,
    potential_mapping: dict | None = None,
) -> dict[str, PotcarData]:
    """Read potentials from POTCAR or set it up from a structure."""
    local_potcar = dir_path / 'POTCAR'
    structure = structure or get_poscar_input(dir_path)
    potentials = {}
    if local_potcar.exists():
        potentials = MultiPotcarIo.read(str(local_potcar)).get_potentials_dict(structure)
        potentials = {kind: potentials[kind] for kind in potentials}
    elif potential_family:
        potentials = PotcarData.get_potcars_from_structure(structure, potential_family, mapping=potential_mapping)
    else:
        raise InputValidationError('no POTCAR found in remote folder and potential_family was not passed')

    return potentials


def get_kpoints_input(dir_path: Path, structure: orm.StructureData | None = None) -> orm.KpointsData:
    """Create a node that contains the KPOINTS content."""
    structure = structure or get_poscar_input(dir_path)
    with open(str(dir_path / 'KPOINTS'), 'r', encoding='utf8') as handler:
        kpoints_parser = KpointsParser(handler=handler, raise_errors=True)
    return get_kpoints_node(kpoints_parser.kpoints, structure.cell)


def get_chgcar_input(dir_path: Path) -> ChargedensityData:
    """Include CHGCAR as input"""
    node = ChargedensityData(str(dir_path / 'CHGCAR'))
    return node


def get_wavecar_input(dir_path: Path) -> WavefunData:
    """Include WAVECAR as input"""
    node = WavefunData(str(dir_path / 'WAVECAR'))

    return node
