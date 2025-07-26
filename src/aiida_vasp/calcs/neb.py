"""
Module for settings up NEB calculations
"""

import os
from pathlib import Path
from typing import Any

from aiida import orm
from aiida.common import CalcInfo
from aiida.common.exceptions import InputValidationError
from aiida.common.folders import Folder

from aiida_vasp.calcs.vasp import VaspCalculation, ordered_unique_symbols
from aiida_vasp.data.chargedensity import ChargedensityData
from aiida_vasp.data.potcar import PotcarData
from aiida_vasp.data.wavefun import WavefunData
from aiida_vasp.parsers.content_parsers.poscar import PoscarParser


class VaspNEBCalculation(VaspCalculation):
    """
    NEB calculations using VASP

    Calculations for performing NEB calculations.
    NEB calculations requires standard VASP inputs, but POSCAR are placed in
    folder names 00, 01, 02... N for N-1 number of images.

    Input frames should be placed under the ``neb_images`` input namespace as a dictionary like::

      {
          'image_00': structure_1,
          'image_01': structure_2
          ....
      }

    Output of individual frames are placed in the corresponding namespace under the same convention.
    """

    # Use stdout as the name for the diverted stdout to be consistent with those in the image folders
    _VASP_OUTPUT = 'stdout'
    _ALWAYS_RETRIEVE_LIST = ['OUTCAR', 'vasprun.xml', _VASP_OUTPUT]
    _PER_IMAGE_ALWAYS_RETRIEVE_LIST = ['OUTCAR', 'CONTCAR', _VASP_OUTPUT]
    _query_type_string = 'vasp.neb'
    _plugin_type_string = 'vasp.neb'
    _default_parser = 'vasp.neb'

    @classmethod
    def define(cls, spec: Any) -> None:
        super(VaspNEBCalculation, cls).define(spec)
        # NEB calculation does not have the structure input port
        spec.inputs.pop('structure')

        # Define the inputs.
        # options is passed automatically.
        spec.input(
            'parameters',
            valid_type=orm.Dict,
            help='The VASP input parameters (INCAR).',
        )
        spec.input(
            'dynamics',
            valid_type=orm.Dict,
            help='The VASP parameters related to ionic dynamics, e.g. flags to set the selective dynamics',
            required=False,
        )
        spec.input(
            'initial_structure',
            valid_type=(orm.StructureData, orm.CifData),
            help='The input structure (POSCAR) for initial image.',
        )
        spec.input(
            'final_structure',
            valid_type=(orm.StructureData, orm.CifData),
            help='The input structure (POSCAR) for the final image.',
        )
        spec.input(
            'metadata.options.withmpi',
            valid_type=bool,
            default=True,
            help='Set the calculation to use mpi',
        )
        spec.input_namespace(
            'neb_images',
            valid_type=(orm.StructureData, orm.CifData),
            help='Starting structure for the NEB images',
            dynamic=True,
        )
        # Need namespace on this as it should also accept keys that are of `kind`. These are unknown
        # until execution.
        spec.input_namespace(
            'potential',
            valid_type=PotcarData,
            help='The potentials (POTCAR).',
            dynamic=True,
        )
        spec.input(
            'kpoints',
            valid_type=orm.KpointsData,
            help='The kpoints to use (KPOINTS).',
        )
        spec.input_namespace(
            'charge_density',
            dynamic=True,
            valid_type=ChargedensityData,
            required=False,
            help='The charge density. (CHGCAR)',
        )
        spec.input_namespace(
            'wavefunctions',
            valid_type=WavefunData,
            dynamic=True,
            required=False,
            help='The wave function coefficients. (WAVECAR)',
        )
        spec.input(
            'settings',
            valid_type=orm.Dict,
            required=False,
            help='Additional parameters not related to VASP itself.',
        )
        spec.input(
            'metadata.options.parser_name',
            default=cls._default_parser,
        )

        # Define outputs.
        # remote_folder and retrieved are passed automatically

        # Turn on dynamic outputs to avoid problem when catching is used
        # it appears the catching implementation does not work with calcjobs with namespace outputs
        spec.outputs.dynamic = True
        spec.output_namespace(
            'structure',
            required=True,
            valid_type=orm.StructureData,
            help='NEB images',
            dynamic=True,
        )
        spec.output(
            'misc',
            valid_type=orm.Dict,
            required=True,
            help='Per-image misc output.',
        )
        spec.output(
            'arrays',
            valid_type=orm.ArrayData,
            required=False,
            help='Array data from each image.',
        )
        spec.output_namespace(
            'trajectory',
            valid_type=orm.TrajectoryData,
            dynamic=True,
            required=False,
            help='The NEB trajectory.',
        )

    def prepare_for_submission(self, folder: Folder) -> CalcInfo:
        """
        Add all files to the list of files to be retrieved.

        Notice that we here utilize both the retrieve batch of files, which are always stored after retrieval and
        the temporary retrieve list which is automatically cleared after parsing.
        """
        calcinfo = super().prepare_for_submission(folder)

        nimages = len(self.inputs.neb_images)
        nimage_keys = sorted(list(self.inputs.neb_images.keys()))

        image_folders = []

        # Iterate though each image_folder that needs to be setup
        for i in range(nimages + 2):
            folder_id = f'{i:02d}'
            image_folder = Path(folder.get_abs_path(folder_id))
            image_folder.mkdir()
            poscar = str(image_folder / 'POSCAR')
            if i == 0:
                # Write the initial image
                self.write_neb_poscar(self.inputs.initial_structure, poscar)
            elif i == nimages + 1:
                # Write the final image
                self.write_neb_poscar(self.inputs.final_structure, poscar)
            else:
                # Write NEB images
                img_key = nimage_keys[i - 1]
                self.write_neb_poscar(self.inputs.neb_images[img_key], poscar)
                image_folders.append(folder_id)

                # Link with singlefile WAVECAR/CHGCAR is needed
                if self._need_wavecar():
                    wavecar = self.inputs.wave_functions[img_key]
                    dst = folder_id + '/' + 'WAVECAR'
                    calcinfo.local_copy_list.append((wavecar.uuid, wavecar.filename, dst))

                if self._need_chgcar():
                    chgcar = self.inputs.charge_densities[img_key]
                    dst = folder_id + '/' + 'CHGCAR'
                    calcinfo.local_copy_list.append((chgcar.uuid, chgcar.filename, dst))
        try:
            store = self.inputs.settings.base.attributes.get('ALWAYS_STORE', default=True)
        except AttributeError:
            store = True

        try:
            additional_retrieve_list = self.inputs.settings.base.attributes.get(
                'PER_IMAGE_ADDITIONAL_RETRIEVE_LIST', default=[]
            )
        except AttributeError:
            additional_retrieve_list = []
        try:
            additional_retrieve_temp_list = self.inputs.settings.base.attributes.get(
                'PER_IMAGE_ADDITIONAL_RETRIEVE_TEMPORARY_LIST', default=[]
            )
        except AttributeError:
            additional_retrieve_temp_list = []

        if store:
            calcinfo.retrieve_list.extend(
                image_folder_paths(image_folders, set(self._PER_IMAGE_ALWAYS_RETRIEVE_LIST + additional_retrieve_list))
            )
            calcinfo.retrieve_temporary_list.extend(image_folder_paths(image_folders, additional_retrieve_temp_list))
        else:
            calcinfo.retrieve_temporary_list.extend(
                image_folder_paths(
                    image_folders, set(self._PER_IMAGE_ALWAYS_RETRIEVE_LIST + additional_retrieve_temp_list)
                )
            )
            calcinfo.retrieve_list.extend(image_folder_paths(image_folders, additional_retrieve_list))

        # self.logger.warning('Calcinfo: {}'.format(calcinfo))

        return calcinfo

    def _structure(self) -> orm.StructureData:
        """
        Get the input structure as AiiDa StructureData.

        This is required in order to support CifData as input as well.
        """
        structure = self.inputs.initial_structure
        if not hasattr(structure, 'get_pymatgen'):
            structure = orm.StructureData(ase=structure.get_ase())
        return structure

    def remote_copy_restart_folder(self) -> list[tuple[str, str, str]]:
        """
        Add all files required for restart to the list of files to be copied from the previous calculation.

        For NEB calculations, the CHGCAR and WAVECAR needs to be copied for each of the actual images.
        """
        restart_folder = self.inputs.restart_folder
        computer = self.node.computer
        included = ['CHGCAR', 'WAVECAR']
        nimages = len(self.inputs.neb_images)
        copy_list = []
        for image_id in range(1, nimages + 1):
            fdname = f'{image_id:02d}'

            existing_files = restart_folder.listdir(fdname)
            for name in included:
                if name not in existing_files:
                    # Here we simple issue an warning as the requirement of files will be explicitly checked by
                    # `write_additional` method
                    self.report(f'WARNING: File {name} does not exist in the restart folder.')
                else:
                    copy_list.append(
                        (
                            computer.uuid,
                            os.path.join(restart_folder.get_remote_path(), fdname, name),
                            os.path.join(fdname, name),
                        )
                    )

        return copy_list

    def write_neb_poscar(
        self, structure: orm.StructureData | orm.CifData, dst: str, positions_dof: list | None = None
    ) -> None:
        """
        Write the POSCAR.

        Passes the structures node (StructureData) to the POSCAR parser for
        preparation and writes to dst.

        :param dst: absolute path of the file to write to

        """
        settings = self.inputs.get('settings')
        settings = settings.get_dict() if settings else {}
        poscar_precision = settings.get('poscar_precision', 10)
        if positions_dof is not None:
            options = {'positions_dof': positions_dof}
        else:
            options = None
        poscar_parser = PoscarParser(data=ensure_structure_data(structure), precision=poscar_precision, options=options)
        poscar_parser.write(dst)

    def verify_inputs(self) -> None:
        """
        Verify the order of elements
        """
        super().verify_inputs()
        last_order = None
        last_num_sites = None
        for structure in list(self.inputs.neb_images.values()) + [
            self.inputs.initial_structure,
            self.inputs.final_structure,
        ]:
            # Convert to StructureData from CifData on demand....
            if not hasattr(structure, 'get_pymatgen'):
                structure_data = orm.StructureData(ase=structure.get_ase())
            else:
                structure_data = structure

            num_sites = len(structure_data.sites)

            order = ordered_unique_symbols(structure_data)
            if last_order is not None and order != last_order:
                raise InputValidationError('The input structures have non-euqal element orders.')
            if last_order is not None and num_sites != last_num_sites:
                raise InputValidationError('The input structures have non-equal number of atoms.')
            last_order = order
            last_num_sites = num_sites

    def write_incar(self, dst: str, validate_tags: bool = False) -> None:
        """
        Write the INCAR without tag validation.
        Validation is performed at `parsevasp` level and VTST tags are not included.

        :param dst: absolute path of the file to write to
        """
        super().write_incar(dst, validate_tags=validate_tags)


def image_folder_paths(image_folders: list[str], retrieve_names: list[str]) -> list[list[str | int]]:
    """
    Return a list of folders paths to be retrieved
    """
    retrieve_list = []
    for key in retrieve_names:
        for fdname in image_folders:
            # Skip the stdout for the first image - it is diverted to the root folder
            if int(fdname) == 1 and key == VaspNEBCalculation._VASP_OUTPUT:  # pylint: disable=protected-access
                continue
            # Need to use the tuple format to keep the sub directory structure
            retrieve_list.append([fdname + '/' + key, '.', 2])
    return retrieve_list


def ensure_structure_data(structure: orm.StructureData | orm.CifData) -> orm.StructureData:
    """
        Get the input structure as AiiDA StructureData.

    This is required in order to support CifData as input as well.
    """
    if not hasattr(structure, 'get_pymatgen'):
        structure = orm.StructureData(ase=structure.get_ase())
    return structure
