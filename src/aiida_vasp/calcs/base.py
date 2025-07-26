"""
Base calculation class.

Base and meta classes for VASP calculation classes.
"""

from __future__ import annotations

# pylint: disable=abstract-method,invalid-metaclass,ungrouped-imports
# explanation: pylint wrongly complains about Node not implementing query
import os
from typing import TYPE_CHECKING, Any

from aiida import orm
from aiida.common import CalcInfo, CodeInfo, ValidationError
from aiida.engine import CalcJob

if TYPE_CHECKING:
    from aiida.common.folders import Folder
    from aiida.engine.processes.calcjobs.calcjob import CalcJobProcessSpec


class VaspCalcBase(CalcJob):
    """
    Base class of all calculations utilizing VASP.

    * Defines internal parameters common to all vasp calculations.
    * provides a basic, extendable implementation of _prepare_for_submission
    * provides hooks, so subclasses can extend the behavior without having to reimplement common functionality
    """

    _default_parser = 'vasp.vasp'

    @classmethod
    def define(cls, spec: CalcJobProcessSpec) -> None:
        super(VaspCalcBase, cls).define(spec)
        spec.input(
            'restart_folder',
            valid_type=orm.RemoteData,
            help='A remote folder to restart from if need be',
            required=False,
        )

    @classmethod
    def max_retrieve_list(cls) -> list[str]:
        """Return a list of all possible output objects from a VASP run."""
        retrieve_list = [
            'CHG',
            'CHGCAR',
            'AECCAR0',
            'AECCAR1',
            'AECCAR2',
            'ELFCAR',
            'PARCHG',  # Density related
            'CONTCAR',
            'XDATCAR',
            'PCDAT',  # Structure related
            'DOSCAR',
            'EIGENVAL',
            'PROCAR',  # Electronic structure related
            'IBZKPT',  # Irreducible k-points
            'LOCPOT',  # Potential related
            'BSEFATBAND',  # Eigenvectors of the BSE matrix
            'WAVECAR',
            'WAVEDER',
            'PROOUT',
            'TMPCAR',
            'W*.tmp',
            'WFULL*.tmp',  # Wave function related properties
            'wannier90*',  # Wannier90 related
            'OSZICAR',  # Convergence related
            'REPORT',  # Output of molecular dynamics runs
            'STOPCAR',  # Controlled stopping
            'vasprun.xml',
            'OUTCAR',
        ]
        return retrieve_list

    def prepare_for_submission(self, folder: Folder) -> CalcInfo:  # pylint: disable=arguments-differ
        """
        Writes the four minimum outputs: INCAR, POSCAR, POTCAR, KPOINTS.

        Delegates the construction and writing / copying to write_<object> methods.
        That way, subclasses can use any form of input nodes and just
        have to implement the write_xxx method accordingly.
        Subclasses can extend by calling the super method and if necessary
        modifying it's output CalcInfo before returning it.
        """
        # write input objects
        incar = folder.get_abs_path('INCAR')
        structure = folder.get_abs_path('POSCAR')
        potentials = folder.get_abs_path('POTCAR')
        kpoints = folder.get_abs_path('KPOINTS')

        remote_copy_list = []

        self.verify_inputs()
        if self._is_restart():
            remote_copy_list.extend(self.remote_copy_restart_folder())
        self.write_incar(incar)
        self.write_poscar(structure)
        self.write_potcar(potentials)
        self.write_kpoints(kpoints)

        # calcinfo
        calcinfo = CalcInfo()
        calcinfo.uuid = self.uuid
        calcinfo.retrieve_list = self.max_retrieve_list()
        codeinfo = CodeInfo()
        codeinfo.code_uuid = self.inputs.code.uuid
        codeinfo.code_pk = self.inputs.code.pk
        calcinfo.codes_info = [codeinfo]
        calcinfo.remote_copy_list = remote_copy_list
        # here we need to do the charge density and wave function copy
        # as we need access to the calcinfo
        calcinfo.local_copy_list = []
        self.write_additional(folder, calcinfo)

        return calcinfo

    def remote_copy_restart_folder(self) -> list[tuple[str, str, str]]:
        """Add all objects required for restart to the list of objects to be copied from the previous calculation."""
        restart_folder = self.inputs.restart_folder
        computer = self.node.computer
        included = ['CHGCAR', 'WAVECAR']
        included += self.inputs.settings.base.attributes.get('ADDITIONAL_REMOTE_COPY_LIST', [])
        existing_objects = restart_folder.listdir()
        to_copy = []
        for name in included:
            if name not in existing_objects:
                # Here we simple issue an warning as the requirement of objects will be explicitly checked by
                # `write_additional` method
                self.report(f'WARNING: Object {name} does not exist in the restart folder.')
            else:
                to_copy.append(name)
        copy_list = [(computer.uuid, os.path.join(restart_folder.get_remote_path(), name), '.') for name in to_copy]
        return copy_list

    def verify_inputs(self) -> bool:
        """
        Hook to be extended by subclasses with checks for input nodes.

        Is called once before submission.
        """
        self.check_restart_folder()
        return True

    def check_restart_folder(self) -> None:
        restart_folder = self.inputs.get('restart_folder', None)
        if restart_folder:
            if not self.node.computer.pk == restart_folder.computer.pk:
                raise ValidationError('Calculation can not be restarted on another computer')

    def _is_restart(self) -> bool:
        restart_folder = self.inputs.get('restart_folder', None)
        return bool(restart_folder)

    def store(self, *args: Any, **kwargs: Any) -> None:
        """Adds a _prestore subclass hook for operations that should be done just before storing."""
        self._prestore()
        super().store(*args, **kwargs)

    def _prestore(self) -> None:
        """Subclass hook for updating attributes etc, just before storing."""
        return

    def write_additional(self, folder: Folder, calcinfo: CalcInfo) -> None:  # pylint: disable=unused-argument,
        """Subclass hook to write additional input objects."""
        return
