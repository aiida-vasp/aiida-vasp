"""
Base calculation class.

-----------------------
Base and meta classes for VASP calculation classes.
"""
# pylint: disable=abstract-method,invalid-metaclass,ungrouped-imports
# explanation: pylint wrongly complains about Node not implementing query
import os

from aiida.engine import CalcJob
from aiida.common import CalcInfo, CodeInfo, ValidationError

from aiida_vasp.utils.aiida_utils import get_data_class


class VaspCalcBase(CalcJob):
    """
    Base class of all calculations utilizing VASP.

    * Defines internal parameters common to all vasp calculations.
    * provides a basic, extendable implementation of _prepare_for_submission
    * provides hooks, so subclasses can extend the behavior without having to reimplement common functionality
    """

    _default_parser = 'vasp.vasp'

    @classmethod
    def define(cls, spec):
        super(VaspCalcBase, cls).define(spec)
        spec.input(
            'restart_folder',
            valid_type=get_data_class('core.remote'),
            help='A remote folder to restart from if need be',
            required=False,
        )

    @classmethod
    def max_retrieve_list(cls):
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
            'OUTCAR'
        ]
        return retrieve_list

    def prepare_for_submission(self, folder):  # pylint: disable=arguments-differ
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

    def remote_copy_restart_folder(self):
        """Add all objects required for restart to the list of objects to be copied from the previous calculation."""
        restart_folder = self.inputs.restart_folder
        computer = self.node.computer
        included = ['CHGCAR', 'WAVECAR']
        existing_objects = restart_folder.listdir()
        to_copy = []
        for name in included:
            if name not in existing_objects:
                # Here we simple issue an warning as the requirement of objects will be explicitly checked by
                # `write_additional` method
                self.report('WARNING: Object {} does not exist in the restart folder.'.format(name))
            else:
                to_copy.append(name)
        copy_list = [(computer.uuid, os.path.join(restart_folder.get_remote_path(), name), '.') for name in to_copy]
        return copy_list

    def verify_inputs(self):
        """
        Hook to be extended by subclasses with checks for input nodes.

        Is called once before submission.
        """
        self.check_restart_folder()
        return True

    def check_restart_folder(self):
        restart_folder = self.inputs.get('restart_folder', None)
        if restart_folder:
            if not self.node.computer.pk == restart_folder.computer.pk:
                raise ValidationError('Calculation can not be restarted on another computer')

    def _is_restart(self):
        restart_folder = self.inputs.get('restart_folder', None)
        is_restart = False
        if restart_folder:
            is_restart = True
        return is_restart

    def store(self, *args, **kwargs):
        """Adds a _prestore subclass hook for operations that should be done just before storing."""
        self._prestore()
        super().store(*args, **kwargs)

    def _prestore(self):  # pylint: disable=no-self-use
        """Subclass hook for updating attributes etc, just before storing."""
        return

    def write_additional(self, folder, calcinfo):  # pylint: disable=no-self-use, unused-argument,
        """Subclass hook to write additional input objects."""
        return
