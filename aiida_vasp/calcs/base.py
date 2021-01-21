"""
Base calculation class.

-----------------------
Base and meta classes for VASP calculation classes.
"""
# pylint: disable=abstract-method,invalid-metaclass,ungrouped-imports
# explanation: pylint wrongly complains about Node not implementing query
import os
from pathlib import Path

from aiida.engine import CalcJob
from aiida.common import CalcInfo, CodeInfo, ValidationError
from aiida.common.folders import SandboxFolder

from aiida_vasp.utils.aiida_utils import get_data_class, get_data_node, cmp_get_transport


class VaspCalcBase(CalcJob):
    """
    Base class of all calculations utilizing VASP.

    * Defines internal parameters common to all vasp calculations.
    * provides a basic, extendable implementation of _prepare_for_submission
    * provides hooks, so subclasses can extend the behaviour without having to reimplement common functionality
    """

    _default_parser = 'vasp.vasp'

    @classmethod
    def define(cls, spec):
        super(VaspCalcBase, cls).define(spec)
        spec.input('restart_folder', valid_type=get_data_class('remote'), help='A remote folder to restart from if need be', required=False)

    @classmethod
    def max_retrieve_list(cls):
        """Return a list of all possible output files from a VASP run."""
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
            'STOPCAR',  # Controlled stopping file
            'vasprun.xml',
            'OUTCAR'
        ]
        return retrieve_list

    def prepare_for_submission(self, tempfolder):  # pylint: disable=arguments-differ
        """
        Writes the four minimum output files, INCAR, POSCAR, POTCAR, KPOINTS.

        Delegates the construction and writing / copying to write_<file> methods.
        That way, subclasses can use any form of input nodes and just
        have to implement the write_xxx method accordingly.
        Subclasses can extend by calling the super method and if necessary
        modifying it's output CalcInfo before returning it.
        """
        # write input files
        incar = tempfolder.get_abs_path('INCAR')
        structure = tempfolder.get_abs_path('POSCAR')
        potentials = tempfolder.get_abs_path('POTCAR')
        kpoints = tempfolder.get_abs_path('KPOINTS')

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
        self.write_additional(tempfolder, calcinfo)

        return calcinfo

    def remote_copy_restart_folder(self):
        """Add all files required for restart to the list of files to be copied from the previous calculation."""
        restart_folder = self.inputs.restart_folder
        computer = self.node.computer
        included = ['CHGCAR', 'WAVECAR']
        existing_files = restart_folder.listdir()
        to_copy = []
        for name in included:
            if name not in existing_files:
                # Here we simple issue an warning as the requirement of files will be explicitly checked by
                # `write_additional` method
                self.report('WARNING: File {} does not exist in the restart folder.'.format(name))
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
        super(VaspCalcBase, self).store(*args, **kwargs)

    def _prestore(self):  # pylint: disable=no-self-use
        """Subclass hook for updating attributes etc, just before storing."""
        return

    def write_additional(self, tempfolder, calcinfo):  # pylint: disable=no-self-use, unused-argument,
        """Subclass hook to write additional input files."""
        return

    @classmethod
    def immigrant(cls, code, remote_path, **kwargs):
        """
        Create an immigrant with appropriate inputs from a code and a remote path on the associated computer.

        More inputs are required to pass resources information, if the POTCAR file is missing from the folder
        or if additional settings need to be passed, e.g. parser instructions.

        :param code: a Code instance for the code originally used.
        :param remote_path: The directory on the code's computer in which the simulation was run.
        :param resources: dict. The resources used during the run (defaults to 1 machine, 1 process).
        :param potcar_spec: dict. If the POTCAR file is not present anymore, this allows to pass a family and
            mapping to find the right POTCARs.
        :param settings: dict. Used for non-default parsing instructions, etc.
        """

        from aiida_vasp.calcs import immigrant as imgr  # pylint: disable=import-outside-toplevel
        remote_path = Path(remote_path)
        proc_cls = imgr.VaspImmigrant
        builder = proc_cls.get_builder()
        builder.code = code
        options = {'max_wallclock_seconds': 1, 'resources': {'num_machines': 1, 'num_mpiprocs_per_machine': 1}}
        metadata = kwargs.get('metadata', {'options': options})
        options = metadata.get('options', options)
        max_wallclock_seconds = options.get('max_wallclock_seconds', 1)
        resources = options.get('resources', {'num_machines': 1, 'num_mpiprocs_per_machine': 1})
        builder.metadata['options']['max_wallclock_seconds'] = max_wallclock_seconds  # pylint: disable=no-member
        builder.metadata['options']['resources'] = resources  # pylint: disable=no-member
        settings = kwargs.get('settings', {})
        settings.update({'import_from_path': str(remote_path)})
        builder.settings = get_data_node('dict', dict=settings)
        with cmp_get_transport(code.computer) as transport:
            with SandboxFolder() as sandbox:
                sandbox_path = Path(sandbox.abspath)
                transport.get(str(remote_path / 'INCAR'), str(sandbox_path))
                transport.get(str(remote_path / 'POSCAR'), str(sandbox_path))
                transport.get(str(remote_path / 'POTCAR'), str(sandbox_path), ignore_nonexisting=True)
                transport.get(str(remote_path / 'KPOINTS'), str(sandbox_path))
                builder.parameters = imgr.get_incar_input(sandbox_path)
                builder.structure = imgr.get_poscar_input(sandbox_path)
                builder.potential = imgr.get_potcar_input(sandbox_path,
                                                          potential_family=kwargs.get('potential_family'),
                                                          potential_mapping=kwargs.get('potential_mapping'))
                builder.kpoints = imgr.get_kpoints_input(sandbox_path)
                cls._immigrant_add_inputs(transport, remote_path=remote_path, sandbox_path=sandbox_path, builder=builder, **kwargs)
        return proc_cls, builder

    @classmethod
    def _immigrant_add_inputs(cls, transport, remote_path, sandbox_path, builder, **kwargs):
        pass
