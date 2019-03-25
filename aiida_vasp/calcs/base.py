# pylint: disable=abstract-method,invalid-metaclass,ungrouped-imports
# explanation: pylint wrongly complains about Node not implementing query
"""Base and meta classes for VASP calculations"""
import os
import six
from py import path as py_path  # pylint: disable=no-name-in-module,no-member

from aiida.engine import CalcJob
from aiida.plugins import DataFactory
from aiida.common import CalcInfo, CodeInfo, ValidationError
from aiida.common.folders import SandboxFolder

from aiida_vasp.utils.aiida_utils import get_data_class, get_data_node, cmp_get_transport


class VaspCalcBase(CalcJob):
    """
    Base class of all calculations utilizing VASP.

    * Defines internal parameters common to all vasp calculations.
    * provides a basic, extendable implementation of _prepare_for_submission
    * provides hooks, so subclasses can extend the behaviour without

    having to reimplement common functionality
    """

    _INPUT_FILE_NAME = 'INCAR'
    _OUTPUT_FILE_NAME = 'OUTCAR'
    _default_parser = 'vasp.vasp'

    @classmethod
    def define(cls, spec):
        super(VaspCalcBase, cls).define(spec)
        spec.input('restart_folder', valid_type=get_data_class('remote'), help='A remote folder to restart from if need be', required=False)

    @classmethod
    def max_retrieve_list(cls):
        """Return a list of all possible output files from a VASP run."""
        retrieve_list = [
            'CHG', 'CHGCAR', 'CONTCAR', 'DOSCAR', 'EIGENVAL', 'ELFCAR', 'IBZKPT', 'LOCPOT', 'OSZICAR', 'OUTCAR', 'PCDAT', 'PROCAR',
            'PROOUT', 'STOPCAR', 'TMPCAR', 'WAVECAR', 'XDATCAR', ['wannier90*', '.', 0], 'vasprun.xml'
        ]
        return retrieve_list

    def prepare_for_submission(self, tempfolder):
        """
        Writes the four minimum output files, INCAR, POSCAR, POTCAR, KPOINTS.

        Delegates the construction and writing / copying to write_<file> methods.
        That way, subclasses can use any form of input nodes and just
        have to implement the write_xxx method accordingly.
        Subclasses can extend by calling the super method and if neccessary
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
        codeinfo.code_uuid = self.get_code().uuid
        codeinfo.code_pk = self.get_code().pk
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
        computer = self.get_computer()
        excluded = ['INCAR', '_aiidasubmit.sh', '.aiida']
        copy_list = [(computer.uuid, os.path.join(restart_folder.get_remote_path(), name), '.')
                     for name in restart_folder.listdir()
                     if name not in excluded]
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
            previous_calc = restart_folder.get_inputs(node_type=CalcJob)[0]
            if not self.get_computer().pk == previous_calc.get_computer().pk:
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

    def _prestore(self):
        """Subclass hook for updating attributes etc, just before storing"""
        pass

    def write_additional(self, tempfolder, calcinfo):
        """Subclass hook to write additional input files."""
        pass

    @classmethod
    def immigrant(cls, code, remote_path, **kwargs):

        return NotImplemented('The immigrant is not yet ported to comply with AiiDA beta. In fact, '
                              'we will most likely wait until an immigrant function is present '
                              'in AiiDA core.')
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
        from aiida_vasp.calcs import immigrant as imgr
        remote_path = py_path.local(remote_path)
        proc_cls = imgr.VaspImmigrantJobProcess.build(cls)
        builder = proc_cls.get_builder()
        builder.code = code
        builder.options.max_wallclock_seconds = 1  # pylint: disable=no-member
        builder.options.resources = kwargs.get('resources', {'num_machines': 1, 'num_mpiprocs_per_machine': 1})  # pylint: disable=no-member
        settings = kwargs.get('settings', {})
        settings.update({'import_from_path': remote_path.strpath})
        builder.settings = get_data_node('dict', dict=settings)
        with cmp_get_transport(code.get_computer()) as transport:
            with SandboxFolder() as sandbox:
                sandbox_path = py_path.local(sandbox.abspath)
                transport.get(remote_path.join('INCAR').strpath, sandbox_path.strpath)
                transport.get(remote_path.join('POSCAR').strpath, sandbox_path.strpath)
                transport.get(remote_path.join('POTCAR').strpath, sandbox_path.strpath, ignore_nonexisting=True)
                transport.get(remote_path.join('KPOINTS').strpath, sandbox_path.strpath)
                builder.parameters = imgr.get_incar_input(sandbox_path)
                builder.structure = imgr.get_poscar_input(sandbox_path)
                builder.potential = imgr.get_potcar_input(sandbox_path, potcar_spec=kwargs.get('potcar_spec', None))
                builder.kpoints = imgr.get_kpoints_input(sandbox_path)
                cls._immigrant_add_inputs(transport, remote_path=remote_path, sandbox_path=sandbox_path, builder=builder, **kwargs)
        return proc_cls, builder

    @classmethod
    def _immigrant_add_inputs(cls, transport, remote_path, sandbox_path, builder, **kwargs):
        pass
