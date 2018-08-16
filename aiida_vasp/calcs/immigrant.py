# pylint: disable=abstract-method
# explanation: pylint wrongly complains about (aiida) Node not implementing query
"""VaspImmigrant calculation: Immigrate externally run VASP calculations into AiiDA."""
import os

from py import path as py_path  # pylint: disable=no-name-in-module,no-member
from aiida.common.folders import SandboxFolder
from aiida.common.exceptions import InputValidationError
from aiida.common.lang import override
from aiida.work import JobProcess

from aiida_vasp.calcs.vasp import VaspCalculation
from aiida_vasp.data.potcar import PotcarData
from aiida_vasp.io.incar import IncarParser
from aiida_vasp.io.kpoints import KpParser
from aiida_vasp.io.poscar import PoscarParser
from aiida_vasp.io.potcar import MultiPotcarIo
from aiida_vasp.utils.aiida_utils import get_data_node, cmp_get_transport


def get_quantity_node(parser, quantity):
    return parser.get_quantity(quantity, None)[quantity]


class VaspImmigrantJobProcess(JobProcess):
    """
    JobProcess subclass for importing non-aiida VASP runs.

    Simulate running the VaspCalculation up to the point where it can be retrieved and parsed,
    then hand over control to the runner for the rest.
    """

    @override
    def run(self):
        import plumpy
        from aiida.common.datastructures import calc_states
        from aiida.common.links import LinkType
        from aiida.work.job_processes import RETRIEVE_COMMAND
        from aiida.orm.calculation.job import _input_subfolder

        _ = super(VaspImmigrantJobProcess, self).run()

        state = self.calc.get_state()

        if state != calc_states.TOSUBMIT:
            raise RuntimeError('immigrant calculations should always begin in NEW state')

        def return_empty_list():
            return []

        setattr(self.calc, '_get_retrieve_list', self.calc.max_retrieve_list)
        setattr(self.calc, '_get_retrieve_singlefile_list', return_empty_list)
        setattr(self.calc, '_get_retrieve_temporary_list', return_empty_list)

        settings = self.calc.get_inputs_dict().get('settings', None)
        settings = settings.get_dict() if settings else {}
        remote_path = settings.get('import_from_path', None)
        if not remote_path:
            raise InputValidationError('immigrant calculations need an input "settings" containing a key "import_from_path"!')
        self.calc._set_state(calc_states.SUBMITTING)  # pylint: disable=protected-access
        self.calc._set_attr('remote_workdir', remote_path)  # pylint: disable=protected-access
        remotedata = get_data_node('remote', computer=self.calc.get_computer(), remote_path=remote_path)
        remotedata.add_link_from(self.calc, label='remote_folder', link_type=LinkType.CREATE)
        remotedata.store()

        remote_path = py_path.local(remote_path)
        with self.calc.get_computer().get_transport() as transport:
            raw_input_folder = self.calc.folder.get_subfolder(_input_subfolder, create=True)
            transport.get(remote_path.join('INCAR').strpath, raw_input_folder.abspath)
            transport.get(remote_path.join('KPOINTS').strpath, raw_input_folder.abspath)
            transport.get(remote_path.join('POSCAR').strpath, raw_input_folder.abspath)

        self.calc._set_state(calc_states.COMPUTED)  # pylint: disable=protected-access
        return plumpy.Wait(msg='Waiting to retrieve', data=RETRIEVE_COMMAND)


def get_incar_input(dir_path):
    return get_quantity_node(IncarParser(file_path=dir_path.join('INCAR').strpath), 'incar')


def get_poscar_input(dir_path):
    return get_quantity_node(PoscarParser(file_path=dir_path.join('POSCAR').strpath), 'poscar-structure')


def get_potcar_input(dir_path, structure=None, potcar_spec=None):
    """Read potentials from a POTCAR file or POSCAR/structure plus potcar spec {'family': ..., 'map', ...}."""
    local_potcar = dir_path.join('POTCAR')
    structure = structure or get_poscar_input(dir_path)
    potentials = {}
    if local_potcar.exists():
        potentials = MultiPotcarIo.read(local_potcar.strpath).get_potentials_dict(structure)
        potentials = {(kind,): potcar for kind, potcar in potentials}
    elif potcar_spec:
        potentials = PotcarData.get_potcars_from_structure(structure, potcar_spec['family'], potcar_spec['map'])
    else:
        raise InputValidationError('no POTCAR found in remote folder and potcar_spec was not passed')

    return potentials


def get_kpoints_input(dir_path, structure=None):
    structure = structure or get_poscar_input(dir_path)
    kpoints = get_quantity_node(KpParser(file_path=dir_path.join('KPOINTS').strpath), 'kpoints-kpoints')
    kpoints.set_cell_from_structure(structure)
    return kpoints


def get_immigrant_with_builder(code, remote_path, resources=None, potcar_spec=None, settings=None):
    """
    Create an immigrant with appropriate inputs from a code and a remote path on the associated computer.

    More inputs are required to pass resources information, if the POTCAR file is missing from the folder
    or if additional settings need to be passed, e.g. parser instructions.

    :param code: a Code instance for the code originally used.
    :param remote_path: The directory on the code's computer in which the simulation was run.
    :param resources: dict. The resources used during the run (defaults to 1 machine, 1 process).
    :param potcar_spec: dict. If the POTCAR file is not present anymore, this allows to pass a family and mapping to find the right POTCARs.
    :param settings: dict. Used for non-default parsing instructions, etc.
    """
    proc_cls = VaspImmigrantJobProcess.build(VaspCalculation)
    builder = proc_cls.get_builder()
    builder.code = code
    builder.options.resources = resources or {'num_machines': 1, 'num_mpiprocs_per_machine': 1}
    settings = settings or {}
    settings.update({'import_from_path': remote_path.strpath})
    builder.settings = get_data_node('parameter', dict=settings)
    with cmp_get_transport(code.get_computer()) as transport:
        with SandboxFolder() as sandbox:
            sandbox_path = py_path.local(sandbox.abspath)
            transport.get(remote_path.join('INCAR').strpath, sandbox_path.strpath)
            transport.get(remote_path.join('POSCAR').strpath, sandbox_path.strpath)
            transport.get(remote_path.join('POTCAR').strpath, sandbox_path.strpath, ignore_nonexisting=True)
            transport.get(remote_path.join('KPOINTS').strpath, sandbox_path.strpath)
            builder.parameters = get_incar_input(sandbox_path)
            builder.structure = get_poscar_input(sandbox_path)
            builder.potential = get_potcar_input(sandbox_path, potcar_spec=potcar_spec)
            builder.kpoints = get_kpoints_input(sandbox_path)
    return proc_cls, builder
