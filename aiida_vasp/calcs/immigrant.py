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
from aiida_vasp.utils.aiida_utils import get_data_node


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


# def get_incar_input(dir_path):
#     incar_dict = IncarParser(file_path=dir_path.join('INCAR').strpath).incar
#     return get_data_node('parameter', dict=incar_dict)

# def get_poscar_input(dir_path):
#     return PoscarParser(file_path=dir_path.join('POSCAR').strpath).structure

# def get_immigrant_builder(code, remote_path, resources=None):
#     builder = VaspImmigrantJobProcess.get_builder()
#     builder.code = code
#     builder.options.resources = resources or {'num_machines': 1, 'num_mpiprocs_per_machine': 1}
#     builder.settings = get_data_node('parameter', dict={'import_from_path': remote_path})
#     with code.get_computer().get_transport() as transport:
#         with SandboxFolder() as sandbox:
#             sandbox_path = py_path.local(sandbox.abspath)
#             transport.get(remote_path.join('INCAR').strpath, sandbox_path.strpath)
#             transport.get(remote_path.join('POSCAR').strpath, sandbox_path.strpath)
#             builder.parameters = get_incar_input(sandbox_path)
#             builder.structure = get_poscar_input(sandbox_path)
#     return builder


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
                Such an instance can be obtained from aiida_vasp.utils.aiida_utils.cmp_get_transport().
        :param remote_workdir: optional only if ``remote_workdir`` has already been set previously using ``set_remote_workdir()``
                Absolute path to the directory where the calculation was run. All input and output files are expected to be found
                at the given path by the given transport instance.
        :param settings_dict: An optional dictionary with a dictionary for the ``settings`` input.
        :param potcar_spec: A dictionary {'family': <potcar family name>, 'map': <dict: kind_name -> potcar full_name>}

        Example (requires aiida v1.0.0a2+)::

            localhost = Computer.get('localhost')  # this computer
            immigrant = VaspImmigrant(computer=localhost, resources={...})
            remote_workdir = '/<path-to-repo>/aiida_vasp/test_data/phonondb'
            settings_dict = {'parser_settings': {'add_structure': True}}
            potcar_spec = {'family': 'PBE54', map={'P': 'P', 'S': 'S', 'Zn': 'Zn'}}
            with localhost.get_transport() as open_transport:
                immigrant.create_input_nodes(open_transport, remote_workdir, settings_dict, potcar_spec)
            immigrant.get_inputs_dict()
            ## outputs
            ## {'parameters': <ParameterData ...>, 'structure': <...>, ...}
        """
        self.test_open_transport(open_transport)
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

    def prepare_for_retrieval_and_parsing(self, open_transport):  # pylint:disable=invalid-name
        """
        Emulate running the calculation setting its state for parsing.

        :param open_transport: An open instance of the transport used by the calculation's computer.
                Such an instance can be obtained from aiida_vasp.utils.aiida_utils.cmp_get_transport().

        Example (requires aiida v1.0.0a2+)::

            comp = Computer.get(...)
            immigrant = VaspImmigrant(remote_folder=..., computer=comp, resources={...})
            with comp.get_transport() as open_transport:
                immigrant.create_input_nodes(open_transport, ...)
                immigrant.prepare_for_retrieval_and_parsing(open_transport)
        """
        from aiida.common.datastructures import calc_states
        from aiida.common.links import LinkType
        from aiida.orm.calculation.job import _input_subfolder

        self.test_open_transport(open_transport)

        self._set_attr('retrieve_list', self.max_retrieve_list())
        self.store_all()

        remote_path = py_path.local(self._get_remote_workdir())
        raw_input_folder = self.folder.get_subfolder(_input_subfolder, create=True)
        open_transport.get(remote_path.join('INCAR').strpath, raw_input_folder.abspath)
        open_transport.get(remote_path.join('KPOINTS').strpath, raw_input_folder.abspath)
        open_transport.get(remote_path.join('POSCAR').strpath, raw_input_folder.abspath)

        self._set_state(calc_states.SUBMITTING)

        remotedata = get_data_node('remote', computer=self.get_computer(), remote_path=self._get_remote_workdir())
        remotedata.add_link_from(self, label='remote_folder', link_type=LinkType.CREATE)
        remotedata.store()

        self._set_state(calc_states.COMPUTED)

    def test_open_transport(self, open_transport):
        """Make sure the transport is open and fits the computer."""
        if not self.get_computer():
            raise InputValidationError('The computer must be set first!')
        if not isinstance(open_transport, self.get_computer().get_transport_class()):
            raise InputValidationError('The transport passed does not match the transport type of the computer!')
        if not open_transport._is_open:  # pylint: disable=protected-access
            raise InputValidationError('The transport passed is not open!')
