"""
Immigrant calculation.

----------------------
Enables the immigration of  externally run VASP calculations into AiiDA.
"""
# pylint: disable=abstract-method, import-outside-toplevel, cyclic-import
# explanation: pylint wrongly complains about (aiida) Node not implementing query
from pathlib import Path

from aiida.common import InputValidationError
from aiida.common.lang import override
from aiida.common.links import LinkType
from aiida.common.folders import SandboxFolder
from aiida.common.extendeddicts import AttributeDict

from aiida_vasp.calcs.vasp import VaspCalculation
from aiida_vasp.data.potcar import PotcarData
from aiida_vasp.parsers.file_parsers.incar import IncarParser
from aiida_vasp.parsers.file_parsers.kpoints import KpointsParser
from aiida_vasp.parsers.file_parsers.poscar import PoscarParser
from aiida_vasp.parsers.file_parsers.potcar import MultiPotcarIo
from aiida_vasp.parsers.file_parsers.chgcar import ChgcarParser
from aiida_vasp.parsers.file_parsers.wavecar import WavecarParser
from aiida_vasp.utils.aiida_utils import get_data_node
from aiida_vasp.utils.aiida_utils import cmp_get_transport

# _IMMIGRANT_EXTRA_KWARGS = """
# vasp.vasp specific kwargs:

# :param use_chgcar: bool, if True, read the CHGCAR file (has to exist) and convert it to an input node.
# :param use_wavecar: bool, if True, read the WAVECAR file (has to exist) and convert it to an input node.
# """

# @update_docstring('immigrant', _IMMIGRANT_EXTRA_KWARGS, append=True)


class VaspImmigrant(VaspCalculation):
    """Parse VASP output files stored in a specified directory.

    Simulate running the VaspCalculation up to the point where it can be
    retrieved and parsed, then hand over control to the runner for the rest.

    Usage examples
    --------------
    Immigrant calculation can be perfomed as follows.

    ::

       code = Code.get_from_string('vasp@local')
       folder = '/home/username/vasp-calc-dir'
       settings = {'parser_settings': {'add_energies': True,
                                       'add_forces': True,
                                       'electronic_step_energies': True}}
       VaspImmigrant = CalculationFactory('vasp.immigrant')
       builder = VaspImmigrant.get_builder_from_folder(code,
                                                       folder,
                                                       settings=settings)
       submit(builder)

    Instead of ``builder``, inputs dict is obtained similarly as

    ::

       code = Code.get_from_string('vasp@local')
       folder = '/home/username/vasp-calc-dir'
       settings = {'parser_settings': {'add_energies': True,
                                       'add_forces': True,
                                       'electronic_step_energies': True}}
       VaspImmigrant = CalculationFactory('vasp.immigrant')
       inputs = VaspImmigrant.get_inputs_from_folder(code,
                                                     folder,
                                                     settings=settings)
       submit(VaspImmigrant, **inputs)

    Note
    ----
    The defaul metadata is set automatically as follows::

       {'options': {'max_wallclock_seconds': 1,
        'resources': {'num_machines': 1, 'num_mpiprocs_per_machine': 1}}}

    Specific scheduler may require setting ``resources`` differently
    (e.g., sge ``'parallel_env'``).

    ``get_inputs_from_folder`` and ``get_builder_from_folder`` accept several
    kwargs, see the docstring of ``get_inputs_from_folder``.

    """

    @override
    def run(self):
        import plumpy
        from aiida.engine.processes.calcjobs.tasks import RETRIEVE_COMMAND

        _ = super(VaspImmigrant, self).run()

        # Make sure the retrieve list is set (done in presubmit so we need to call that also)
        with SandboxFolder() as folder:
            self.presubmit(folder)

        settings = self.inputs.get('settings', None)
        settings = settings.get_dict() if settings else {}
        remote_path = settings.get('import_from_path', None)
        if not remote_path:
            raise InputValidationError('immigrant calculations need an input "settings" containing a key "import_from_path"!')
        self.node.set_remote_workdir(remote_path)  # pylint: disable=protected-access
        remotedata = get_data_node('remote', computer=self.node.computer, remote_path=remote_path)
        remotedata.add_incoming(self.node, link_type=LinkType.CREATE, link_label='remote_folder')
        remotedata.store()

        return plumpy.Wait(msg='Waiting to retrieve', data=RETRIEVE_COMMAND)

    @classmethod
    def get_inputs_from_folder(cls, code, remote_path, **kwargs):
        """
        Create inputs to launch immigrant from a code and a remote path on the associated computer.

        If POTCAR does not exist, the provided ``potential_family`` and
        ``potential_mapping`` are used to link potential to inputs. In this
        case, at least ``potential_family`` has to be provided. Unless
        ``potential_mapping``, this mapping is generated from structure, i.e.,

        ::

            potential_mapping = {element: element for element in structure.get_kind_names()}

        :param code: a Code instance for the code originally used.
        :param remote_path: Directory or folder name where VASP inputs and outputs are stored.
        :param settings: dict. This is used as the input port of VaspCalculation.
        :param potential_family: str. This will be obsolete at v3.0.
        :param potential_mapping: dict. This will be obsolete at v3.0.
        :param use_wavecar: bool. Try to read WAVECAR.
        :param use_chgcar bool. Try to read CHGCAR.
        """

        inputs = AttributeDict()
        inputs.code = code
        options = {'max_wallclock_seconds': 1, 'resources': {'num_machines': 1, 'num_mpiprocs_per_machine': 1}}
        inputs.metadata = AttributeDict()
        inputs.metadata.options = options
        settings = kwargs.get('settings', {})
        settings.update({'import_from_path': str(remote_path)})
        inputs.settings = get_data_node('dict', dict=settings)
        _remote_path = Path(remote_path)
        with cmp_get_transport(code.computer) as transport:
            with SandboxFolder() as sandbox:
                sandbox_path = Path(sandbox.abspath)
                transport.get(str(_remote_path / 'INCAR'), str(sandbox_path))
                transport.get(str(_remote_path / 'POSCAR'), str(sandbox_path))
                transport.get(str(_remote_path / 'POTCAR'), str(sandbox_path), ignore_nonexisting=True)
                transport.get(str(_remote_path / 'KPOINTS'), str(sandbox_path))
                inputs.parameters = get_incar_input(sandbox_path)
                inputs.structure = get_poscar_input(sandbox_path)
                try:
                    inputs.potential = get_potcar_input(sandbox_path,
                                                        potential_family=kwargs.get('potential_family'),
                                                        potential_mapping=kwargs.get('potential_mapping'))
                except InputValidationError:
                    pass
                inputs.kpoints = get_kpoints_input(sandbox_path)
                cls._add_inputs(transport, _remote_path, sandbox_path, inputs, **kwargs)

        return inputs

    @classmethod
    def get_builder_from_folder(cls, code, remote_path, **kwargs):
        """
        Create an immigrant builder from a code and a remote path on the associated computer.
        See more details in the docstring of ``get_inputs_from_folder``.
        """

        inputs = cls.get_inputs_from_folder(code, remote_path, **kwargs)
        builder = cls.get_builder()
        for key, val in inputs.items():
            builder[key] = val
        return builder

    @classmethod
    def _add_inputs(cls, transport, remote_path, sandbox_path, inputs, **kwargs):
        """Add some more inputs"""
        add_wavecar = kwargs.get('use_wavecar') or bool(inputs.parameters.get_dict().get('istart', 0))
        add_chgcar = kwargs.get('use_chgcar') or inputs.parameters.get_dict().get('icharg', -1) in [1, 11]
        if add_chgcar:
            transport.get(str(remote_path / 'CHGCAR'), str(sandbox_path))
            inputs.charge_density = get_chgcar_input(sandbox_path)
        if add_wavecar:
            transport.get(str(remote_path / 'WAVECAR'), str(sandbox_path))
            inputs.wavefunctions = get_wavecar_input(sandbox_path)


def get_incar_input(dir_path):
    incar = IncarParser(file_path=str(dir_path / 'INCAR')).incar
    return get_data_node('dict', dict=incar)


def get_poscar_input(dir_path):
    return PoscarParser(file_path=str(dir_path / 'POSCAR')).structure


def get_potcar_input(dir_path, structure=None, potential_family=None, potential_mapping=None):
    """Read potentials from a POTCAR file or set it up from a structure."""
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


def get_kpoints_input(dir_path, structure=None):
    structure = structure or get_poscar_input(dir_path)
    kpoints = KpointsParser(file_path=str(dir_path / 'KPOINTS')).kpoints
    kpoints.set_cell_from_structure(structure)
    return kpoints


def get_chgcar_input(dir_path):
    return ChgcarParser(file_path=str(dir_path / 'CHGCAR')).chgcar


def get_wavecar_input(dir_path):
    return WavecarParser(file_path=str(dir_path / 'WAVECAR')).wavecar
