# pylint: disable=abstract-method
# explanation: pylint wrongly complains about (aiida) Node not implementing query
"""VaspImmigrant calculation: Immigrate externally run VASP calculations into AiiDA."""
from py import path as py_path  # pylint: disable=no-name-in-module,no-member
from aiida.common import InputValidationError
from aiida.common.lang import override
from aiida.engine import JobCalc

from aiida_vasp.data.potcar import PotcarData
from aiida_vasp.parsers.file_parsers.incar import IncarParser
from aiida_vasp.parsers.file_parsers.kpoints import KpointsParser
from aiida_vasp.parsers.file_parsers.poscar import PoscarParser
from aiida_vasp.parsers.file_parsers.potcar import MultiPotcarIo
from aiida_vasp.parsers.file_parsers.chgcar import ChgcarParser
from aiida_vasp.parsers.file_parsers.wavecar import WavecarParser
from aiida_vasp.utils.aiida_utils import get_data_node


class VaspImmigrantJobProcess(JobCalc):
    """
    JobCalc subclass for importing non-aiida VASP runs.

    Simulate running the VaspCalculation up to the point where it can be retrieved and parsed,
    then hand over control to the runner for the rest.
    """

    @override
    def run(self):
        import plumpy
        from aiida.common.datastructures import calc_states
        from aiida.common.links import LinkType
        from aiida.engine.job_processes import RETRIEVE_COMMAND
        from aiida.engine.calculation.job import _input_subfolder

        _ = super(VaspImmigrantJobProcess, self).run()

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
            if 'wavefunctions' in self.inputs:
                transport.get(remote_path.join('WAVECAR').strpath, raw_input_folder.abspath)
            if 'charge_density' in self.inputs:
                transport.get(remote_path.join('CHGCAR').strpath, raw_input_folder.abspath)

        self.calc._set_state(calc_states.COMPUTED)  # pylint: disable=protected-access
        return plumpy.Wait(msg='Waiting to retrieve', data=RETRIEVE_COMMAND)


def get_incar_input(dir_path):
    incar = IncarParser(file_path=dir_path.join('INCAR').strpath).incar
    return get_data_node('dict', dict=incar)


def get_poscar_input(dir_path):
    return PoscarParser(file_path=dir_path.join('POSCAR').strpath).structure


def get_potcar_input(dir_path, structure=None, potcar_spec=None):
    """Read potentials from a POTCAR file or POSCAR/structure plus potcar spec {'family': ..., 'map', ...}."""
    local_potcar = dir_path.join('POTCAR')
    structure = structure or get_poscar_input(dir_path)
    potentials = {}
    if local_potcar.exists():
        potentials = MultiPotcarIo.read(local_potcar.strpath).get_potentials_dict(structure)
        potentials = {(kind,): potentials[kind] for kind in potentials}
    elif potcar_spec:
        potentials = PotcarData.get_potcars_from_structure(structure, potcar_spec['family'], potcar_spec['map'])
    else:
        raise InputValidationError('no POTCAR found in remote folder and potcar_spec was not passed')

    return potentials


def get_kpoints_input(dir_path, structure=None):
    structure = structure or get_poscar_input(dir_path)
    kpoints = KpointsParser(file_path=dir_path.join('KPOINTS').strpath).kpoints
    kpoints.set_cell_from_structure(structure)
    return kpoints


def get_chgcar_input(dir_path):
    return ChgcarParser(file_path=dir_path.join('CHGCAR').strpath).chgcar


def get_wavecar_input(dir_path):
    return WavecarParser(file_path=dir_path.join('WAVECAR').strpath).wavecar
