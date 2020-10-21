"""
Immigrant calculation.

----------------------
Enables the immigration of  externally run VASP calculations into AiiDA.
"""
# pylint: disable=abstract-method, import-outside-toplevel, cyclic-import
# explanation: pylint wrongly complains about (aiida) Node not implementing query
from aiida.common import InputValidationError
from aiida.common.lang import override
from aiida.common.links import LinkType

from aiida_vasp.calcs.vasp import VaspCalculation
from aiida_vasp.data.potcar import PotcarData
from aiida_vasp.parsers.file_parsers.incar import IncarParser
from aiida_vasp.parsers.file_parsers.kpoints import KpointsParser
from aiida_vasp.parsers.file_parsers.poscar import PoscarParser
from aiida_vasp.parsers.file_parsers.potcar import MultiPotcarIo
from aiida_vasp.parsers.file_parsers.chgcar import ChgcarParser
from aiida_vasp.parsers.file_parsers.wavecar import WavecarParser
from aiida_vasp.utils.aiida_utils import get_data_node


class VaspImmigrant(VaspCalculation):
    """
    CalcJob subclass for importing non-aiida VASP runs.

    Simulate running the VaspCalculation up to the point where it can be retrieved and parsed,
    then hand over control to the runner for the rest.
    """

    @override
    def run(self):
        import plumpy
        from aiida.engine.processes.calcjobs.tasks import RETRIEVE_COMMAND
        from aiida.common.folders import SandboxFolder

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
