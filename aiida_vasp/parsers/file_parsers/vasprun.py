"""
vasprun parser.

---------------
The file parser that handles the parsing of vasprun.xml files.
"""
# pylint: disable=too-many-public-methods
import numpy as np

from parsevasp.vasprun import Xml
from parsevasp.kpoints import Kpoint
from parsevasp import constants as parsevaspct
from aiida_vasp.parsers.file_parsers.parser import BaseFileParser, SingleFile

DEFAULT_OPTIONS = {
    'quantities_to_parse': [
        'structure', 'eigenvalues', 'dos', 'bands', 'kpoints', 'occupancies', 'trajectory', 'energies', 'projectors', 'dielectrics',
        'born_charges', 'hessian', 'dynmat', 'forces', 'stress', 'total_energies', 'maximum_force', 'maximum_stress'
    ],
    'energy_type': ['energy_no_entropy']
}


class VasprunParser(BaseFileParser):
    """Interface to parsevasp's xml parser."""

    PARSABLE_ITEMS = {
        'structure': {
            'inputs': [],
            'name': 'structure',
            'prerequisites': [],
            'alternatives': ['poscar-structure']
        },
        'eigenvalues': {
            'inputs': [],
            'name': 'eigenvalues',
            'prerequisites': [],
            'alternatives': ['eigenval-eigenvalues']
        },
        'dos': {
            'inputs': [],
            'name': 'dos',
            'prerequisites': [],
            'alternatives': ['doscar-dos']
        },
        'kpoints': {
            'inputs': [],
            'name': 'kpoints',
            'prerequisites': [],
            'alternatives': ['kpoints-kpoints']
        },
        'occupancies': {
            'inputs': [],
            'name': 'occupancies',
            'prerequisites': [],
        },
        'trajectory': {
            'inputs': [],
            'name': 'trajectory',
            'prerequisites': [],
        },
        'energies': {
            'inputs': [],
            'name': 'energies',
            'prerequisites': [],
        },
        'total_energies': {
            'inputs': [],
            'name': 'total_energies',
            'prerequisites': [],
        },
        'projectors': {
            'inputs': [],
            'name': 'projectors',
            'prerequisites': [],
        },
        'dielectrics': {
            'inputs': [],
            'nodeName': 'dielectrics',
            'prerequisites': [],
        },
        'stress': {
            'inputs': [],
            'name': 'stress',
            'prerequisites': [],
        },
        'forces': {
            'inputs': [],
            'nodeName': 'forces',
            'prerequisites': [],
        },
        'born_charges': {
            'inputs': [],
            'name': 'born_charges',
            'prerequisites': [],
        },
        'hessian': {
            'inputs': [],
            'name': 'hessian',
            'prerequisites': [],
        },
        'dynmat': {
            'inputs': [],
            'name': 'dynmat',
            'prerequisites': [],
        },
        'fermi_level': {
            'inputs': [],
            'name': 'fermi_level',
            'prerequisites': [],
        },
        'maximum_force': {
            'inputs': [],
            'name': 'maximum_force',
            'prerequisites': []
        },
        'maximum_stress': {
            'inputs': [],
            'name': 'maximum_stress',
            'prerequisites': []
        },
    }

    def __init__(self, *args, **kwargs):
        super(VasprunParser, self).__init__(*args, **kwargs)
        self._xml = None
        self.init_with_kwargs(**kwargs)

    def _init_with_file_path(self, path):
        """Init with a filepath."""
        self._parsed_data = {}
        self.parsable_items = self.__class__.PARSABLE_ITEMS
        self._data_obj = SingleFile(path=path)

        # Since vasprun.xml can be fairly large, we will parse it only
        # once and store the parsevasp Xml object.
        try:
            self._xml = Xml(file_path=path, k_before_band=True, logger=self._logger)
        except SystemExit:
            self._logger.warning('Parsevasp exited abruptly. Returning None.')
            self._xml = None

    def _init_with_data(self, data):
        """Init with SingleFileData."""
        self.parsable_items = self.__class__.PARSABLE_ITEMS
        self._init_with_file_path(data.get_file_abs_path())

    def _parse_file(self, inputs):

        # Since all quantities will be returned by properties, we can't pass
        # inputs as a parameter, so we store them in self._parsed_data
        for key, value in inputs.items():
            self._parsed_data[key] = value

        quantities_to_parse = DEFAULT_OPTIONS.get('quantities_to_parse')
        if self.settings is not None and self.settings.quantities_to_parse:
            quantities_to_parse = self.settings.quantities_to_parse

        result = {}

        if self._xml is None:
            # parsevasp threw an exception, which means vasprun.xml could not be parsed.
            for quantity in quantities_to_parse:
                if quantity in self.parsable_items:
                    result[quantity] = None
            return result

        for quantity in quantities_to_parse:
            if quantity in self.parsable_items:
                result[quantity] = getattr(self, quantity)

        return result

    @property
    def eigenvalues(self):
        """Fetch eigenvalues from parsevasp."""

        # fetch eigenvalues
        eigenvalues = self._xml.get_eigenvalues()

        if eigenvalues is None:
            # eigenvalues not present
            self._vasp_parser.exit_status = self._vasp_parser.exit_codes.ERROR_NOT_ABLE_TO_PARSE_QUANTITY
            return None

        eigen = []
        eigen.append(eigenvalues.get('total'))

        if eigen[0] is None:
            # spin decomposed?
            eigen[0] = eigenvalues.get('up')
            eigen.append(eigenvalues.get('down'))

        if eigen[0] is None:
            # safety, should not really happen?
            self._vasp_parser.exit_status = self._vasp_parser.exit_codes.ERROR_NOT_ABLE_TO_PARSE_QUANTITY
            return None

        return eigen

    @property
    def occupancies(self):
        """Fetch occupancies from parsevasp."""

        # fetch occupancies
        occupancies = self._xml.get_occupancies()

        if occupancies is None:
            # occupancies not present, should not really happen?
            self._vasp_parser.exit_status = self._vasp_parser.exit_codes.ERROR_NOT_ABLE_TO_PARSE_QUANTITY
            return None

        occ = []
        occ.append(occupancies.get('total'))

        if occ[0] is None:
            # spin decomposed
            occ[0] = occupancies.get('up')
            occ.append(occupancies.get('down'))

        if occ[0] is None:
            # should not really happen
            self._vasp_parser.exit_status = self._vasp_parser.exit_codes.ERROR_NOT_ABLE_TO_PARSE_QUANTITY
            return None

        return occ

    @property
    def kpoints(self):
        """Fetch the kpoints from parsevasp an store in KpointsData."""

        kpts = self._xml.get_kpoints()
        kptsw = self._xml.get_kpointsw()
        kpoints_data = None
        if (kpts is not None) and (kptsw is not None):
            # create a KpointsData object and store k-points
            kpoints_data = {}
            kpoints_data['mode'] = 'explicit'
            kpoints_data['points'] = []
            for kpt, kptw in zip(kpts, kptsw):
                kpoints_data['points'].append(Kpoint(kpt, weight=kptw))

        return kpoints_data

    @property
    def structure(self):
        """
        Fetch a given structure.

        Which structure to fetch is controlled by inputs.

        eFL: Need to clean this so that we can set different
        structures to pull from the outside. Could be usefull not
        pulling the whole trajectory.

        Currently defaults to the last structure.

        """

        return self.last_structure

    @property
    def last_structure(self):
        """
        Fetch the structure.

        After or at the last recorded ionic step from parsevasp.

        """

        last_lattice = self._xml.get_lattice('final')
        if last_lattice is None:
            self._vasp_parser.exit_status = self._vasp_parser.exit_codes.ERROR_NOT_ABLE_TO_PARSE_QUANTITY
            return None
        return _build_structure(last_lattice)

    @property
    def final_structure(self):
        """
        Fetch the structure.

        After or at the last recorded ionic step from parsevasp. Should in
        principle be the same as the method above.

        """

        return self.last_structure

    @property
    def last_forces(self):
        """
        Fetch forces.

        After or at the last recorded ionic step from parsevasp.

        """

        force = self._xml.get_forces('final')
        return force

    @property
    def final_forces(self):
        """
        Fetch forces.

        After or at the last recorded ionic step from parsevasp.

        """

        return self.last_forces

    @property
    def forces(self):
        """
        Fetch forces.

        This container should contain all relevant forces.
        Currently, it only contains the final forces, which can be obtain
        by the id `final_forces`.

        """

        final_forces = self.final_forces
        forces = {'final': final_forces}

        return forces

    @property
    def maximum_force(self):
        """Fetch the maximum force of at the last ionic run."""

        forces = self.final_forces
        if forces is None:
            self._vasp_parser.exit_status = self._vasp_parser.exit_codes.ERROR_NOT_ABLE_TO_PARSE_QUANTITY
            return None
        norm = np.linalg.norm(forces, axis=1)
        return np.amax(np.abs(norm))

    @property
    def last_stress(self):
        """
        Fetch stess.

        After or at the last recorded ionic step from parsevasp.

        """

        stress = self._xml.get_stress('final')
        return stress

    @property
    def final_stress(self):
        """
        Fetch stress.

        After or at the last recorded ionic step from parsevasp.

        """

        return self.last_stress

    @property
    def stress(self):
        """
        Fetch stress.

        This container should contain all relevant stress.
        Currently, it only contains the final stress, which can be obtain
        by the id `final_stress`.

        """

        final_stress = self.final_stress
        stress = {'final': final_stress}
        return stress

    @property
    def maximum_stress(self):
        """Fetch the maximum stress of at the last ionic run."""

        stress = self.final_stress
        if stress is None:
            self._vasp_parser.exit_status = self._vasp_parser.exit_codes.ERROR_NOT_ABLE_TO_PARSE_QUANTITY
            return None
        norm = np.linalg.norm(stress, axis=1)
        return np.amax(np.abs(norm))

    @property
    def trajectory(self):
        """
        Fetch unitcells, positions, species, forces and stress.

        For all calculation steps from parsevasp.

        """

        unitcell = self._xml.get_unitcell('all')
        positions = self._xml.get_positions('all')
        species = self._xml.get_species()
        forces = self._xml.get_forces('all')
        stress = self._xml.get_stress('all')
        # make sure all are sorted, first to last calculation
        # (species is constant)
        unitcell = sorted(unitcell.items())
        positions = sorted(positions.items())
        forces = sorted(forces.items())
        stress = sorted(stress.items())
        # convert to numpy
        unitcell = np.asarray([item[1] for item in unitcell])
        positions = np.asarray([item[1] for item in positions])
        forces = np.asarray([item[1] for item in forces])
        stress = np.asarray([item[1] for item in stress])
        # Aiida wants the species as symbols, so invert
        elements = _invert_dict(parsevaspct.elements)
        symbols = np.asarray([elements[item].title() for item in species])

        if (unitcell is not None) and (positions is not None) and \
           (species is not None) and (forces is not None) and \
           (stress is not None):
            trajectory_data = {}

            keys = ('cells', 'positions', 'symbols', 'forces', 'stress', 'steps')
            stepids = np.arange(unitcell.shape[0])

            for key, data in zip(keys, (unitcell, positions, symbols, forces, stress, stepids)):
                trajectory_data[key] = data
            return trajectory_data

        self._vasp_parser.exit_status = self._vasp_parser.exit_codes.ERROR_NOT_ABLE_TO_PARSE_QUANTITY
        return None

    @property
    def total_energies(self):
        """Fetch the total energies after the last ionic run."""

        energies = self.energies
        if energies is None:
            self._vasp_parser.exit_status = self._vasp_parser.exit_codes.ERROR_NOT_ABLE_TO_PARSE_QUANTITY
            return None
        # fetch the type of energies that the user wants to extract
        settings = self._parsed_data.get('settings', DEFAULT_OPTIONS)
        energies_dict = {}
        for etype in settings.get('energy_type', DEFAULT_OPTIONS['energy_type']):
            energies_dict[etype] = energies[etype][-1]

        return energies_dict

    @property
    def energies(self):
        return self._energies(nosc=True)

    @property
    def energies_sc(self):
        """
        Fetch the total energies.

        Store in ArrayData for all self-consistent electronic steps.

        """

        # raise error due to lack of knowledge if
        # the Aiida data structure support for instance
        # lists of ndarrays.
        raise NotImplementedError
        #return self.energies(nosc = False)

    def _energies(self, nosc):
        """Fetch the total energies for all calculations (i.e. ionic steps)."""

        # fetch the type of energies that the user wants to extract
        settings = self._parsed_data.get('settings', DEFAULT_OPTIONS)

        enrgy = {}
        for etype in settings.get('energy_type', DEFAULT_OPTIONS['energy_type']):

            # this returns a list, not an ndarray due to
            # the posibility of returning the energies for all
            # self consistent steps, which contain a different
            # number of elements, not supported by Numpy's std.
            # arrays
            enrgies = self._xml.get_energies(status='all', etype=etype, nosc=nosc)
            if enrgies is None:
                self._vasp_parser.exit_status = self._vasp_parser.exit_codes.ERROR_NOT_ABLE_TO_PARSE_QUANTITY
                return None
            # should be a list, but convert to ndarray, here
            # staggered arrays are not a problem
            # two elements for a static run, both are similar,
            # only take the last
            if len(enrgies) == 2:
                enrgies = enrgies[-1:]
            enrgy[etype] = np.asarray(enrgies)

        return enrgy

    @property
    def projectors(self):
        """Fetch the projectors."""

        proj = self._xml.get_projectors()
        if proj is None:
            self._vasp_parser.exit_status = self._vasp_parser.exit_codes.ERROR_NOT_ABLE_TO_PARSE_QUANTITY
            return None
        projectors = {}
        prj = []
        try:
            prj.append(proj['total'])
        except KeyError:
            try:
                prj.append(proj['up'])
                prj.append(proj['down'])
            except KeyError:
                self._logger.error('Did not detect any projectors. Returning.')
        if len(prj) == 1:
            projectors['projectors'] = prj[0]
        else:
            projectors['projectors'] = np.asarray(prj)

        return projectors

    @property
    def dielectrics(self):
        """Fetch the dielectric function."""

        diel = self._xml.get_dielectrics()
        if diel is None:
            self._vasp_parser.exit_status = self._vasp_parser.exit_codes.ERROR_NOT_ABLE_TO_PARSE_QUANTITY
            return None
        dielectrics = {}
        energy = diel.get('energy')
        idiel = diel.get('imag')
        rdiel = diel.get('real')
        epsilon = diel.get('epsilon')
        epsilon_ion = diel.get('epsilon_ion')
        if energy is not None:
            dielectrics['ediel'] = energy
        if idiel is not None:
            dielectrics['rdiel'] = rdiel
        if rdiel is not None:
            dielectrics['idiel'] = idiel
        if epsilon is not None:
            dielectrics['epsilon'] = epsilon
        if epsilon_ion is not None:
            dielectrics['epsilon_ion'] = epsilon_ion

        return dielectrics

    @property
    def born_charges(self):
        """Fetch the Born effective charges."""

        brn = self._xml.get_born()
        if brn is None:
            self._vasp_parser.exit_status = self._vasp_parser.exit_codes.ERROR_NOT_ABLE_TO_PARSE_QUANTITY
            return None
        born = {'born_charges': brn}
        return born

    @property
    def hessian(self):
        """Fetch the Hessian matrix."""

        hessian = self._xml.get_hessian()
        if hessian is None:
            self._vasp_parser.exit_status = self._vasp_parser.exit_codes.ERROR_NOT_ABLE_TO_PARSE_QUANTITY
            return None
        hess = {'hessian': hessian}
        return hess

    @property
    def dynmat(self):
        """Fetch the dynamical eigenvectors and eigenvalues."""

        dynmat = self._xml.get_dynmat()
        if dynmat is None:
            self._vasp_parser.exit_status = self._vasp_parser.exit_codes.ERROR_NOT_ABLE_TO_PARSE_QUANTITY
            return None
        dyn = {}
        dyn['dynvec'] = dynmat['eigenvectors']
        dyn['dyneig'] = dynmat['eigenvalues']
        return dyn

    @property
    def dos(self):
        """Fetch the total density of states."""

        dos = self._xml.get_dos()
        if dos is None:
            self._vasp_parser.exit_status = self._vasp_parser.exit_codes.ERROR_NOT_ABLE_TO_PARSE_QUANTITY
            return None
        densta = {}
        # energy is always there, regardless of
        # total, spin or partial
        energy = dos['total']['energy']
        densta['energy'] = energy
        tdos = None
        pdos = None
        upspin = dos.get('up')
        downspin = dos.get('down')
        total = dos.get('total')
        if (upspin is not None) and (downspin is not None):
            tdos = np.stack((upspin['total'], downspin['total']))
            if (upspin['partial'] is not None) and \
               (downspin['partial'] is not None):
                pdos = np.stack((upspin['partial'], downspin['partial']))
        else:
            tdos = total['total']
            pdos = total['partial']
        densta['tdos'] = tdos
        if pdos is not None:
            densta['pdos'] = pdos

        return densta

    @property
    def fermi_level(self):
        """Fetch Fermi level."""

        return self._xml.get_fermi_level()


def _build_structure(lattice):
    """Builds a structure according to AiiDA spec."""
    structure_dict = {}
    structure_dict['unitcell'] = lattice['unitcell']
    structure_dict['sites'] = []

    # AiiDA wants the species as symbols, so invert
    elements = _invert_dict(parsevaspct.elements)
    for pos, specie in zip(lattice['positions'], lattice['species']):
        site = {}
        site['position'] = np.dot(pos, lattice['unitcell'])
        site['symbol'] = elements[specie].title()
        site['kind_name'] = elements[specie].title()
        structure_dict['sites'].append(site)

    return structure_dict


def _invert_dict(dct):
    return dct.__class__(map(reversed, dct.items()))
