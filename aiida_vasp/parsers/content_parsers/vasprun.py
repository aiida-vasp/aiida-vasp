"""
vasprun.xml parser.

---------------
Contains the parsing interfaces to parsevasp used to parse vasprun.xml.
"""
# pylint: disable=too-many-public-methods, protected-access
import sys
import numpy as np

from parsevasp.vasprun import Xml
from parsevasp.kpoints import Kpoint
from parsevasp import constants as parsevaspct
from aiida_vasp.parsers.content_parsers.parser import BaseFileParser


class VasprunParser(BaseFileParser):
    """The parser interface that enables parsing of vasprun.xml.

    The parser is triggered by using the keys listed in PARSABLE_QUANTITIES.

    """

    DEFAULT_OPTIONS = {
        'quantities_to_parse': [
            'structure',
            'eigenvalues',
            'dos',
            'bands',
            'kpoints',
            'occupancies',
            'trajectory',
            'energies',
            'projectors',
            'dielectrics',
            'born_charges',
            'hessian',
            'dynmat',
            'forces',
            'stress',
            'total_energies',
            'maximum_force',
            'maximum_stress',
            'version',
        ],
        'energy_type': ['energy_extrapolated'],
        'electronic_step_energies': False
    }

    PARSABLE_QUANTITIES = {
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
            'name': 'dielectrics',
            'prerequisites': [],
        },
        'stress': {
            'inputs': [],
            'name': 'stress',
            'prerequisites': [],
        },
        'forces': {
            'inputs': [],
            'name': 'forces',
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
        'version': {
            'inputs': [],
            'name': 'version',
            'prerequisites': [],
        }
    }

    def _init_from_handler(self, handler):
        """Initialize using a file like handler."""

        try:
            self._content_parser = Xml(file_handler=handler, k_before_band=True, logger=self._logger)
        except SystemExit:
            self._logger.warning('Parsevasp exited abnormally.')

    @property
    def version(self):
        """Fetch the VASP version from parsevasp and return it as a string object."""

        # fetch version
        version = self._content_parser.get_version()

        if version is None:
            return None

        return version

    @property
    def eigenvalues(self):
        """Fetch eigenvalues from parsevasp."""

        # Fetch eigenvalues
        eigenvalues = self._content_parser.get_eigenvalues()

        if eigenvalues is None:
            return None

        return eigenvalues

    @property
    def occupancies(self):
        """Fetch occupancies from parsevasp."""

        # Fetch occupancies
        occupancies = self._content_parser.get_occupancies()

        if occupancies is None:
            # occupancies not present, should not really happen?
            return None

        return occupancies

    @property
    def kpoints(self):
        """Fetch the kpoints from parsevasp an prepare for consumption by the NodeComposer."""

        kpts = self._content_parser.get_kpoints()
        kptsw = self._content_parser.get_kpointsw()
        # k-points in XML is always in reciprocal if spacing methods have been used
        # but what about explicit/regular
        cartesian = False
        kpoints_data = None
        if (kpts is not None) and (kptsw is not None):
            # Create a dictionary and store k-points that can be consumed by the NodeComposer
            kpoints_data = {}
            kpoints_data['mode'] = 'explicit'
            kpoints_data['cartesian'] = cartesian
            kpoints_data['points'] = kpts
            kpoints_data['weights'] = kptsw

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

        last_lattice = self._content_parser.get_lattice('last')
        if last_lattice is None:
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

        force = self._content_parser.get_forces('last')
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
            return None
        norm = np.linalg.norm(forces, axis=1)
        return np.amax(np.abs(norm))

    @property
    def last_stress(self):
        """
        Fetch stess.

        After or at the last recorded ionic step from parsevasp.

        """

        stress = self._content_parser.get_stress('last')
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
            return None
        norm = np.linalg.norm(stress, axis=1)
        return np.amax(np.abs(norm))

    @property
    def trajectory(self):
        """
        Fetch unitcells, positions, species, forces and stress.

        For all calculation steps from parsevasp.

        """

        unitcell = self._content_parser.get_unitcell('all')
        positions = self._content_parser.get_positions('all')
        species = self._content_parser.get_species()
        forces = self._content_parser.get_forces('all')
        stress = self._content_parser.get_stress('all')
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
        symbols = np.asarray([elements[item].title() for item in species.tolist()])

        if (unitcell is not None) and (positions is not None) and \
           (species is not None) and (forces is not None) and \
           (stress is not None):
            trajectory_data = {}

            keys = ('cells', 'positions', 'symbols', 'forces', 'stress', 'steps')
            stepids = np.arange(unitcell.shape[0])

            for key, data in zip(keys, (unitcell, positions, symbols, forces, stress, stepids)):
                trajectory_data[key] = data
            return trajectory_data

        return None

    @property
    def total_energies(self):
        """Fetch the total energies after the last ionic run."""
        energies = self.energies
        if energies is None:
            return None
        energies_dict = {}
        for etype in self._settings.get('energy_type', self.DEFAULT_OPTIONS['energy_type']):
            energies_dict[etype] = energies[etype][-1]

        return energies_dict

    @property
    def energies(self):
        """Fetch the total energies."""
        # Check if we want total energy entries for each electronic step.
        electronic_step_energies = self._settings.get('electronic_step_energies', self.DEFAULT_OPTIONS['electronic_step_energies'])

        return self._energies(nosc=not electronic_step_energies)

    def _energies(self, nosc):
        """
        Fetch the total energies for all energy types, calculations (ionic steps) and electronic steps.

        The returned dict from the parser contains the total energy types as a key (plus the _final, which is
        the final total energy ejected by VASP after the closure of the electronic steps). The energies can then
        be found in the flattened ndarray where the key `electronic_steps` indicate how many electronic steps
        there is per ionic step. Using the combination, one can rebuild the electronic step energy per ionic step etc.

        """
        etype = self._settings.get('energy_type', self.DEFAULT_OPTIONS['energy_type'])
        energies = self._content_parser.get_energies(status='all', etype=etype, nosc=nosc)
        if energies is None:
            return None

        return energies

    @property
    def projectors(self):
        """Fetch the projectors."""

        proj = self._content_parser.get_projectors()
        if proj is None:
            return None
        projectors = {}
        prj = []
        try:
            prj.append(proj['total'])  # pylint: disable=unsubscriptable-object
        except KeyError:
            try:
                prj.append(proj['up'])  # pylint: disable=unsubscriptable-object
                prj.append(proj['down'])  # pylint: disable=unsubscriptable-object
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

        diel = self._content_parser.get_dielectrics()
        if diel is None:
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

        brn = self._content_parser.get_born()
        if brn is None:
            return None
        born = {'born_charges': brn}
        return born

    @property
    def hessian(self):
        """Fetch the Hessian matrix."""

        hessian = self._content_parser.get_hessian()
        if hessian is None:
            return None
        hess = {'hessian': hessian}
        return hess

    @property
    def dynmat(self):
        """Fetch the dynamical eigenvectors and eigenvalues."""

        dynmat = self._content_parser.get_dynmat()
        if dynmat is None:
            return None
        dyn = {}
        dyn['dynvec'] = dynmat['eigenvectors']  # pylint: disable=unsubscriptable-object
        dyn['dyneig'] = dynmat['eigenvalues']  # pylint: disable=unsubscriptable-object
        return dyn

    @property
    def dos(self):
        """Fetch the total density of states."""

        dos = self._content_parser.get_dos()
        if dos is None:
            return None
        densta = {}
        # energy is always there, regardless of
        # total, spin or partial
        energy = dos['total']['energy']  # pylint: disable=unsubscriptable-object
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

        return self._content_parser.get_fermi_level()

    @property
    def run_status(self):
        """Fetch run_status information"""
        info = {}
        # First check electronic convergence by comparing executed steps to the
        # maximum allowed number of steps (NELM).
        energies = self._content_parser.get_energies('last', nosc=False)
        parameters = self._content_parser.get_parameters()
        info['finished'] = not self._content_parser_truncated
        # Only set to true for untruncated run to avoid false positives
        if energies is None:
            info['electronic_converged'] = False
        elif energies.get('electronic_steps')[0] < parameters['nelm'] and not self._content_parser_truncated:
            info['electronic_converged'] = True
        else:
            info['electronic_converged'] = False

        # Then check the ionic convergence by comparing executed steps to the
        # maximum allowed number of steps (NSW).
        energies = self._content_parser.get_energies('all', nosc=True)
        if energies is None:
            info['ionic_converged'] = False
        else:
            if len(energies.get('electronic_steps')) < parameters['nsw'] and not self._content_parser_truncated:
                info['ionic_converged'] = True
            else:
                info['ionic_converged'] = False
        # Override if nsw is 0 - no ionic steps are performed
        if parameters['nsw'] < 1:
            info['ionic_converged'] = None

        return info


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
