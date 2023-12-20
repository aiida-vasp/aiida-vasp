"""
The vasprun.xml parser interface.

---------------------------------
Contains the parsing interfaces to ``parsevasp`` used to parse ``vasprun.xml`` content.
"""
# pylint: disable=abstract-method, too-many-public-methods
import numpy as np
from parsevasp import constants as parsevaspct
from parsevasp.vasprun import Xml

from aiida_vasp.parsers.content_parsers.base import BaseFileParser
from aiida_vasp.utils.compare_bands import get_band_properties


class VasprunParser(BaseFileParser):
    """The parser interface that enables parsing of ``vasprun.xml`` content.

    The parser is triggered by using the keys listed in ``PARSABLE_QUANTITIES``.

    """

    OPEN_MODE = 'rb'
    DEFAULT_SETTINGS = {
        'quantities_to_parse': [
            'structure',
            'eigenvalues',
            'dos',
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
            'band_properties',
            'version',
        ],
        'energy_type': ['energy_extrapolated'],
        'electronic_step_energies':
        False
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
        'band_properties': {
            'inputs': [],
            'name': 'band_properties',
            'prerequisites': [],
        },
        'version': {
            'inputs': [],
            'name': 'version',
            'prerequisites': [],
        }
    }

    # Mapping of the energy names to those returned by parsevasp.vasprunl.Xml
    ENERGY_MAPPING = {
        'energy_extrapolated': 'energy_extrapolated_final',
        'energy_free': 'energy_free_final',
        'energy_no_entropy': 'energy_no_entropy_final',
        'energy_extrapolated_electronic': 'energy_extrapolated',
        'energy_free_electronic': 'energy_free',
        'energy_no_entropy_electronic': 'energy_no_entropy',
    }

    ENERGY_MAPPING_VASP5 = {
        'energy_extrapolated': 'energy_no_entropy_final',
        'energy_free': 'energy_free_final',
        # Not that energy_extrapolated_final parsed is the entropy term
        'energy_no_entropy': 'energy_extrapolated_final',
        'energy_extrapolated_electronic': 'energy_extrapolated',
        'energy_free_electronic': 'energy_free',
        'energy_no_entropy_electronic': 'energy_no_entropy',
    }

    def _init_from_handler(self, handler):
        """Initialize using a file like handler."""

        self.overflow = False
        try:
            self._content_parser = Xml(file_handler=handler, k_before_band=True, logger=self._logger)
        except SystemExit as exception:
            if exception.code == 509:
                # Xml might be fine but overflow is detected
                self.overflow = True
                self._logger.warning('Parsevasp exited abnormally due to overflow in XML file.')
            else:
                self._logger.warning('Parsevasp exited abnormally.')

    @property
    def version(self):
        """Fetch the VASP version from ``parsevasp`` and return it as a string object."""

        # fetch version
        version = self._content_parser.get_version()

        if version is None:
            return None

        return version

    @property
    def eigenvalues(self):
        """Fetch eigenvalues."""

        # Fetch eigenvalues
        eigenvalues = self._content_parser.get_eigenvalues()

        if eigenvalues is None:
            return None

        return eigenvalues

    @property
    def occupancies(self):
        """Fetch occupancies."""

        # Fetch occupancies
        occupancies = self._content_parser.get_occupancies()

        if occupancies is None:
            # occupancies not present, should not really happen?
            return None

        return occupancies

    @property
    def kpoints(self):
        """Fetch the kpoints an prepare for consumption by the NodeComposer."""

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

        After or at the last recorded ionic step.

        """

        last_lattice = self._content_parser.get_lattice('last')
        if last_lattice is None:
            return None
        return _build_structure(last_lattice)

    @property
    def final_structure(self):
        """
        Fetch the structure.

        After or at the last recorded ionic step. Should in
        principle be the same as the method above.

        """

        return self.last_structure

    @property
    def last_forces(self):
        """
        Fetch forces.

        After or at the last recorded ionic step.

        """

        force = self._content_parser.get_forces('last')
        return force

    @property
    def final_forces(self):
        """
        Fetch forces.

        After or at the last recorded ionic step.

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

        After or at the last recorded ionic step.

        """

        stress = self._content_parser.get_stress('last')
        return stress

    @property
    def final_stress(self):
        """
        Fetch stress.

        After or at the last recorded ionic step.

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

        For all calculation steps.

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
        for etype in self._settings.get('energy_type', self.DEFAULT_SETTINGS['energy_type']):
            energies_dict[etype] = energies[etype][-1]
            # Also return the raw electronic steps energy
            energies_dict[etype + '_electronic'] = energies[etype + '_electronic'][-1]

        return energies_dict

    @property
    def energies(self):
        """Fetch the total energies."""
        # Check if we want total energy entries for each electronic step.
        electronic_step_energies = self._settings.get(
            'electronic_step_energies', self.DEFAULT_SETTINGS['electronic_step_energies']
        )

        return self._energies(nosc=not electronic_step_energies)

    def _energies(self, nosc):
        """
        Fetch the total energies for all energy types, calculations (ionic steps) and electronic steps.

        The returned dict from the parser contains the total energy types as a key (plus the _final, which is
        the final total energy ejected by VASP after the closure of the electronic steps). The energies can then
        be found in the flattened ndarray where the key `electronic_steps` indicate how many electronic steps
        there is per ionic step. Using the combination, one can rebuild the electronic step energy per ionic step etc.

        Because the VASPrun parser returns both the electronic step energies (at the end of each cycles) and the ionic step
        energies (_final), we apply a mapping to recovery the naming such that the ionic step energies do not have the suffix,
        but the electronic step energies do.
        """

        etype = self._settings.get('energy_type', self.DEFAULT_SETTINGS['energy_type'])
        # Create a copy
        etype = list(etype)
        etype_orig = list(etype)

        # Apply mapping and request the correct energies from the parsing results
        # VASP 5 has a bug where the energy_no_entropy is not included in the XML output - we have to calculate it here
        if self.version.startswith('5'):
            # For energy_no_entropy needs to be calculated here
            if 'energy_no_entropy' in etype_orig:
                etype.append('energy_free')
                etype.append('energy_extrapolated')

            # energy extrapolated is stored as energy_no_entropy for the ionic steps
            if 'energy_extrapolated' in etype_orig:
                etype.append('energy_no_entropy')

            # Remove duplicates
            etype = list(set(etype))
            energies = self._content_parser.get_energies(status='all', etype=etype, nosc=nosc)
            # Here we must calculate the true `energy_no_entropy`
            if 'energy_no_entropy' in etype_orig:
                # The energy_extrapolated_final is the entropy term itself in VASP 5
                # Store the calculated energy_no_entropy under 'energy_extrapolated_final',
                # which is then recovered as `energy_no_entropy` later
                energies['energy_extrapolated_final'
                         ] = energies['energy_free_final'] - energies['energy_extrapolated_final']
        else:
            energies = self._content_parser.get_energies(status='all', etype=etype, nosc=nosc)

        if energies is None:
            return None

        # Apply mapping - those with `_final` has the suffix removed and those without has `_electronic` added
        mapped_energies = {}
        mapping = self.ENERGY_MAPPING_VASP5 if self.version.startswith('5') else self.ENERGY_MAPPING
        # Reverse the mapping - now key is the name of the original energies output
        revmapping = {value: key for key, value in mapping.items()}
        for key, value in energies.items():
            # Apply mapping if needed
            if key in revmapping:
                if revmapping[key].replace('_electronic', '') in etype_orig:
                    mapped_energies[revmapping[key]] = value
            else:
                mapped_energies[key] = value

        return mapped_energies

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
        info['finished'] = not self._content_parser.truncated
        # Only set to true for untruncated run to avoid false positives
        if energies is None:
            info['electronic_converged'] = False
        elif energies.get('electronic_steps')[0] < parameters['nelm'] and not self._content_parser.truncated:
            info['electronic_converged'] = True
        else:
            info['electronic_converged'] = False

        # Then check the ionic convergence by comparing executed steps to the
        # maximum allowed number of steps (NSW).
        energies = self._content_parser.get_energies('all', nosc=True)
        if energies is None:
            info['ionic_converged'] = False
        else:
            if len(energies.get('electronic_steps')) < parameters['nsw'] and not self._content_parser.truncated:
                info['ionic_converged'] = True
            else:
                info['ionic_converged'] = False
        # Override if nsw is 0 - no ionic steps are performed
        if parameters['nsw'] < 1:
            info['ionic_converged'] = None

        return info

    @property
    def band_properties(self):
        """Fetch key properties of the electronic structure."""

        eigenvalues = self.eigenvalues
        occupancies = self.occupancies
        if eigenvalues is None:
            return None

        # Convert dict to index in numpy array
        if 'total' in eigenvalues:
            eig = np.array(eigenvalues['total'])
            occ = np.array(occupancies['total'])
        else:
            eig = np.array([eigenvalues['up'], eigenvalues['down']])
            occ = np.array([occupancies['up'], occupancies['down']])

        return get_band_properties(eig, occ)


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
