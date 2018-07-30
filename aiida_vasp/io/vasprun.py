"""Tools for parsing vasprun.xml files."""
import operator
import numpy as np

from parsevasp.vasprun import Xml
from parsevasp import constants as parsevaspct
from aiida_vasp.io.parser import BaseFileParser
from aiida_vasp.utils.aiida_utils import get_data_class, get_data_node

DEFAULT_OPTIONS = {
    'quantities_to_parse': [
        'parameters', 'structure', 'bands', 'dos', 'kpoints', 'occupations', 'trajectory', 'energies', 'projectors', 'dielectrics',
        'born_charges', 'hessian', 'dynmat'
    ],
    'energy_type': ['energy_no_entropy'],
    'output_params': []
}


class ExtendedXml(Xml):
    """
    Extension of parsevasp's Xml class.

    To keep parser interfaces similar.

    """

    @property
    def path(self):
        return self._file_path

    def write(self, dst):
        """Copy vasprun.xml to destination."""
        import shutil
        shutil.copyfile(self._file_path, dst)


class VasprunParser(BaseFileParser):
    """Interface to parsevasp's xml parser."""

    PARSABLE_ITEMS = {
        'structure': {
            'inputs': [],
            'parsers': ['vasprun.xml'],
            'nodeName': 'structure',
            'prerequisites': [],
            'alternatives': ['poscar-structure']
        },
        'bands': {
            'inputs': [],
            'parsers': ['vasprun.xml'],
            'nodeName': 'bands',
            'prerequisites': [],
            'alternatives': ['eigenval-bands']
        },
        'dos': {
            'inputs': [],
            'parsers': ['vasprun.xml'],
            'nodeName': 'dos',
            'prerequisites': [],
            'alternatives': ['doscar-dos']
        },
        'kpoints': {
            'inputs': [],
            'parsers': ['vasprun.xml'],
            'nodeName': 'kpoints',
            'prerequisites': [],
            'alternatives': ['kpoints-kpoints']
        },
        'occupations': {
            'inputs': [],
            'parsers': ['vasprun.xml'],
            'nodeName': 'occupations',
            'prerequisites': [],
            #'alternatives': ['eigenval-occupations']
        },
        'trajectory': {
            'inputs': [],
            'parsers': ['vasprun.xml'],
            'nodeName': 'trajectory',
            'prerequisites': [],
            #'alternatives': ['xdatcar-trajectory']
        },
        'energies': {
            'inputs': [],
            'parsers': ['vasprun.xml'],
            'nodeName': 'energies',
            'prerequisites': [],
            'alternatives': ['outcar-energies']
        },
        'projectors': {
            'inputs': [],
            'parsers': ['vasprun.xml'],
            'nodeName': 'projectors',
            'prerequisites': [],
            #'alternatives': ['procar-projectors']
        },
        'dielectrics': {
            'inputs': [],
            'parsers': ['vasprun.xml'],
            'nodeName': 'dielectrics',
            'prerequisites': [],
            #'alternatives': ['outcar-dielectrics']
        },
        'final_stress': {
            'inputs': [],
            'parsers': ['vasprun.xml'],
            'nodeName': '',
            'prerequisites': [],
            #'alternatives': ['outcar-final_stress']
        },
        'final_forces': {
            'inputs': [],
            'parsers': ['vasprun.xml'],
            'nodeName': '',
            'prerequisites': [],
            #'alternatives': ['outcar-final_stress']
        },
        'final_structure': {
            'inputs': [],
            'parsers': ['vasprun.xml'],
            'nodeName': '',
            'prerequisites': [],
            #'alternatives': ['outcar-final_structure']
        },
        'born_charges': {
            'inputs': [],
            'parsers': ['vasprun.xml'],
            'nodeName': 'born_charges',
            'prerequisites': [],
            #'alternatives': ['outcar-born_charges']
        },
        'hessian': {
            'inputs': [],
            'parsers': ['vasprun.xml'],
            'nodeName': 'hessian',
            'prerequisites': [],
            #'alternatives': ['outcar-hessian']
        },
        'dynmat': {
            'inputs': [],
            'parsers': ['vasprun.xml'],
            'nodeName': 'dynmat',
            'prerequisites': [],
            #'alternatives': ['outcar-dynmat']
        },
        'parameters': {
            'inputs': [],
            'parsers': ['vasprun.xml'],
            'nodeName': 'parameters',
            'prerequisites': [],
            'alternatives': ['outcar-parameters']
        }
    }

    def __init__(self, *args, **kwargs):
        super(VasprunParser, self).__init__(*args, **kwargs)
        self.init_with_kwargs(**kwargs)
        self._settings = {}

    def _init_with_file_path(self, path):
        """Init with a filepath."""
        self._parsed_data = {}
        self._parsable_items = self.__class__.PARSABLE_ITEMS

        # Since vasprun.xml can be fairly large, we will parse it only
        # once and store the parsevasp Xml object instead of the filepath.
        # Also, make sure the k-point index runs before the band index as
        # per BandsData spec.
        try:
            self._data_obj = ExtendedXml(file_path=path, k_before_band=True)
        except SystemExit:
            self._logger.warning("Parsevasp exited abruptly. Returning None.")
            self._data_obj = None

    def _init_with_data(self, data):
        """Init with singleFileData."""
        self._parsable_items = self.__class__.PARSABLE_ITEMS
        self._init_with_file_path(data.get_file_abs_path())

    def _parse_file(self, inputs):

        # Since all quantities will be returned by properties, we can't pass
        # inputs as a parameter, so we store them in self._parsed_data
        for key, value in inputs.items():
            self._parsed_data[key] = value

        settings = inputs.get('settings', DEFAULT_OPTIONS)
        if not settings:
            settings = DEFAULT_OPTIONS
        self._settings = settings

        quantities_to_parse = settings.get('quantities_to_parse', DEFAULT_OPTIONS['quantities_to_parse'])
        result = {}
        for quantity in quantities_to_parse:
            if quantity in self._parsable_items:
                result[quantity] = getattr(self, quantity)

        return result

    @property
    def bands(self):
        """
        Return a BandsData node.

        Contains the bandstructure parsed from vasprun.xml.

        """

        # fetch eigenvalues and occupancies
        eigenvalues = self.eigenvalues
        occupations = self.occupations_bands

        if eigenvalues is None:
            # did not find any eigenvalues
            return None

        # generate Aiida specific BandsData for storage
        band_data = get_data_class('array.bands')()

        # put everything into BandData and KpointsData
        band_data.set_kpointsdata(self.kpoints)
        band_data.set_bands(eigenvalues, occupations=occupations)

        return band_data

    @property
    def eigenvalues(self):
        """Fetch eigenvalues from parsevasp."""

        # fetch eigenvalues
        eigenvalues = self._data_obj.get_eigenvalues()

        if eigenvalues is None:
            # eigenvalues not present
            return None

        eigen = []
        eigen.append(eigenvalues.get("total"))

        if eigen[0] is None:
            # spin decomposed?
            eigen[0] = eigenvalues.get("up")
            eigen.append(eigenvalues.get("down"))

        if eigen[0] is None:
            # safety, should not really happen?
            return None

        return eigen

    @property
    def occupations_bands(self):
        """Fetch occupations from parsevasp."""

        # fetch occupations
        occupations = self._data_obj.get_occupancies()

        if occupations is None:
            # occupations not present, should not really happen?
            return None

        occ = []
        occ.append(occupations.get("total"))

        if occ[0] is None:
            # spin decomposed
            occ[0] = occupations.get("up")
            occ.append(occupations.get("down"))

        if occ[0] is None:
            # should not really happen
            return None

        return occ

    @property
    def occupations(self):
        """Fetch occupations from parsevasp."""

        # fetch occupations
        occupations = self._data_obj.get_occupancies()

        if occupations is None:
            # occupations not present, should not really happen?
            return None

        array_data = get_data_class('array')()

        total = occupations.get('total')
        upspin = occupations.get('up')
        downspin = occupations.get('down')
        if total is not None:
            # we have total
            array_data.set_array('total', total)
        elif upspin is not None:
            # we have spin decomposed
            array_data.set_array('up', upspin)
            if downspin is None:
                self._logger.error("Serious error, detected spin up, but no spin down " "channel. This should not happen. Continuing.")
            array_data.set_array('down', downspin)
        else:
            # safety, should not really happen?
            return None

        return array_data

    @property
    def parameters(self):
        """Assemble the 'output_params' node."""

        parameters = {}
        outcar_parameters = self._parsed_data.get('ocp_parameters')
        if outcar_parameters is not None:
            parameters.update(outcar_parameters)
        for quantity in self._settings.get('output_params', DEFAULT_OPTIONS['output_params']):
            parameters[quantity] = getattr(self, quantity)

        output_parameters = get_data_node('parameter', dict=parameters)

        return output_parameters

    @property
    def kpoints(self):
        """Fetch the kpoints from parsevasp an store in KpointsData."""

        kpts = self._data_obj.get_kpoints()
        kptsw = self._data_obj.get_kpointsw()
        kpoints_data = None
        if (kpts is not None) and (kptsw is not None):
            # create a KpointsData object and store k-points
            kpoints_data = get_data_class('array.kpoints')()
            kpoints_data.set_kpoints(kpts, weights=kptsw)

        return kpoints_data

    @property
    def structure(self):
        """
        Fetch a given structure and store as StructureData.

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

        After or at the last recorded ionic step from parsevasp
        and store as StructureData.

        """

        last_lattice = self._data_obj.get_lattice("final")
        if last_lattice is not None:
            return _build_structure(last_lattice)
        return None

    @property
    def final_structure(self):
        """
        Fetch the structure.

        After or at the last recorded ionic step from parsevasp
        and store as StructureData. Should in
        principle be the same as the method above.

        """

        return self.last_structure

    @property
    def last_forces(self):
        """
        Fetch forces.

        After or at the last recorded ionic step from parsevasp
        and store as ArrayData.

        """

        frs = self._data_obj.get_forces("final")
        if frs is None:
            return None

        return frs

    @property
    def final_forces(self):
        """
        Fetch forces.

        After or at the last recorded
        ionic step from parsevasp and store as ArrayData.

        """

        return self.last_forces

    @property
    def last_stress(self):
        """
        Fetch stess.

        After or at the last recorded ionic step from parsevasp
        and store as ArrayData.

        """

        strs = self._data_obj.get_stress("final")
        if strs is None:
            return None
        stress = get_data_class('array')()
        stress.set_array('forces', strs)
        return stress

    @property
    def final_stress(self):
        """
        Fetch stress.

        After or at the last recorded ionic step from parsevasp
        and store as ArrayData.

        """

        return self.last_stress

    @property
    def trajectory(self):
        """
        Fetch unitcells, positions, species, forces and stress.

        For all calculation steps from parsevasp and store as TrajectoryData.

        """

        unitcell = self._data_obj.get_unitcell("all")
        positions = self._data_obj.get_positions("all")
        species = self._data_obj.get_species()
        forces = self._data_obj.get_forces("all")
        stress = self._data_obj.get_stress("all")
        # make sure all are sorted, first to last calculation
        # (species is constant)
        unitcell = sorted(unitcell.items())
        positions = sorted(positions.items())
        forces = sorted(forces.items())
        stress = sorted(stress.items())
        # convert to numpy
        unitcell = np.asarray(map(operator.itemgetter(1), unitcell))
        positions = np.asarray(map(operator.itemgetter(1), positions))
        forces = np.asarray(map(operator.itemgetter(1), forces))
        stress = np.asarray(map(operator.itemgetter(1), stress))
        # Aiida wants the species as symbols, so invert
        elements = _invert_dict(parsevaspct.elements)
        symbols = np.asarray([elements[item].title() for item in species])

        if (unitcell is not None) and (positions is not None) and \
           (species is not None) and (forces is not None) and \
           (stress is not None):
            array_node = get_data_class('array')()
            trajectory_node = get_data_class('array.trajectory')()

            keys = ('cells', 'positions', 'symbols', 'forces', 'stress')

            trajectory_node.set_trajectory(stepids=np.arange(unitcell.shape[0]), cells=unitcell, symbols=symbols, positions=positions)

            for key, data in zip(keys, (unitcell, positions, symbols, forces, stress)):
                array_node.set_array(key, data)
                trajectory_node.set_array(key, data)
            return trajectory_node, array_node

        return None

    # @property
    # def trajectory_full(self):
    #     """Fetch unitcells, positions, species, forces and stress
    #     for all calculation steps from parsevasp and store as ArrayData.

    #     """

    #     unitcell = self._data_obj.get_unitcell("all")
    #     positions = self._data_obj.get_positions("all")
    #     species = self._data_obj.get_species()
    #     forces = self._data_obj.get_forces("all")
    #     stress = self._data_obj.get_stress("all")
    #     make sure all are sorted, first to last calculation
    #     (species is constant)
    #     unitcell = sorted(unitcell.items())
    #     positions = sorted(positions.items())
    #     forces = sorted(forces.items())
    #     stress = sorted(stress.items())
    #     convert to numpy
    #     unitcell = np.asarray(map(operator.itemgetter(1),unitcell))
    #     positions = np.asarray(map(operator.itemgetter(1),positions))
    #     forces = np.asarray(map(operator.itemgetter(1),forces))
    #     stress = np.asarray(map(operator.itemgetter(1),stress))
    #     Aiida wants the species as symbols, so invert
    #     elements = self._invert_dict(parsevaspct.elements)
    #     symbols = np.asarray([elements[item].title() for item in species])

    #     if (unitcell is not None) and (positions is not None) and \
    #        (species is not None) and (forces is not None) and \
    #        (stress is not None):
    #         array_node = get_data_class('array')()

    #         keys = ('cells', 'positions', 'symbols', 'forces', 'stress')

    #         for key, data in zip(keys, (unitcell,
    #                                     positions,
    #                                     symbols,
    #                                     forces,
    #                                     stress)):
    #             array_node.set_array(key, data)
    #         return array_node
    #     else:
    #         return None

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

    @property
    def energies(self, nosc=True):
        """
        Fetch the total energies.

        Store as ArrayData for all calculations (i.e. ionic steps).

        """

        # create a ArrayData object
        enrgy = get_data_class('array')()

        # fetch the type of energies that the user wants to extract
        settings = self._parsed_data.get('settings', DEFAULT_OPTIONS)
        for etype in settings.get('energy_type', DEFAULT_OPTIONS['energy_type']):

            # this returns a list, not an ndarray due to
            # the posibility of returning the energies for all
            # self consistent steps, which contain a different
            # number of elements, not supported by Numpy's std.
            # arrays
            enrgies = self._data_obj.get_energies(status="all", etype=etype, nosc=nosc)
            if enrgies is None:
                return None
            enrgy = get_data_class('array')()
            # should be a list, but convert to ndarray, here
            # staggered arrays are not a problem
            # two elements for a static run, both are similar,
            # only take the last
            if len(enrgies) == 2:
                enrgies = enrgies[-1:]
            enrgy.set_array(etype, np.asarray(enrgies))

        return enrgy

    @property
    def projectors(self):
        """
        Fetch the projectors.

        Store as ArrayData.

        """

        proj = self._data_obj.get_projectors()
        if proj is None:
            return None
        projectors = get_data_class('array')()
        prj = []
        try:
            prj.append(proj["total"])
        except KeyError:
            try:
                prj.append(proj["up"])
                prj.append(proj["down"])
            except KeyError:
                self._logger.error("Did not detect any projectors. " "Returning.")
        if len(prj) == 1:
            projectors.set_array('projectors', prj[0])
        else:
            projectors.set_array('projectors', np.asarray(prj))
        return projectors

    @property
    def dielectrics(self):
        """
        Fetch the dielectric function.

        Store as ArrayData.

        """

        diel = self._data_obj.get_dielectrics()
        if diel is None:
            return None
        dielectrics = get_data_class('array')()
        energy = diel["energy"]
        dielectrics.set_array('ediel', energy)
        dielectrics.set_array('rdiel', diel["real"])
        dielectrics.set_array('idiel', diel["imag"])
        return dielectrics

    @property
    def born_charges(self):
        """
        Fetch the Born effective charges.

        Store as ArrayData.

        """

        brn = self._data_obj.get_born()
        if brn is None:
            return None
        born = get_data_class('array')()
        born.set_array('born_charges', brn)
        return born

    @property
    def hessian(self):
        """
        Fetch the Hessian matrix.

        Store as ArrayData.

        """

        hessian = self._data_obj.get_hessian()
        if hessian is None:
            return None
        hess = get_data_class('array')()
        hess.set_array('hessian', hessian)
        return hess

    @property
    def dynmat(self):
        """
        Fetch the dynamical eigenvectors and eigenvalues.

        Store as ArrayData.

        """

        dynmat = self._data_obj.get_dynmat()
        if dynmat is None:
            return None
        dyn = get_data_class('array')()
        dyn.set_array('dynvec', dynmat["eigenvectors"])
        dyn.set_array('dyneig', dynmat["eigenvalues"])
        return dyn

    @property
    def dos(self):
        """
        Fetch the total density of states.

        Store as ArrayData.

        """

        dos = self._data_obj.get_dos()
        if dos is None:
            return None
        densta = get_data_class('array')()
        # energy is always there, regardless of
        # total, spin or partial
        energy = dos["total"]["energy"]
        densta.set_array('energy', energy)
        tdos = None
        pdos = None
        upspin = dos.get("up")
        downspin = dos.get("down")
        total = dos.get("total")
        if (upspin is not None) and (downspin is not None):
            tdos = np.stack((upspin["total"], downspin["total"]))
            if (upspin["partial"] is not None) and \
               (downspin["partial"] is not None):
                pdos = np.stack((upspin["partial"], downspin["partial"]))
        else:
            tdos = total["total"]
            pdos = total["partial"]
        densta.set_array('tdos', tdos)
        if pdos is not None:
            densta.set_array('pdos', pdos)

        return densta

    @property
    def fermi_level(self):
        """Fetch Fermi level."""

        return self._data_obj.get_fermi_level()


def _build_structure(lattice):
    """Builds a structure according to Aiida spec."""

    structure_cls = get_data_class('structure')
    unitcell = lattice["unitcell"]
    structure = structure_cls(cell=unitcell)
    # Aiida wants the species as symbols, so invert
    elements = _invert_dict(parsevaspct.elements)
    for pos, specie in zip(lattice["positions"], lattice["species"]):
        structure.append_atom(position=np.dot(pos, unitcell), symbols=elements[specie].title())
    return structure


def _invert_dict(dct):
    return dct.__class__(map(reversed, dct.items()))
