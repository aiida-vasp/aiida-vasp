"""Tools for parsing vasprun.xml files."""
import numpy as np

from parsevasp.xml import Xml
from aiida_vasp.io.parser import BaseFileParser
from aiida_vasp.utils.aiida_utils import get_data_class

DEFAULT_OPTIONS = {
    'quantities_to_parse': ['occupations', 'vrp_pdos', 'vrp_tdos'],
    'output_params': ['energies', 'forces', 'efermi'],
}


class VasprunParser(BaseFileParser):
    """Wrapper for parsevasps Xml class parsing vasprun.xml files."""

    PARSABLE_ITEMS = {
        'occupations': {
            'inputs': [],
            'parsers': ['vasprun.xml'],
            'nodeName': 'intermediate_data',
            'prerequisites': []
        },
        'vrp_pdos': {
            'inputs': [],
            'parsers': ['vasprun.xml'],
            'nodeName': 'intermediate_data',
            'prerequisites': []
        },
        'vrp_tdos': {
            'inputs': [],
            'parsers': ['vasprun.xml'],
            'nodeName': 'intermediate_data',
            'prerequisites': []
        },
        'vrp_parameters': {
            'inputs': ['ocp_parameters'],
            'parsers': ['vasprun.xml'],
            'nodeName': 'intermediate_data',
            'prerequisites': []
        },
    }

    def __init__(self, *args, **kwargs):
        super(VasprunParser, self).__init__(*args, **kwargs)
        self.init_with_kwargs(**kwargs)

    def _init_with_path(self, filepath):
        """Init with a filepath."""
        self._parsed_data = {}
        self._parsable_items = VasprunParser.PARSABLE_ITEMS

        # Since vasprun.xml can be fairly large, we will parse it only
        # once and store the parsevasp Xml object instead of the filepath.
        try:
            self._data_obj = Xml(None, file_path=filepath)
        except SystemExit:
            self._data_obj = None

    @property
    def _parsed_object(self):
        """The vasprunParser will most likely not write a vasprun.xml file."""

        return None

    def _parse_file(self, inputs):

        # Since all qunatiies will be returned by properties, we can't pass
        # inputs as a parameter, so we store them in self._parsed_data
        for key, value in inputs.iteritems():
            self._parsed_data[key] = value

        settings = inputs.get('settings', DEFAULT_OPTIONS)
        if not settings:
            settings = DEFAULT_OPTIONS

        quantities_to_parse = settings.get('quantities_to_parse', DEFAULT_OPTIONS['quantities_to_parse'])

        result = {}
        for quantity in quantities_to_parse:
            if quantity in self._parsable_items:
                result[quantity] = getattr(self, quantity)

        return result

    @property
    def vrp_bands(self):
        """Return a BandsData node containing the bandstructure parsed from vasprun.xml."""

        eigenvalues = self.eigenvalues
        occupations = self.occupations

        if eigenvalues is None and occupations is None:
            return None

        bands = get_data_class('array.bands')()
        kpoints, weights = self.kpoints_and_weights

        bands.set_kpointsdata(kpoints, weights=weights, cartesian=False)
        bands.set_bands(eigenvalues, occupations=occupations)

        return bands

    @property
    def kpoints_and_weights(self):
        """Fetch the kpoints and weights."""
        return self._data_obj.get_kpoints(), self._data_obj.get_kpointsw()

    @property
    def eigenvalues(self):
        """Get eigenvalues from vasprun.xml."""
        eigenvalues = self._data_obj.get_eigenvalues()

        if eigenvalues is None:
            return None

        eigen = []
        eigen.append(eigenvalues.get("total"))

        if eigen[0] is None:
            eigen[0] = eigenvalues.get("up")
            eigen.append(eigenvalues.get("down"))

        if eigen[0] is None:
            return None

        return eigen

    @property
    def occupations(self):
        """Get occupations from vasprun.xml."""

        occupations = self._data_obj.get_occupancies()

        if occupations is None:
            return None

        occ = []
        occ.append(occupations.get("total"))

        if occ[0] is None:
            occ[0] = occupations.get("up")
            occ.append(occupations.get("down"))

        if occ[0] is None:
            return None

        return occ

    @property
    def efermi(self):
        """Get the Fermi energy from vasprun.xml."""
        return self._data_obj.get_fermi_level()

    @property
    def energies(self, nosc=True):
        """Fetch the total energies."""

        # energy without entropy
        etype = "energy_no_entropy"
        energies = self._data_obj.get_energies(status="all", etype=etype, nosc=nosc)

        if energies is None:
            return None

        # two elements for a static run, both are similar,
        # only take the last
        if len(energies) == 2:
            return np.asarray(energies[-1:])

        return np.asarray(energies)

    @property
    def stress(self):
        """Fetch last recorded stress."""
        stress = self._data_obj.get_stress("final")
        return stress

    @property
    def forces(self):
        """Fetch the last recorded forces"""
        forces = self._data_obj.get_forces("final")
        return forces

    @property
    def vrp_tdos(self):
        """Get the total dos from vasprun.xml."""

        dos = self._data_obj.get_dos()

        if dos is None:
            return np.array([])

        tdos = []
        tdos.append(dos.get("energy"))

        total = dos.get('total')
        dosup = dos.get('up')
        down = dos.get('down')

        if dosup is None and down is None:
            tdos.append(total['total'])
        else:
            tdos.append(dosup['total'])
            tdos.append(down['total'])

        return np.array(tdos)

    @property
    def vrp_pdos(self):
        """The partial DOS array"""
        dos = self._data_obj.get_dos()

        if dos is None:
            return np.array([])

        pdos = []
        pdos.append(dos.get("energy"))

        total = dos.get('total')
        dosup = dos.get('up')
        down = dos.get('down')

        if dosup is None and down is None:
            pdos.append(total['partial'])
        else:
            pdos.append(dosup['partial'])
            pdos.append(down['partial'])

        return np.array(pdos)

    @property
    def vrp_parameters(self):
        """Assemble the 'output_parameters' node."""
        parameters = {}
        parameters.update(self._parsed_data.get('ocp_parameters'))

        settings = self._parsed_data.get('settings', DEFAULT_OPTIONS)
        for quantity in settings.get('output_params', DEFAULT_OPTIONS['output_params']):
            parameters[quantity] = getattr(self, quantity)

        return parameters
