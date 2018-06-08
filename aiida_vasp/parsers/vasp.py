#encoding: utf-8
"""AiiDA Parser for a aiida_vasp.VaspCalculation"""

from aiida_vasp.io.doscar import DosParser
from aiida_vasp.io.eigenval import EigParser
from aiida_vasp.io.kpoints import KpParser
from aiida_vasp.io.outcar import OutcarParser
from aiida_vasp.io.vasprun import VasprunParser
from aiida_vasp.io.chgcar import ChgcarParser
from aiida_vasp.io.wavecar import WavecarParser
from aiida_vasp.io.poscar import PoscarParser
from aiida_vasp.parsers.base import BaseParser
from aiida_vasp.utils.delegates import delegate
from aiida.orm.data.parameter import ParameterData
from aiida.orm.data.structure import StructureData
from aiida.orm.data.array.kpoints import KpointsData
from aiida.orm.data.array.bands import BandsData
from aiida.orm.data.array.trajectory import TrajectoryData
from aiida.orm.data.array import ArrayData

LINKNAME_DICT = {
    'parameters': 'output_parameters',
    'kpoints': 'output_kpoints',
    'structure': 'output_structure',
    'trajectory': 'output_trajectory',
    'trajectory_full': 'output_trajectory_arr',
    'bands': 'output_bands',
    'dos': 'output_dos',
    'energies': 'output_energies',
    'projectors': 'output_projectors',
    'born_charges': 'output_born_charges',
    'dielectrics': 'output_dielectrics',
    'hessian': 'output_hessian',
    'dynmat': 'output_dynmat',
    'chgcar': 'chgcar',
    'wavecar': 'wavecar',
}

ALLOWED_TYPES = {
    'parameters': ParameterData,
    'kpoints': KpointsData,
    'structure': StructureData,
    'trajectory': TrajectoryData,
    'trajectory_full': ArrayData,
    'bands': BandsData,
    'dos': ArrayData,
    'energies': ArrayData,
    'projectors': ArrayData,
    'born_charges': ArrayData,
    'dielectrics': ArrayData,
    'hessian': ArrayData,
    'dynmat': ArrayData
}

DEFAULT_OPTIONS = {
    'add_trajectory': False,
    'add_bands': False,
    'add_chgcar': False,
    'add_dos': False,
    'add_kpoints': True,
    'add_energies': True,
    'add_parameters': True,
    'add_structure': True,
    'add_projectors': False,
    'add_born_charges': False,
    'add_dielectrics': False,
    'add_hessian': False,
    'add_dynmat': False,
    'add_wavecar': False,
    'should_parse_DOSCAR': False,
    'should_parse_EIGENVAL': False,
    'should_parse_IBZKPT': False,
    'should_parse_OUTCAR': False,
    'should_parse_CONTCAR': False,
    'should_parse_vasprun.xml': True,
    'should_parse_WAVECAR': False,
    'should_parse_CHGCAR': False
}

PARSABLE_FILES = {
    'DOSCAR': {
        'parser_class': DosParser,
        'is_critical': False,
        'status': 'Unknown'
    },
    'EIGENVAL': {
        'parser_class': EigParser,
        'is_critical': False,
        'status': 'Unknown'
    },
    'IBZKPT': {
        'parser_class': KpParser,
        'is_critical': False,
        'status': 'Unknown'
    },
    'OUTCAR': {
        'parser_class': OutcarParser,
        'is_critical': False,
        'status': 'Unknown'
    },
    'vasprun.xml': {
        'parser_class': VasprunParser,
        'is_critical': True,
        'status': 'Unknown'
    },
    'CHGCAR': {
        'parser_class': ChgcarParser,
        'is_critical': False,
        'status': 'Unknown'
    },
    'WAVECAR': {
        'parser_class': WavecarParser,
        'is_critical': False,
        'status': 'Unknown'
    },
    'CONTCAR': {
        'parser_class': PoscarParser,
        'is_critical': False,
        'status': 'Unknown'
    },
}


class VaspParser(BaseParser):
    """
    Parses all Vasp calculations.

    This particular class manages all the specific file parsers in
    aiida_vasp.io. The parser will check which quantities to parse and which nodes to add
    to the calculation based on the 'parser_settings' card in the 'settings' ParameterData of the
    corresponding VaspCalculation.

    Parser Settings usage:

    Parser settings can be passed through the input node `settings` as follows::

        settings = ParameterData(dict={
            'parser_settings': {
                ...
            }
        })

    Valid keys for `parser_settings` are:

    * `add_<quantity>`, where quantity is one of:

        'parameters': Parameterdata node containing various quantities from OUTCAR and vasprun.xml.
        'structure':  (Default) StructureData node parsed from CONTCAR
        'bands':      Band structure node parsed from EIGENVAL.
        'dos':        ArrayData node containing the DOS parsed from DOSCAR.
        'kpoints':    KpointsData node parsed from IBZKPT.
        'wavecar':    FileData node containing the WAVECAR file.
        'chgcar':     FileData node containing the CHGCAR file.

        Set any of these to True, if the corresponding output node should be added.

    * `output_params`: A list of things that should be added to the 'output_parameters'
                       node, available options are:
                       DEFAULT: ['energies', 'forces']
    """

    def __init__(self, calc):
        super(VaspParser, self).__init__(calc)

        self.out_folder = None

        self._parsers = {}
        self._parsable_quantities = {}

        # Gather all parsable items as defined in the file parsers.
        for filename, value in PARSABLE_FILES.iteritems():
            self._parsable_quantities.update(value['parser_class'].PARSABLE_ITEMS)
            self._parsers[filename] = None

        self._settings = DEFAULT_OPTIONS
        calc_settings = self._calc.get_inputs_dict().get('settings')
        if calc_settings:
            self._settings.update(calc_settings.get_dict().get('parser_settings', DEFAULT_OPTIONS))

        self._quantities_to_parse = []
        self._output_nodes = {}
        self._requested_quantities = []

        self._check_and_validate_settings()

    def parse_with_retrieved(self, retrieved):

        self.check_state()
        self.out_folder = self.get_folder(retrieved)
        
        if not self.out_folder:
            return self.result(success=False)

        # Get all specialised file parsers. Warnings will be issued if a file should be parsed and
        # the corresponding files do not exist.
        success = self._set_file_parsers()

        if not success:
            # A critical file i.e. OUTCAR does not exist. Abort parsing.
            return self.result(success=False)

        # Parse all implemented quantities in the quantities_to_parse list.
        print("quantities_to_parse:", self._quantities_to_parse)
        while self._quantities_to_parse:
            quantity = self._quantities_to_parse.pop(0)
            self._output_nodes.update(self.get_quantity(quantity, self._settings))
            
        # Add output nodes if the corresponding data exists.
        for key, value in self._output_nodes.iteritems():
            print("KEY, VALUE: ", key, value)
            if key != self._parsable_quantities[key]['nodeName']:
                # this is just an intermediate result and should not be added as a node.
                continue

            if value:
                self._set_node(key, value)

        return self.result(success=True)

    def _check_and_validate_settings(self):
        """Check the settings and set which files should be parsed based on the input."""

        import copy
        new_settings = copy.deepcopy(self._settings)

        calc_settings = self._calc.get_inputs_dict().get('settings')
        calc_parser_settings = calc_settings.get_dict().get('parser_settings')

        for key, value in self._settings.iteritems():
            if not key.startswith('add_'):
                # only keys starting with 'add_' will change the behaviour of the parser so get the next one.
                continue
            if not value:
                # The quantity should not be added, so the corresponding files do not have to be parsed.
                continue
            quantity = key[4:]
            if quantity not in self._parsable_quantities:
                self.logger.warning('{0} has been requested by setting add_{0}'.format(quantity) +
                                    ' however it has not been implemented. Please check the docstrings' +
                                    ' in  aiida_vasp.parsers.vasp.py for valid input.')
                continue

            # Found a node, which should be added, add itself and all it's perequisites to the quantities to parse.
            for quantity in self._parsable_quantities[quantity]['prerequisites'] + [quantity]:
                if quantity in self._quantities_to_parse:
                    continue
                self._quantities_to_parse.append(quantity)
                # Flag all the required files for being parsed.
                for filename in self._parsable_quantities[quantity]['parsers']:
                    k = 'should_parse_' + filename
                    new_settings[k] = value
                    # make sure we rewrite if this was overridden for the
                    # calc settings
                    # eFL: THIS IS SUPER NASTY AND SHOULD BE REPLACED ASAP
                    #try:
                    #    new_settings[k] = calc_parser_settings[k]
                    #except KeyError:
                    #    pass

        self._settings = new_settings

    def _set_file_parsers(self):
        """
        Set the specific file parsers for OUTCAR, DOSCAR, EIGENVAL and vasprun.xml.

        Return False if a critical file is missing, which will abort the parsing.
        """

        for file_name, value in PARSABLE_FILES.iteritems():
            if not self._settings['should_parse_' + file_name]:
                # Based on the settings, we should not parse this file. continue with the next one.
                continue
            if self._parsers[file_name] is not None:
                # This fileParser has already been set, continue with the next one.
                continue

            # We should parse this file and the parser has not been set yet.
            file_to_parse = self.get_file(file_name)
            if not file_to_parse:
                self._parsers[file_name] = None
                if value['is_critical']:
                    self.logger.error('{} not found, look at the scheduler output for troubleshooting.'.format(file_name))
                    return False

                self.logger.warning('{0} not found, but should be parsed.'.format(file_name))
            else:
                # The file should be parsed and has been found
                self._parsers[file_name] = value['parser_class'](self, file_path=file_to_parse)

        # All critical files have been found, so we can safely return True.
        return True

    # pylint: disable=unused-argument, no-self-use
    @delegate()
    def get_quantity(self, quantity, settings):
        """Delegate to which the FileParsers will subscribe their get_quantity methods when they are initialised."""
        return {quantity: None}

    def get_inputs(self, quantity):
        """
        Return a quantity required as input for another quantity.

        This method will be called by the FileParsers in order to get a required input quantity
        from self._output_nodes. If the quantity is not in the dictionary the VaspParser will
        try to parse it. If a quantiy has been requested this way two times, parsing will be
        aborted because there is a cyclic dependency of the parsable items.
        """

        if quantity in self._requested_quantities:
            raise RuntimeError('{0} has been requested for parsing a second time.'.format(quantity) +
                               ' There is probably a cycle in the prerequisites of the parsable_items in the single FileParsers.')

        # This is the first time this quantity has been requested, keep track of it.
        self._requested_quantities.append(quantity)

        if quantity not in self._output_nodes:
            # The quantity is not in the output_nodes. Try to parse it
            self._output_nodes.update(self.get_quantity(quantity, self._settings))

        # parsing the quantity without requesting it a second time was succesful, so remove it from requested_quantities.
        self._requested_quantities.remove(quantity)

        # since the quantity has already been parsed as an input now, we don't have to parse it a second time later.
        if quantity in self._quantities_to_parse:
            self._quantities_to_parse.remove(quantity)

        return self._output_nodes.get(quantity)

    def _set_node(self, node_name, node):
        """Wrapper for self.add_node, checking whether the Node is None and using the correct linkname"""

        if node is not None:
            self.add_node(LINKNAME_DICT[node_name], node)
