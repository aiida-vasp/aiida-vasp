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
from aiida_vasp.utils.extended_dicts import DictWithAttributes

LINKNAME_DICT = {
    'parameters': 'output_parameters',
    'kpoints': 'output_kpoints',
    'structure': 'output_structure',
    'array': 'output_array',
    'trajectory': 'output_trajectory',
    'bands': 'output_band',
    'dos': 'output_dos',
    'chgcar': 'chgcar',
    'wavecar': 'wavecar',
    'born_charges': 'born_charges',
}

DEFAULT_OPTIONS = {
    'add_bands': False,
    'add_chgcar': False,
    'add_dos': False,
    'add_kpoints': False,
    'add_parameters': True,
    'add_structure': True,
    'add_wavecar': False,
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
        'is_critical': True,
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

QUANTITY_ATTR_DEFAULTS = {
    'alternatives': [],
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
    """

    def __init__(self, calc):
        super(VaspParser, self).__init__(calc)

        self.out_folder = None

        self._parsers = {}
        self._parsable_quantities = {}

        self._settings = DEFAULT_OPTIONS
        calc_settings = self._calc.get_inputs_dict().get('settings')
        if calc_settings:
            self._settings.update(calc_settings.get_dict().get('parser_settings', DEFAULT_OPTIONS))

        self._quantities_to_parse = []
        self._output_nodes = {}

        # this list is for bookkeeping, to check whether a quantity has been requested
        # twice during the parsing cycle.
        self._requested_quantities = []

    def parse_with_retrieved(self, retrieved):

        def missing_critical_file():
            for file_name, value_dict in PARSABLE_FILES.iteritems():
                if file_name not in self.out_folder.get_folder_list() and value_dict['is_critical']:
                    return True
            return False

        self.check_state()
        self.out_folder = self.get_folder(retrieved)

        if not self.out_folder:
            return self.result(success=False)

        if missing_critical_file():
            # A critical file i.e. OUTCAR does not exist. Abort parsing.
            return self.result(success=False)

        # Get the _parsable_quantities from the FileParsers.
        self._set_parsable_quantities()

        # Set the quantities to parse list. Warnings will be issued if a quantity should be parsed and
        # the corresponding files do not exist.
        self._check_and_validate_settings()

        self._set_file_parsers()

        # Parse all implemented quantities in the quantities_to_parse list.
        while self._quantities_to_parse:
            quantity = self._quantities_to_parse.pop(0)
            self._output_nodes.update(self.get_quantity(quantity, self._settings))

        # Add output nodes if the corresponding data exists.
        for key, value in self._output_nodes.iteritems():

            if key != self._parsable_quantities[key]['nodeName']:
                # this is just an intermediate result and should not be added as a node.
                continue

            if value:
                self._set_node(key, value)

        return self.result(success=True)

    def _set_parsable_quantities(self):
        """Set the parsable_quantities dictionary based on parsable_items obtained from the FileParsers."""

        def have_all(item_list, available_items):
            """Check whether all items are in item_list."""
            missing_items = []
            for item in item_list:
                if item not in available_items:
                    missing_items.append(item)
            return missing_items

        # Gather all parsable items as defined in the file parsers.
        for filename, value in PARSABLE_FILES.iteritems():
            self._parsers[filename] = None
            for quantity, quantity_dict in value['parser_class'].PARSABLE_ITEMS.iteritems():

                # Create quantity objects.
                self._parsable_quantities[quantity] = DictWithAttributes(quantity_dict, QUANTITY_ATTR_DEFAULTS)
                self._parsable_quantities[quantity].name = quantity

                # Check whether all files required for parsing this quantity have been retrieved and store it.
                missing_files = have_all(quantity_dict.get('parsers', []), self.out_folder.get_folder_list())
                self._parsable_quantities[quantity].have_files = not missing_files
                self._parsable_quantities[quantity].missing_files = missing_files

    def _check_and_validate_settings(self):
        """Check the settings and set which files should be parsed based on the input."""

        def add_quantity(quantity_to_add):
            """Check, whether a quantity or it's alternatives can be added."""
            for item in [quantity_to_add] + self._parsable_quantities[quantity_to_add].alternatives:
                if self._parsable_quantities[item].have_files:
                    self._quantities_to_parse.append(item)
                    return True
            return False

        for key, value in self._settings.iteritems():
            if not key.startswith('add_'):
                # only keys starting with 'add_' will change the behaviour of the parser so get the next one.
                continue
            if not value:
                # The quantity should not be added, so the corresponding files do not have to be parsed.
                continue
            quantity = key[4:]
            if quantity not in self._parsable_quantities:
                self.logger.warning(
                    '{0} has been requested by setting add_{0}'.format(quantity) +
                    ' however it has not been implemented. Please check the docstrings' + ' in aiida_vasp.parsers.vasp.py for valid input.')
                continue

            # Found a node, which should be added, add itself and all it's perequisites to the quantities to parse.
            # if all files required for this quantity have been retrieved. If there are alternatives for this quantity
            # also try those.
            success = add_quantity(quantity)

            if not success:
                # Neither the quantity or it's alternatives could be added to the quantities_to_parse.
                # Gather a list of all the missing files and issue a warning.
                missing_files = []
                for quant in [quantity] + self._parsable_quantities[quantity].alternatives:
                    for missing_file in quant.missing_files:
                        missing_files.append(missing_file)

                missing_files = ", ".join(missing_files)
                self.logger.warning('{0} has been requested,'.format(quantity) + ' however the following files required for parsing' +
                                    ' have not been retrieved: {0}.'.format(missing_files))

    def _set_file_parsers(self):
        """
        Set the specific file parsers for OUTCAR, DOSCAR, EIGENVAL and vasprun.xml.

        Return False if a critical file is missing, which will abort the parsing.
        """

        for quantity in self._quantities_to_parse:
            for filename in self._parsable_quantities[quantity]['parsers']:
                if self._parsers[filename] is not None:
                    # This parser has already been checked.
                    continue
                file_to_parse = self.get_file(filename)
                self._parsers[filename] = PARSABLE_FILES[filename]['parser_class'](file_to_parse, calc_parser_cls=self)

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
