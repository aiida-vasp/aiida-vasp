#encoding: utf-8
"""AiiDA Parser for a aiida_vasp.VaspCalculation"""

from aiida_vasp.io.doscar import DosParser
from aiida_vasp.io.eigenval import EigParser
from aiida_vasp.io.kpoints import KpParser
from aiida_vasp.io.outcar import OutcarParser
from aiida_vasp.io.vasprun import VasprunParser
from aiida_vasp.io.chgcar import ChgcarParser
from aiida_vasp.io.wavecar import WavecarParser
from aiida_vasp.parsers.base import BaseParser

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
    'should_parse_DOSCAR': False,
    'should_parse_EIGENVAL': False,
    'should_parse_IBZKPT': False,
    'should_parse_OUTCAR': True,
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
}


class VaspParser(BaseParser):
    """
    Parses all Vasp calculations. This particular class manages all the specific file parsers in
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

        # Gather all parsable items as defined in the file parsers.
        for filename in PARSABLE_FILES:
            self._parsable_quantities.update(filename['parser_class'].PARSABLE_ITEMS)
            self._parsers[filename] = None

        self._settings = DEFAULT_OPTIONS
        self._settings.update(self._calc.inp.settings.get_dict().get('parser_settings', DEFAULT_OPTIONS))

        self._check_and_validate_settings()

        self._quantities_to_parse = []
        self._output_nodes = {}

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

        # Parse all implemented quantities in the nodesToAdd list, if they should be parsed. The list
        # might get dynamically updated during the loop.
        while self._quantities_to_parse:
            quantity = self._quantities_to_parse.pop(0)
            if self._settings['add_' + quantity]:
                if not self._check_prerequisites(quantity):
                    continue
                for parser in self._parsable_quantities['quantity']['parsers']:
                    parser.get_quantities(quantity, self._output_nodes)

        # Add output nodes if the corresponding data exists.
        for key, value in self._output_nodes.iteritems():

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
                    new_settings['should_parse_' + filename] = value

        self._settings = new_settings

    def _set_file_parsers(self):
        """
        Set the specific file parsers for OUTCAR, DOSCAR, EIGENVAL and vasprun.xml.
        Return False if a critical file is missing, which will abort the parsing.
        """

        for key, value in PARSABLE_FILES.iteritems():
            if not self._settings['should_parse_' + key]:
                continue
            if self._parsers[key]:
                continue

            # We should parse this file and the parser has not been set yet.
            file_to_parse = self.get_file(key)
            if not file_to_parse:
                self._parsers[key] = None
                if value['is_critical']:
                    self.logger.error('{} not found, look at the scheduler output for troubleshooting.'.format(key))
                    return False

                self.logger.warning('{0} not found, but should be parsed.'.format(key))
            else:
                # The file should be parsed and has been found
                self._parsers[key] = value['parser_class'](file_to_parse, key)

        # All critical files have been found, so we can safely return True.
        return True

    def _check_prerequisites(self, quantity):
        """Check whether the prerequesites of a given quantity have been met. If not either
           requeue or prevent this quantity from being parsed."""

        prerequesites = self._parsable_quantities[quantity]['prerequesites']
        for preq in prerequesites:
            if preq in self._output_nodes:
                # requirement met, check the next one
                continue

            # Requirement not met yet.
            if preq in self._quantities_to_parse:
                # The prerequesite is in the queue, requeue this quantity and return
                self._quantities_to_parse.append(quantity)
                return False

            # The prerequesite is not met and also not in the queue. Don't parse this quantity.
            return False

        # All requirements have been met
        return True

    def _should_parse_quantity(self, quantity):
        """
        Check whether a quantity should be parsed based on whether all required parsers are
        available.
        """
        for filename in self._parsable_quantities[quantity]['parsers']:
            if not self._parsers[filename]:
                return False

        return self._settings.get('add_' + self._parsable_quantities[quantity]['nodename'])

    def _set_node(self, node_name, node):
        """Wrapper for self.add_node, checking whether the Node is None and using the correct linkname"""

        if node is not None:
            self.add_node(LINKNAME_DICT[node_name], node)
