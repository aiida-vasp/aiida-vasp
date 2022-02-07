"""
VASP parser.

------------
The main driver routine for parsing VASP related objects. The parser is
modular and contains several modules:

- ``node_composer`` handles the quantity composition of nodes
- ``quantity`` the actual quantity to parse and what object parsers to use to obtain it
- ``settings`` general parser settings
- ``manager`` takes the quantity definitions and executes the actual parsing needed
"""
#encoding: utf-8
# pylint: disable=no-member
import traceback

from aiida.common.exceptions import NotExistent
from aiida_vasp.parsers.base import BaseParser
from aiida_vasp.parsers.quantity import ParsableQuantities
from aiida_vasp.parsers.settings import ParserSettings, ParserDefinitions
from aiida_vasp.parsers.node_composer import NodeComposer

DEFAULT_SETTINGS = {
    'add_trajectory': False,
    'add_bands': False,
    'add_charge_density': False,
    'add_dos': False,
    'add_kpoints': False,
    'add_energies': False,
    'add_misc': True,
    'add_structure': False,
    'add_projectors': False,
    'add_born_charges': False,
    'add_dielectrics': False,
    'add_hessian': False,
    'add_dynmat': False,
    'add_wavecar': False,
    'add_forces': False,
    'add_stress': False,
    'add_site_magnetization': False,
    'critical_notifications': {
        'add_brmix': True,
        'add_cnormn': True,
        'add_denmp': True,
        'add_dentet': True,
        'add_edddav_zhegv': True,
        'add_eddrmm_zhegv': True,
        'add_edwav': True,
        'add_fexcp': True,
        'add_fock_acc': True,
        'add_non_collinear': True,
        'add_not_hermitian': True,
        #add_psmaxn': True,
        'add_pzstein': True,
        'add_real_optlay': True,
        'add_rhosyg': True,
        'add_rspher': True,
        'add_set_indpw_full': True,
        'add_sgrcon': True,
        'add_no_potimm': True,
        'add_magmom': True,
    }
}


class VaspParser(BaseParser):
    """
    Parses all VASP calculations.

    This particular class manages all the specific parsers located in
    aiida_vasp.parsers.content_parsers. The parser will check which quantities to parse
    and which nodes to add to the calculation based on the 'parser_settings' card in
    the 'settings' Dict of the corresponding VaspCalculation.

    Parser Settings usage:

    Parser settings can be passed through the input node `settings` as follows::

        settings = Dict(dict={
            'parser_settings': {
                ...
            }
        })

    Valid keys for `parser_settings` are:

    * `add_<quantity>`, where quantity is one of:

        'misc': (Default) Parameterdata node containing various quantities from OUTCAR and vasprun.xml.
        'structure':  (Default) StructureData node parsed from CONTCAR
        'bands':      Band structure node parsed from EIGENVAL.
        'dos':        ArrayData node containing the DOS parsed from DOSCAR.
        'kpoints':    KpointsData node parsed from IBZKPT.
        'wavecar':    SinglefileData node containing the WAVECAR.
        'charge_density':     SinglefileData node containing the CHGCAR.
       If the value is set to ``False`` the quantity will not be returned.

    * `critical_notifications`: A dictionary of critical errors to be checked with items like `'add_<key>': True`, similiar
      to the `add_<quantity>` syntax described above.

    * `output_params`: A list of quantities, that should be added to the 'misc' node.

    * `object_parser_set`: String (DEFAULT = 'default').

        By this option the default set of object parsers can be chosen. See settings.py
        for available options.

    * `ignore_all_errors`: If set to `True`, will skip checks for critical error messages. Defaults to `False`.

    Additional object parsers can be added to the VaspParser by using

        VaspParser.add_object_parser(parser_name, parser_definition_dict),

    where the 'parser_definition_dict' should contain the 'parser_class' and the
    'is_critical' flag. Keep in mind adding an additional object parsers after 'parse_with_retrieved'
    is called, will only have an effect when parsing a second time.
    """

    def __init__(self, node):
        super(VaspParser, self).__init__(node)

        try:
            calc_settings = self.node.inputs.settings
        except NotExistent:
            calc_settings = None

        parser_settings = None
        if calc_settings:
            parser_settings = calc_settings.get_dict().get('parser_settings')

        self._definitions = ParserDefinitions()
        self._settings = ParserSettings(parser_settings, default_settings=DEFAULT_SETTINGS, vasp_parser_logger=self.logger)
        self._parsable_quantities = ParsableQuantities(vasp_parser_logger=self.logger)

    def add_parser_definition(self, name, parser_dict):
        """Add the definition of an oobject parser to self._definitions."""
        self._definitions.add_parser_definition(name, parser_dict)

    def add_parsable_quantity(self, quantity_name, quantity_dict):
        """Add a single parsable quantity to the _parsable_quantities."""
        self._parsable_quantities.add_parsable_quantity(quantity_name, quantity_dict)

    def add_custom_node(self, node_name, node_dict):
        """Add a custom node to the settings."""
        self._settings.add_output_node(node_name, node_dict)

    def parse(self, **kwargs):  # pylint: disable=too-many-return-statements
        """The function that triggers the parsing of a calculation."""

        error_code = self._compose_retrieved_content(kwargs)
        if error_code is not None:
            return error_code

        for name, value_dict in self._definitions.parser_definitions.items():
            if name not in self._retrieved_content.keys() and value_dict['is_critical']:
                self.logger.error('Missing content: {} which is tagged as critical by the parser'.format(name))
                return self.exit_codes.ERROR_CRITICAL_MISSING_OBJECT
        self._parsable_quantities.setup(retrieved_content=self._retrieved_content.keys(),
                                        parser_definitions=self._definitions.parser_definitions,
                                        quantity_names_to_parse=self._settings.quantity_names_to_parse)

        # Update the parser settings to make sure that the quantities that have been requested from
        # the collection of the nodes are included. Quantities already present in settings are preserved.
        self._settings.update_quantities_to_parse(self._parsable_quantities.quantity_keys_to_parse)

        # Parse the quantities from retrived objects
        parsed_quantities, failed_to_parse_quantities = self._parse_quantities()

        # Compose the output nodes using the parsed quantities
        requested_nodes = self._settings.output_nodes_dict
        equivalent_quantity_keys = dict(self._parsable_quantities.equivalent_quantity_keys)
        composed_nodes = NodeComposer(requested_nodes, equivalent_quantity_keys, parsed_quantities, logger=self.logger)
        for link_name, node in composed_nodes.successful.items():
            self.out(link_name, node)

        nodes_failed_to_create = composed_nodes.failed

        # Check for execution related errors
        exit_code = self._check_vasp_errors(parsed_quantities)
        if exit_code is not None:
            return exit_code

        # Deal with missing quantities
        if failed_to_parse_quantities:
            return self.exit_codes.ERROR_NOT_ABLE_TO_PARSE_QUANTITY.format(quantity=', '.join(failed_to_parse_quantities))

        # Deal with missing node/nodes
        if nodes_failed_to_create:
            return self.exit_codes.ERROR_NOT_ABLE_TO_CREATE_NODE.format(nodes=', '.join(nodes_failed_to_create))

        return self.exit_codes.NO_ERROR

    def _parse_quantities(self):
        """
        This method dispatch the parsing to object parsers

        :returns: A tuple of parsed quantities dictionary and a list of quantities failed to obtain due to exceptions
        """
        parsed_quantities = {}
        # A dictionary for catching instantiated object parser objects
        object_parser_instances = {}
        failed_to_parse_quantities = []
        for quantity_key in self._parsable_quantities.quantity_keys_to_parse:
            name = self._parsable_quantities.quantity_keys_to_content[quantity_key]
            object_parser_cls = self._definitions.parser_definitions[name]['parser_class']
            # If a parsed object has been instantiated, use it.
            if object_parser_cls in object_parser_instances:
                parser = object_parser_instances[object_parser_cls]
            else:
                try:
                    # The next line may except for ill-formated object
                    with self._get_handler(name) as handler:
                        parser = object_parser_cls(settings=self._settings.settings, handler=handler)
                except Exception:  # pylint: disable=broad-except
                    parser = None
                    failed_to_parse_quantities.append(quantity_key)
                    self.logger.warning('Cannot instantiate {}, exception {}:'.format(object_parser_cls, traceback.format_exc()))

                object_parser_instances[object_parser_cls] = parser

            if parser is None:
                # If the parser cannot be instantiated, add the quantity to a list of unavailable ones
                failed_to_parse_quantities.append(quantity_key)
                continue
            exception = None
            try:
                # The next line may still except for ill-formated object - some parser load all data at
                # instantiation time, the others may not.
                parsed_quantity = parser.get_quantity(quantity_key)
            except Exception:  # pylint: disable=broad-except
                parsed_quantity = None
                exception = traceback.format_exc()

            if parsed_quantity is not None:
                parsed_quantities[quantity_key] = parsed_quantity
            else:
                self.logger.warning('Parsing {} from {} failed, exception: {}'.format(quantity_key, parser, exception))
                failed_to_parse_quantities.append(quantity_key)

        return parsed_quantities, failed_to_parse_quantities

    @property
    def parser_settings(self):
        """The `parser_settings` dictionary passed"""
        return self._settings._settings  # pylint: disable=protected-access

    @property
    def _check_ionic_convergence(self):
        """
        Wether to check the ionic convergence

        This can be customised using flag in the settings of the calculation

        Usage::

          builder.settings = Dict(dict={
              'CHECK_IONIC_CONVERGENCE': True
          })

        The default is `True` so a calculation that has ran for NSW steps is treated
        as not converged.
        """

        if 'settings' in self.node.inputs:
            settings = self.node.inputs.settings.get_dict()
        else:
            settings = {}
        return settings.get('CHECK_IONIC_CONVERGENCE', True)

    def _check_vasp_errors(self, quantities):
        """
        Detect simple vasp execution problems and returns the exit_codes to be set
        """

        if 'run_status' not in quantities:
            return self.exit_codes.ERROR_DIAGNOSIS_OUTPUTS_MISSING
        run_status = quantities['run_status']

        # Return errors related to execution and convergence problems.
        # Note that the order is important here - if a calculation is not finished, we cannot
        # comment on wether properties are converged are not.
        if run_status['finished'] is False:
            return self.exit_codes.ERROR_DID_NOT_FINISH

        if run_status['electronic_converged'] is False:
            return self.exit_codes.ERROR_ELECTRONIC_NOT_CONVERGED

        # Check the ionic convergence issues
        if run_status['ionic_converged'] is False:
            if self._check_ionic_convergence:
                return self.exit_codes.ERROR_IONIC_NOT_CONVERGED
            self.logger.warning('The ionic relaxation is not converged, but the calculation is treated as successful.')

        # Check for the existence of critical warnings
        if 'notifications' in quantities:
            notifications = quantities['notifications']
            ignore_all = self.parser_settings.get('ignore_all_errors', False)
            if not ignore_all:
                composer = NotificationComposer(notifications,
                                                quantities,
                                                self.node.inputs,
                                                self.exit_codes,
                                                parser_settings=self._settings)
                exit_code = composer.compose()
                if exit_code is not None:
                    return exit_code
        else:
            self.logger.warning('WARNING: missing notification output for VASP warnings and errors.')

        return None


class NotificationComposer:
    """Compose errors codes based on the notifications"""

    def __init__(self, notifications, parsed_quantities, inputs, exit_codes, parser_settings):
        """
        Composed error codes based on the notifications

        Some of the errors need to have additional properties inspected before they can be emitted,
        as they might be trigged in a harmless way.

        To add new checkers, one needs to implement a property with the name of the error for this class and
        contains the code for checking. This property should return the exit_code or return None. The property
        is inspected if its name is in the list critical notifications.

        :param notification: The list of parsed notifications from the stream parser.
        :param parsed_quantities: The dictionary of parsed quantities.
        :param inputs: The dictionary of the input nodes.
        :param exit_codes: The dictionary of the exit codes from the parser.
        :param ignored: A list of critical notification that are allowed to have
        """
        self.notifications = notifications
        self.notifications_dict = {item['name']: item['message'] for item in self.notifications}
        self.parsed_quantities = parsed_quantities
        self.inputs = inputs
        self.exit_codes = exit_codes
        self.parser_settings = parser_settings

    def compose(self):
        """
        Compose the exit codes

        Retruns None if no exit code should be emitted, otherwise emit the error code.
        """
        for critical in self.parser_settings.critical_notifications_to_check:
            # Check for any special handling
            if hasattr(self, critical):
                output = getattr(self, critical)
                if output:
                    return output
            # No special handling, just check if it exists
            elif critical in self.notifications_dict:
                return self.exit_codes.ERROR_VASP_CRITICAL_ERROR.format(error_message=self.notifications_dict[critical])
        return None

    @property
    def brmix(self):
        """Check if BRMIX should be emitted"""
        if not 'brmix' in self.notifications_dict:
            return None

        # If NELECT is set explicitly for the calculation then this is not an critical error
        if 'parameters' in self.inputs and 'nelect' in self.inputs['parameters'].get_dict():
            return None

        return self.exit_codes.ERROR_VASP_CRITICAL_ERROR.format(error_message=self.notifications_dict['brmix'])
