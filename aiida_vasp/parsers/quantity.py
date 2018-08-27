"""Classes representing quantities and a list of quantities."""

from aiida_vasp.utils.extended_dicts import DictWithAttributes

NODES = {
    'parameters': 'output_parameters',
    'kpoints': 'output_kpoints',
    'structure': 'output_structure',
    'trajectory': 'output_trajectory',
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


class ParsableQuantity(DictWithAttributes):
    """Container class for parsable quantities."""

    def __init__(self, name, init, files_list):
        # assign default values for optional attributes
        self.alternatives = []

        super(ParsableQuantity, self).__init__(init)
        self.name = name

        # Check whether the file required for parsing this quantity have been retrieved.
        missing_files = []

        if files_list is None or self.file_name not in files_list:
            missing_files.append(self.file_name)
        self.missing_files = missing_files

        # Check whether this quantity represents a node
        self.is_node = self.nodeName in NODES


class ParsableQuantities(object):
    """
    A Database of parsable quantities.

    - Get the parsable quantities from the FileParsers and initialise them.
    - Provide ways to get parsable quantities.
    """

    def __init__(self, vasp_parser=None):
        self._quantities = {}
        self.additional_quantities = {}

        self._vasp_parser = vasp_parser

    def add_parsable_quantity(self, quantity_name, quantity_dict, retrieved_files=None):
        self._quantities[quantity_name] = ParsableQuantity(quantity_name, quantity_dict, retrieved_files)

    def remove_parsable_quantity(self, quantity_name):
        _ = self._quantities.pop(quantity_name, None)

    def get_equivalent_quantities(self, quantity_name):
        """Get a list of equivalent quantities."""
        quantity = self._quantities[quantity_name]
        return [quantity] + [self._quantities[item] for item in quantity.alternatives]

    def get_by_name(self, quantity_name):
        """Get a quantity by name."""
        return self._quantities.get(quantity_name)

    def get_missing_files(self, quantity_name):
        """Return a list with all missing files for a quantity."""
        missing_files = []
        for quantity in self.get_equivalent_quantities(quantity_name):
            for missing_file in quantity.missing_files:
                missing_files.append(missing_file)

        return missing_files

    def setup(self):
        """Set the parsable_quantities dictionary based on parsable_items obtained from the FileParsers."""

        retrieved = self._vasp_parser.out_folder.get_folder_list()

        # check uniqueness and add parsable quantities
        self._check_uniqueness_add_parsable(retrieved)

        # check consistency, that the quantity is parsable and
        # alternatives
        self._check_consitency_and_alternatives()

    def _check_uniqueness_add_parsable(self, retrieved):
        """Check uniqueness and add parsable quantities."""

        self._quantities = {}

        # Add all the additional quantities that have been added by the user.
        for key, value in self.additional_quantities.items():
            self.add_parsable_quantity(key, value, retrieved)

        # Gather all parsable items as defined in the file parsers.
        for file_name, value in self._vasp_parser.settings.parser_definitions.items():
            for quantity, quantity_dict in value['parser_class'].PARSABLE_ITEMS.items():
                if quantity in self._quantities:
                    # This quantity has already been added so it is not unique.
                    raise RuntimeError('The quantity {quantity} defined in {filename} has been '
                                       'defined by two FileParser classes. Quantity names must '
                                       'be unique. If both quantities are equivalent, define one '
                                       'as an alternative for the other.'.format(quantity=quantity, filename=file_name))
                # Create quantity objects.
                quantity_dict['file_name'] = file_name
                self.add_parsable_quantity(quantity, quantity_dict, retrieved)

    def _check_consitency_and_alternatives(self):
        """Check the consistency and alternatives."""

        import copy

        # Make a local copy of parsable_quantities, because during the next step
        # dummy quantities for missing quantities might be added.
        parsable_quantities = copy.deepcopy(self._quantities.keys())

        # Check for every quantity, whether all of the prerequisites are parsable, and whether
        # they are alternatives for another quantity.
        for quantity in parsable_quantities:
            value = self._quantities[quantity]
            is_parsable = True
            for prereq in value.prerequisites:
                if self._quantities[prereq].missing_files:
                    is_parsable = False
            value.is_parsable = is_parsable and not value.missing_files

            # Add this quantity to the list of alternatives of another quantity.
            if value.is_alternative is not None:
                if value.is_alternative not in self._quantities:
                    # The quantity which this quantity is an alternative to is not in _parsable_quantities.
                    # Add an unparsable dummy quantity for it.
                    self.add_parsable_quantity(value.is_alternative, {'alternatives': [], 'nodeName': self._quantities[quantity].nodeName})

                if quantity not in self._quantities[value.is_alternative].alternatives:
                    self._quantities[value.is_alternative].alternatives.append(quantity)
