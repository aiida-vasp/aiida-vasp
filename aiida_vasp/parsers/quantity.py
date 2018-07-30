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
        self.parsers = []
        self.alternatives = []

        super(ParsableQuantity, self).__init__(init)
        self.name = name

        # Check whether all files required for parsing this quantity have been retrieved and store it.
        missing_files = self.has_all(files_list)
        if missing_files is None:
            self.has_files = False
            missing_files = []
        else:
            self.has_files = not missing_files
        self.missing_files = missing_files

        # Check whether this quantity represents a node
        self.is_node = self.nodeName in NODES

    def has_all(self, available_items):
        """Check whether all items are in item_list."""
        missing_items = []
        if not self.parsers:
            return None
        if not available_items:
            return self.parsers
        for item in self.parsers:
            if item not in available_items:
                missing_items.append(item)
        return missing_items


class ParsableQuantities(object):
    """
    A Database of parsable quantities.

    - Get the parsable quantities from the FileParsers and initialise them.
    - Provide ways to get parsable quantities.
    """

    def __init__(self, vasp_parser=None):
        self._quantities = {}

        self._vasp_parser = vasp_parser
        self._parsers = None

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

        if self._parsers is None:
            self._parsers = self._vasp_parser.parsers

        retrieved = self._vasp_parser.out_folder.get_folder_list()
        # check uniqueness and add parsable quantities
        self._check_uniqueness_add_parsable(retrieved)

        # check consistency, that the quantity is parsable and
        # alternatives
        self._check_consitency_and_alternatives()

    def _check_uniqueness_add_parsable(self, retrieved):
        """Check uniqueness and add parsable quantities."""

        self._quantities = {}
        # Gather all parsable items as defined in the file parsers.
        for filename, value in self._parsers.get_parsers():
            for quantity, quantity_dict in value['parser_class'].PARSABLE_ITEMS.items():
                if quantity in self._quantities:
                    # This quantity has already been added so it is not unique.
                    raise RuntimeError('The quantity {quantity} defined in {filename} has been '
                                       'defined by two FileParser classes. Quantity names must '
                                       'be unique. If both quantities are equivalent, define one '
                                       'as an alternative for the other.'.format(quantity=quantity, filename=filename))
                # Create quantity objects.
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
                if not self._quantities[prereq].has_files:
                    is_parsable = False
            value.is_parsable = is_parsable and value.has_files

            # Add this quantity to the list of alternatives of another quantity.
            if value.is_alternative is not None:
                if value.is_alternative not in self._quantities:
                    # The quantity which this quantity is an alternative to is not in _parsable_quantities.
                    # Add an unparsable dummy quantity for it.
                    self.add_parsable_quantity(value.is_alternative, {'alternatives': [], 'nodeName': self._quantities[quantity].nodeName})

                if quantity not in self._quantities[value.is_alternative].alternatives:
                    self._quantities[value.is_alternative].alternatives.append(quantity)
