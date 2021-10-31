"""
Base classes for the VASP object parsers.

---------------------------------------
Contains the base classes for the VASP content parsers.
"""
# pylint: disable=import-outside-toplevel
import re
from aiida.common import AIIDA_LOGGER as aiidalogger
from aiida_vasp.utils.delegates import delegate_method_kwargs


class BaseFileParser(object):
    """
    Base class for all the content parsers which parse (read and write) VASP files.

    The actually read, interpret and write to/from files is handled by
    parsevasp and this interface ensures that the parsing framework in AiiDA and
    preparation before submission is successful. The specific content parser
    interfaces are always childs of this class.

    We assume that there are two main paths when parsing of VASP files takes place.
    One is when a file is present and we want to interpret it and convert its
    content into one of the usable data structures in AiiDA.
    The second is when we already have an AiiDA data structure  and want to write a
    file based on its content. In the former case we basically
    read from a file, while in the latter we write to a file.

    The first approach is enabled if we initialize the parser with the argument `handler`.
    This should be a file like handler, i.e. from a context manager telling where the
    content can be located. The respective quantities can then be fetched using the
    `get_quantity` function of the instance with the key representing the quantity.
    The content of the handler is parsed after
    initialization. The valid keys representing fetchable quantities is defined for each
    content parser class using the `PARSABLE_QUANTITIES` class parameter.

    The second approach is enabled if we initialize with an argument `data`. This should be
    a valid AiiDA data structure node. Using the `get_quantity('somekey')` function of the instance
    return the same AiiDA data structure node back as was supplied to `data`.

    One can write the parsed content or the content of the StructureData using the
    function `write` of the class instance.

    Parameters
    ----------
    handler : object, optional
        A file like object, for instance a file handler representing the file or object
        containing content to be parsed. Typically used when parsing completed calculations and is
        also the parameter used during initialization of the content parser, when the `parse`
        function of this AiiDA plugin is executed.
    data : object. optional
        An AiiDA data structure node. Typically used when one later want to write VASP input
        files.
    settings : dict
        Parser settings. Used to set parser settings, e.g. which quantities to compose into nodes etc.
    options : dict
        Parser options. Used to set extra options to the content parsers. For instance for the
        POSCAR/CONTAR parser one set `options.positions_dof` to supply selective tags to enable proper
        construction of selective dynamics POSCAR/CONTCAR from a StructureData. The StructureData does
        not contain this type of possibilities.

    """

    PARSABLE_QUANTITIES = {}

    def __init__(self, *, handler=None, data=None, settings=None, options=None):  # pylint: disable=unused-argument
        super(BaseFileParser, self).__init__()
        # Make sure we only accept initialization with either `handler` or `data`.
        if (handler is not None and data is not None) or (handler is None and data is None):
            raise TypeError('Supply at bare minimum either argument handler or data to initialize parser.')

        # Make sure logger messages in the parser are passed to the AiiDA logger.
        self._logger = aiidalogger.getChild(self.__class__.__name__)

        # A few defaults
        self._exit_code = None
        # What quantities the specific content parser can provide.
        self._parsable_quantities = self.PARSABLE_QUANTITIES
        # The container for the parsed data when the `get_quantity` is executed, i.e. in the node composer
        # at a later stage.
        self._parsed_content = {}
        # The content parser, which will be an instance of one of the parsevasp parser classes.
        self._content_parser = None
        # Content data, which is an AiiDA data structure.
        self._content_data = None
        # Parser settings.
        self._settings = settings
        # Parser options.
        self._options = options

        # Set `handler`, `data` or raise if we supply something else in the parameters
        if handler is not None:
            self._init_from_handler(handler)
        if data is not None:
            self._init_from_data(data)

    @property
    def parsable_quantities(self):
        """Fetch the quantities that this content parser can provide."""
        return self._parsable_quantities

    @property
    def exit_code(self):
        """Fetch exit code."""
        return self._exit_code

    def get_quantity(self, quantity_key):
        """Fetch the required quantity from the content parser.

        Either fetch it from an existing AiiDA data structure, a parsed content dictionary
        if that exists, otherwise parse this specific quantity using the loaded instance,
        which is now a specific content parser.

        Parameters
        ----------
        quantity_key : str
            A string specifying the key of the quantity to be fetched.

        Returns
        -------
        result : object
            If we have initialized the content parser with an AiiDA data structure, we return it in
            its original form. If the `quantity_key` is not find to be parsable by this content
            parser None is returned. Finally, if the content parser was initialized with a file like
            object and the requested `quantity_key` is found to be parsable we return the quantity.

        """

        if self._content_data is not None:
            # If we have already set an AiiDA data structure, return it.
            # This is straightforward in our case as there is for the PoscarParser,
            # KpointsParser a 1:1 mapping between the parser and the AiiDA data (if
            # we ignore conversions between representations etc.).
            return self._content_data

        # We continue assuming we need to parse this quantity
        if quantity_key not in self._parsable_quantities:
            # Check if this quantity can be parsed by this content parser.
            return None

        if self._parsed_content.get(quantity_key) is None:
            # Parsed content does not contain this quantity,
            # most likely not parsed. Parse it now and store.
            self._parsed_content = self._parse_content()

        return self._parsed_content.get(quantity_key)

    def write(self, path):
        """
        Writes VASP content to file.

        Uses the write method defined in this loaded content parser.

        :param path: A string describing the relative path in the submission folder
        to write the file.

        """

        if self._content_parser is None or self._content_data is None:
            # Only write if we have an AiiDA data structure or parser prepared.
            if self._content_parser is None:
                # If we do not have a parser loaded before write, we have an
                # AiiDA data structure. Make sure that is on the form parsevasp expects.
                self._content_parser = self._content_data_to_content_parser()
            # Now a content parser should be ready and its content can be
            # written using parsevasp
            with open(path, 'w') as handler:
                self._content_parser.write(file_handler=handler)
        else:
            raise ValueError('The content parser has not been initialized or no ' 'AiiDA data structure is preparred.')

    def _init_from_handler(self, data):
        raise NotImplementedError('{classname} does not implement a _init_from_handler() '
                                  'method.'.format(classname=self.__class__.__name__))

    def _init_from_data(self, data):
        raise NotImplementedError('{classname} does not implement a _init_from_data() ' 'method.'.format(classname=self.__class__.__name__))

    def _parse_content(self, inputs):
        """Abstract base method to parse content.

        Has to be overwritten by the child class, which typically is a specific content parser.

        :param inputs:

        """

        raise NotImplementedError('{classname} does not implement a _parse_content() ' 'method.'.format(classname=self.__class__.__name__))
