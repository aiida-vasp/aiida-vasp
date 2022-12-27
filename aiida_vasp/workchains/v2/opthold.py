"""
Module containing the OptionHolder class
"""
from typing import List, Tuple

from aiida.common.exceptions import InputValidationError
from aiida.common.extendeddicts import AttributeDict
from aiida.orm import Dict

#pylint:disable=raise-missing-from


class Option(property):
    """
    Base class for descriptors to be used as pre-defined fields for `OptionContainer`.

    The point of using these seemingly complex descriptors is to allow the target `OptionContainer`
    class to support pre-defined properties acting as the fields to be set with the following
    functionalities:

    * Tab completion of the field name
    * Assignment time checks of the correct object types
    * Default values at a per `OptionHolder` subclass level
    * Enforcement of a field being required, e.g. no default value is available.
    * Automatic type conversion where necessary

    Note that the inheritance from 'property' is need for IPython introspect to work, but the usage
    is rather differently than the actual 'property' (as decorators/factories).
    However, the instantiated objects both sever as descriptors in a similar way.
    """

    def __init__(self, docstring, default_value=None, required=False):
        """Initialise an option and passing the docstring"""
        super().__init__()
        self.__doc__ = docstring
        self.required = required
        self.default_value = default_value
        self.name = None

    def __set_name__(self, owner, name):
        """Methods for automatically setting the `name` attribute - works for python 3.6+ only"""
        self.name = name

    def __get__(self, obj, owner=None):
        """Get the stored value"""
        if obj is None:
            return self
        if self.required and self.name not in obj._opt_data:
            raise ValueError(f'Field {self.name} has not been set yet!')

        return obj._opt_data.get(self.name, self.default_value)

    def __set__(self, obj, value):
        obj._opt_data[self.name] = value

    def __delete__(self, obj):
        """Delete the option from the holder dictionary"""
        if self.name in obj._opt_data:
            del obj._opt_data[self.name]


class TypedOption(Option):
    """Class for an option that enforces a specific type"""

    target_type = bool

    def __init__(self, docstring, default_value=None, required=False, enforce_type=False):
        """
        Instantiate an TypedOption field

        If ``enforce_type`` is True, will strictly check the type of the passed value.
        Otherwise, the value will be converted into the target type using the default constructor.
        """
        super().__init__(docstring, default_value, required)
        self.enforce_type = enforce_type

    def __set__(self, obj, value):
        """Setter for setting the option"""
        if self.enforce_type:
            if isinstance(value, self.target_type):
                obj._opt_data[self.name] = value
            else:
                raise ValueError(f'{value} is not a {self.target_type} type!')
        else:
            obj._opt_data[self.name] = self.target_type(value)

    def __get__(self, obj, owner=None):
        if obj is None:
            return self

        raw_value = super().__get__(obj, owner)
        if raw_value is not None:
            return self.target_type(raw_value)
        return None


class ChoiceOption(Option):
    """Option that only allow certain values"""

    def __init__(self, docstring, choices, default_value=None, required=False):
        super().__init__(docstring, default_value, required)
        self.choices = choices

    def __set__(self, obj, value):
        """Setter that sets the field"""
        if value not in self.choices:
            raise ValueError(f'{value} is not a valid choice, choose from: {self.choices}.')
        obj._opt_data[self.name] = value


class BoolOption(TypedOption):
    """Class for an option that accepts bool values"""

    target_type = bool


class FloatOption(TypedOption):
    """Class for an option that accepts float values"""

    target_type = float


class IntOption(TypedOption):
    """Class for an option that accepts integer values"""

    target_type = int


class DictOption(TypedOption):
    """Class for an option that accepts a dictionary"""

    target_type = dict


class ListOption(TypedOption):
    """Class for an option that accepts a list"""

    target_type = list


class ListOrStringOption(Option):
    """Class for an option that accepts a list"""


class StringOption(TypedOption):
    """Class for an option that accepts only string values"""

    def __init__(self, docstring, default_value=None, required=False, enforce_type=True):
        """Instantiate an object, note that we enforce_type by default here."""
        super().__init__(
            docstring,
            default_value=default_value,
            required=required,
            enforce_type=enforce_type,
        )


class OptionContainer:
    """
    Base class for a container of options
    """

    def __init__(self, **kwargs):
        """
        A holder of options

        Arguments:
            kwargs: unpack keyword arguments and set them as the attributes
        """
        self._opt_data = {}

        # the key word arguments are interpreted as the options to be set, an we check
        # If they are valid in the first place
        (
            self.valid_options,
            self.required_options,
        ) = self._get_valid_and_required_options()
        for key, value in kwargs.items():
            if key in self.valid_options:
                setattr(self, key, value)
            else:
                raise ValueError(f'{key} is not a valid option for a {type(self)} instance!')

    def _get_valid_and_required_options(self) -> Tuple[list, list]:
        """
        Return a list of valid option names

        This method extracts a tuple of list of the valid option names and those that are
        marked as `required`.
        """
        options = []
        required = []
        class_dict = vars(type(self))
        for name in dir(self):
            if name.startswith('__'):
                pass
            # Check if the name exists in the class's __dict__
            if name not in class_dict:
                continue
            optobj = class_dict.get(name)

            # Check if the name is an Option
            if isinstance(optobj, Option):
                options.append(name)
                if optobj.required:
                    required.append(name)
        return options, required

    @property
    def _invalid_attributes(self) -> List[str]:
        """Any attribute store inside __dict__ is an invalid option"""
        known = ['_opt_data', 'valid_options', 'required_options']
        invalid = []
        for key in self.__dict__:
            if key not in known:
                invalid.append(key)
        return invalid

    def to_dict(self, check_invalids=True) -> AttributeDict:
        """Return a python dict representation - including all fields"""
        # Check for any additional attributes
        invalid_attrs = self._invalid_attributes
        if check_invalids and invalid_attrs:
            raise ValueError(f'The following attributes are not valid options: {invalid_attrs}')

        outdict = {}
        for key in self.valid_options:
            value = getattr(self, key)
            # Note that 'None' is interpreted specially meaning that the option should not be present
            if value is not None:
                outdict[key] = value
        return AttributeDict(outdict)

    def to_aiida_dict(self):
        """Return an ``aiida.orm.Dict`` presentation"""

        python_dict = self.to_dict()
        return Dict(dict=python_dict)

    def __setitem__(self, key, value) -> None:
        """Set items - we call the setattr method"""
        if key not in self.valid_options:
            raise KeyError(f'{key} is not an valid option for a {type(self)} instance.')
        setattr(self, key, value)

    def to_string(self) -> str:
        """In string format"""
        return (self.to_dict(check_invalids=False).__repr__().replace('AttributeDict', ''))

    def __getitem__(self, key):
        """Set items - we just the getattr method"""
        return getattr(self, key)

    def __repr__(self):
        string = self.to_string()
        string = string.replace('\n', ' ')
        string = string.replace('#', '')
        string = string.strip()
        if len(string) > 60:
            string = string[:60] + '...'
        return f'{type(self).__name__}<{string}>'

    @classmethod
    def validate_dict(cls, input_dict, port=None) -> None:  #pylint:disable=unused-argument
        """
        Vaildate a dictionary/Dict node, this can be used as the validator for
        the Port accepting the inputs
        """
        obj = cls(**input_dict)
        all_options = list(obj.valid_options)  # A copy for all options

        # Are we receiving a Dict object?
        if not isinstance(input_dict, dict):
            input_dict = input_dict.get_dict()

        for key in input_dict.keys():
            if key not in all_options:
                raise InputValidationError(f"Key '{key}' is not a valid option")
            all_options.remove(key)

        # Check for any missing required fields
        missing = [key for key in all_options if key in obj.required_options]
        if missing:
            raise InputValidationError(f'There are missing options: {missing}')

        # Check for any problems with obj.to_dict()
        try:
            obj.to_dict()
        except ValueError as error:
            raise InputValidationError(f'Error during validation: {error.args}')

    @classmethod
    def serialise(cls, value):
        """
        Serialise a dictionary into Dict

        This method can be passed as a `serializer` key word parameter of for the `spec.input` call.
        """
        obj = cls(**value)
        return obj.to_aiida_dict()

    @classmethod
    def setup_spec(cls, spec, port_name, **kwargs) -> None:
        """Setup the spec for this input"""
        # Check if we have 'required' fields
        obj = cls()
        # The constructor is different with/without any required_options
        # If there are any required options, it does not make any sense to have a default for the port.
        if obj.required_options:
            spec.input(
                port_name,
                validator=cls.validate_dict,
                serializer=cls.serialise,
                **kwargs,
            )
        else:
            spec.input(
                port_name,
                validator=cls.validate_dict,
                default=lambda: cls().to_aiida_dict(),
                serializer=cls.serialise,
                **kwargs,
            )

    @classmethod
    def get_description(cls):
        """
        Return a string for the options of a OptionContains in a human-readable format.
        """

        obj = cls()
        template = '{:>{width_name}s}:  {:10s} \n{default:>{width_name2}}: {}'
        entries = []
        for name in obj.valid_options:
            if name not in obj.required_options:
                value = getattr(obj, name)
                # Each entry is name, type, doc, default value
                entries.append([name, getattr(cls, name).__doc__, str(type(value)), value])
            else:
                entries.append([name, getattr(cls, name).__doc__, 'Undefined', 'None (required)'])

        max_width_name = max(len(entry[0]) for entry in entries) + 2

        lines = []
        for entry in entries:
            lines.append(template.format(
                *entry,
                width_name=max_width_name,
                width_name2=max_width_name + 10,
                default='Default',
            ))
        return '\n'.join(lines)
