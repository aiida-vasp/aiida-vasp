"""
Module containing the OptionHolder class
"""

from aiida.common.exceptions import InputValidationError
from aiida.orm import Dict
from pydantic import BaseModel, ValidationError

# pylint:disable=raise-missing-from


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


class OptionContainer(BaseModel):
    """
    Base class for a container of options
    """

    def aiida_dict(self):
        """Return an ``aiida.orm.Dict`` presentation"""

        python_dict = self.model_dump()
        return Dict(dict=python_dict)

    @classmethod
    def validate_dict(cls, input_dict) -> None:  # pylint:disable=unused-argument
        """
        Vaildate a dictionary/Dict node, this can be used as the validator for
        the Port accepting the inputs
        """
        try:
            cls(**input_dict)
        except ValidationError as error:
            raise InputValidationError(f'Input validation failed: {error}') from error
        return True

    @classmethod
    def serialize(cls, python_dict: dict):
        """
        serialize a dictionary into Dict

        This method can be passed as a `serializer` key word parameter of for the `spec.input` call.
        """
        obj = cls(**python_dict)
        return obj.aiida_dict()

    @classmethod
    def get_description(cls):
        """
        Return a string for the options of a OptionContains in a human-readable format.
        """

        obj = cls()
        template = '{:>{width_name}s}:  {:10s} \n{default:>{width_name2}}: {}'
        entries = []
        for name, field in obj.model_fields.items():
            # Each entry is name, type, doc, default value
            entries.append([name, str(field.annotation.__name__), field.description, field.default])
        max_width_name = max(len(entry[0]) for entry in entries) + 2

        lines = []
        for entry in entries:
            lines.append(
                template.format(
                    *entry,
                    width_name=max_width_name,
                    width_name2=max_width_name + 10,
                    default='Default',
                )
            )
        return '\n'.join(lines)
