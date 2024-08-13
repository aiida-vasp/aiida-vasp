"""
Module containing the OptionHolder class
"""

from aiida.orm import Dict
from pydantic import BaseModel, ValidationError

# pylint:disable=raise-missing-from


class OptionContainer(BaseModel):
    """
    Base class for a container of options
    """

    def aiida_dict(self):
        """Return an ``aiida.orm.Dict`` presentation"""

        python_dict = self.model_dump()
        return Dict(dict=python_dict)

    @classmethod
    def aiida_validate(cls, input_dict, namespace=None) -> None:  # pylint:disable=unused-argument
        """
        Validate a dictionary/Dict node, this can be used as the validator for
        the Port accepting the inputs
        """
        if isinstance(input_dict, Dict):
            input_dict = input_dict.get_dict()
        try:
            cls(**input_dict)
        except ValidationError as error:
            return str(error)
        return None

    @classmethod
    def aiida_serialize(cls, python_dict: dict):
        """
        serialize a dictionary into Dict

        This method can be passed as a `serializer` key word parameter of for the `spec.input` call.
        """
        obj = cls(**python_dict)
        return obj.aiida_dict()

    @classmethod
    def aiida_description(cls):
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
