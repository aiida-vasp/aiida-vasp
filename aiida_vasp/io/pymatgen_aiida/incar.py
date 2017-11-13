"""Translate between pymatgen Incar class and Aiida datastructures."""
from pymatgen.io.vasp import Incar

from aiida_vasp.io.pymatgen_aiida.vasprun import get_data_node


class IncarToAiida(object):
    """Adapt between Incar object and AiiDA Nodes."""

    def __init__(self, incar_obj):
        self.incar_obj = incar_obj

    def get_aiida_node(self):
        return get_data_node('parameter', dict=self._normalize_for_aiida(self.incar_obj))

    def get_pymatgen(self):
        return self.incar_obj

    def write_file(self, file_path):
        return self.incar_obj.write_file(file_path)

    @classmethod
    def from_aiida_node(cls, param_node):
        data = param_node.get_dict()
        incar_dict = self._normalize_for_pymatgen(data)
        return cls(Incar.from_dict(incar_dict))

    @classmethod
    def from_dict(cls, input_dict):
        return cls(Incar.from_dict(cls._normalize_for_pymatgen(input_dict)))

    @classmethod
    def from_string(cls, input_string):
        return cls(Incar.from_string(input_string))

    @classmethod
    def from_file(cls, file_path):
        return cls(Incar.from_file(file_path))

    @classmethod
    def _normalize_for_aiida(cls, input_dict):
        return {k.lower(): Incar.proc_value(k, v) for k, v in input_dict.items()}

    @classmethod
    def _normalize_for_pymatgen(cls, input_dict):
        return {k.upper(): Incar.proc_value(k, v) for k, v in input_dict.items()}
