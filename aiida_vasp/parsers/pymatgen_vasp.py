"""A pymatgen based VaspCalculation parser"""
from pymatgen.io.vasp.outputs import Vasprun

from aiida_vasp.parsers.base import BaseParser
from aiida_vasp.data.pymatgen.vasprun import VasprunToAiida


class PymatgenParser(BaseParser):
    """
    Use pymatgen to parse vasp output.

    * Translate pymatgen Vasprun class into output nodes
    * can read broken ``vasprun.xml`` files

    Parameters can be passed to the ``pymatgen.io.vasp.outputs.Vasprun``
    constructor via the settings input node for VaspCalculation (see example).

    Example Usage::

        calc = CalculationFactory('vasp.vasp')
        calc.set_parser_cls(ParserFactory('vasp.pymatgen'))
        calc.use_settings(DataFactory('parameter')(dict={
            'pymatgen_parser': {
                'exception_on_bad_xml': False,
                'parse_dos': False
            }
        }

        # ... set up calc further
        calc.submit()
    """

    def __init__(self, calc):
        super(PymatgenParser, self).__init__(calc)
        self.vasprun_adapter = None

    def parse_with_retrieved(self, retrieved):
        """Parse the calculation"""
        self.check_state()
        super_result = super(PymatgenParser,
                             self).parse_with_retrieved(retrieved)

        if not super_result[0]:
            return super_result

        vasprun_path = self.get_file('vasprun.xml')
        parsed_vasprun = Vasprun(vasprun_path, **self.get_vasprun_options())
        self.vasprun_adapter = VasprunToAiida(parsed_vasprun)

        self.add_node('structure', self.vasprun_adapter.last_structure)
        self.add_node('kpoints', self.vasprun_adapter.actual_kpoints)

    def get_vasprun_options(self):
        settings_input = self._calc.get_inputs_dic().get('settings')
        if not settings_input:
            return {}
        return settings_input.get_dict().get('pymatgen_parser', {})
