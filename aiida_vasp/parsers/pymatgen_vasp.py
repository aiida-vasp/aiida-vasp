"""A pymatgen based VaspCalculation parser"""
import xml
from pymatgen.io.vasp.outputs import Vasprun
from aiida.common.exceptions import ParsingError as AiidaParsingError

from aiida_vasp.parsers.base import BaseParser
from aiida_vasp.data.pymatgen.vasprun import VasprunToAiida


class PymatgenParser(BaseParser):
    """
    Use pymatgen to parse vasp output.

    * Translate ``pymatgen`` ``Vasprun`` class into output nodes
    * can read broken ``vasprun.xml`` files

    Parameters can be passed to the ``pymatgen.io.vasp.outputs.Vasprun``
    constructor via the settings input node for VaspCalculation (see example).

    :note: on broken ``vasprun.xml``

        VASP needs to have completed at least one ionic step and must have written
        the corresponding <calculation>...</calculation> tag to the vasprun.xml output
        or the parser will still fail. This is a limitation of the ``pymatgen``'s ``Vasprun``
        class.

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

    _linkname_outparams = 'output_parameters'

    def __init__(self, calc):
        super(PymatgenParser, self).__init__(calc)
        self.vasprun_adapter = None

    def parse_with_retrieved(self, retrieved):
        """Parse the calculation"""
        self.check_state()
        success, _ = super(PymatgenParser,
                           self).parse_with_retrieved(retrieved)

        if not success:
            return self.result(success)

        parsed_vasprun = self.try_parse_vasprun()
        self.vasprun_adapter = VasprunToAiida(parsed_vasprun)

        self.add_node('structure', self.vasprun_adapter.last_structure)
        self.add_node('kpoints', self.vasprun_adapter.actual_kpoints)
        self.add_node('forces', self.vasprun_adapter.forces)
        self.add_node('output_parameters',
                      self.vasprun_adapter.output_parameters)

        return self.result(success)

    def try_parse_vasprun(self):
        """
        Try to parse the vasprun file

        and give helpful error messages in case something goes wrong
        """
        vasprun_path = self.get_file('vasprun.xml')
        try:
            parsed_vasprun = Vasprun(vasprun_path,
                                     **self.get_vasprun_options())
        except (xml.etree.ElementTree.ParseError,
                xml.etree.cElementTree.ParseError):
            msg_tpl = '{problem}, {solution}'
            problem = 'vasprun.xml contains invalid xml'
            solution_tpl = 'try parsing it anyway using\n{example}'
            example = '\n'.join([
                '', 'vasp_calc.use_settings(ParameterData(dict={',
                '    "pymatgen_parser": {',
                '        "exception_on_bad_xml": False', '    }', '}))', ''
            ])
            msg = msg_tpl.format(
                problem=problem, solution=solution_tpl.format(example=example))
            raise AiidaParsingError(msg)
        except IndexError:
            raise AiidaParsingError(
                'vasprun.xml must contain at least one complete ionic step (<calculation> tag)'
            )
        return parsed_vasprun

    def get_vasprun_options(self):
        settings_input = self._calc.get_inputs_dict().get('settings')
        if not settings_input:
            return {}
        return settings_input.get_dict().get('pymatgen_parser', {})
