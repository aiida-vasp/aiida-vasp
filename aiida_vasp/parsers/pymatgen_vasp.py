"""Provide a pymatgen based VaspCalculation parser."""
import xml

from aiida.common.exceptions import ParsingError as AiidaParsingError
from pymatgen.io.vasp.outputs import Vasprun, Outcar

from aiida_vasp.io.pymatgen_aiida.vasprun import VasprunToAiida, get_data_node
from aiida_vasp.io.pymatgen_aiida.outcar import OutcarToAiida
from aiida_vasp.parsers.base import BaseParser


class PymatgenParser(BaseParser):
    """
    Use pymatgen_aiida to parse vasp output.

    * Translate ``pymatgen_aiida`` ``Vasprun`` class into output nodes
    * can read broken ``vasprun.xml`` files

    Parameters can be passed to the ``pymatgen_aiida.io.vasp.outputs.Vasprun``
    constructor via the settings input node for VaspCalculation (see example).

    generic parser settings (not passed on to pymatgen_aiida):

        * ``parse_dos``: [True] if false, suppress the dos output node
        * ``parse_bands``: [True] if false, suppress the bands output node

    :note: on broken ``vasprun.xml``

        VASP needs to have completed at least one ionic step and must have written
        the corresponding <calculation>...</calculation> tag to the vasprun.xml output
        or the parser will still fail. This is a limitation of the ``pymatgen_aiida``'s ``Vasprun``
        class.

    Example Usage::

        calc = CalculationFactory('vasp.vasp')
        calc.set_parser_cls(ParserFactory('vasp.pymatgen_aiida'))
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
    _linkname_kpoints = 'kpoints'
    _linkname_structure = 'structure'
    _linkname_forces = 'forces'
    _linkname_bands = 'bands'
    _linkname_dos = 'dos'
    _linkname_born_charges = 'born_charges'
    _default_options = {'parse_dos': True, 'parse_bands': True}

    def __init__(self, calc):
        super(PymatgenParser, self).__init__(calc)
        self.vasprun_adapter = None
        self.outcar_adapter = None

    def parse_with_retrieved(self, retrieved):
        """Parse the calculation."""
        self.check_state()
        success, _ = super(PymatgenParser, self).parse_with_retrieved(retrieved)

        if not success:
            return self.result(success)

        parsed_vasprun = self.try_parse_vasprun()
        self.vasprun_adapter = VasprunToAiida(parsed_vasprun, logger=self.logger)

        self.add_node(self._linkname_structure, self.vasprun_adapter.last_structure)
        self.add_node(self._linkname_kpoints, self.vasprun_adapter.actual_kpoints)
        self.add_node(self._linkname_forces, self.vasprun_adapter.forces)
        if self.get_option('parse_bands'):
            self.add_node(self._linkname_bands, self.vasprun_adapter.band_structure)
        if self.get_option('parse_dos'):
            self.add_node(self._linkname_dos, self.vasprun_adapter.tdos)

        parsed_outcar = self.parse_outcar()
        if parsed_outcar:
            self.outcar_adapter = OutcarToAiida(parsed_outcar)

            if parsed_outcar.lepsilon:
                self.add_node(self._linkname_born_charges, self.outcar_adapter.born_charges)

        self.add_node(self.get_linkname_outparams(), self.get_output_parameters())

        return self.result(success)

    def try_parse_vasprun(self):
        """
        Try to parse the vasprun file.

        and give helpful error messages in case something goes wrong
        """
        vasprun_path = self.get_file('vasprun.xml')
        try:
            parsed_vasprun = Vasprun(vasprun_path, **self.get_vasprun_options())
        except (xml.etree.ElementTree.ParseError, xml.etree.cElementTree.ParseError):
            msg_tpl = '{problem}, {solution}'
            problem = 'vasprun.xml contains invalid xml'
            solution_tpl = 'try parsing it anyway using\n{example}'
            example = '\n'.join([
                '', 'vasp_calc.use_settings(ParameterData(dict={', '    "pymatgen_parser": {', '        "exception_on_bad_xml": False',
                '    }', '}))', ''
            ])
            msg = msg_tpl.format(problem=problem, solution=solution_tpl.format(example=example))
            raise AiidaParsingError(msg)
        except IndexError:
            raise AiidaParsingError('vasprun.xml must contain at least one complete ionic step (<calculation> tag)')
        return parsed_vasprun

    def parse_outcar(self):
        """Try to parse the OUTCAR file."""
        outcar_path = self.get_file('OUTCAR')
        if not outcar_path:
            return None
        parsed_outcar = Outcar(outcar_path)
        return parsed_outcar

    def get_vasprun_options(self):
        """Get the options to the pymatgen_aiida Vasprun constructor from settings."""
        settings_input = self._calc.get_inputs_dict().get('settings')
        if not settings_input:
            return {}
        return settings_input.get_dict().get('pymatgen_parser', {})

    def get_options(self):
        """Get the generic parser options."""
        settings_input = self._calc.get_inputs_dict().get('settings')
        if not settings_input:
            return self._default_options
        return settings_input.get_dict().get('parser', self._default_options)

    def get_option(self, option):
        """Get a specific parser option."""
        return self.get_options()[option]

    def get_output_parameters(self):
        """Get the output parameter node."""
        output_parameter_dict = self.vasprun_adapter.output_dict
        if self.outcar_adapter:
            output_parameter_dict.update(self.outcar_adapter.output_dict)
        output_parameters = get_data_node('parameter', dict=output_parameter_dict)
        return output_parameters
