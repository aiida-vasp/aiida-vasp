"""
Parser module for composing output of a VASP calculation.
"""

from dataclasses import dataclass
from typing import List, Type

from aiida import orm
from aiida.parsers.parser import Parser
from pydantic import BaseModel, Field

from aiida_vasp.parsers.content_parsers import *


class ParserSettingsConfig(BaseModel):
    """
    Settings for the VASP parser.
    """

    outputs: List[str] = Field(default_factory=lambda x: list(DEFAULT_OUTPUTS))
    required_outputs: List[str] = Field(default_factory=lambda x: list(DEFAULT_REQUIRED_OUTPUTS))
    critical_notifications: List[str] = Field(default_factory=lambda x: list(DEFAULT_CRITICAL_NOTIFICATIONS))
    parse_ibzkpt: bool = False
    parse_eigenval: bool = False
    parse_doscar: bool = False


class QuantityConfig(BaseModel):
    """Output node configuration"""

    link_name: str
    tyle: str
    quantities: List[str]


@dataclass
class OutputDestination:
    link_name: str
    entrypoint: Type


QUANTITIES_CONFIG = {
    'total_energies': OutputDestination('misc', orm.Dict),
    'maximum_stress': OutputDestination('misc', orm.Dict),
    'notifications': OutputDestination('misc', orm.Dict),
    'run_status': OutputDestination('misc', orm.Dict),
    'run_stats': OutputDestination('misc', orm.Dict),
    'version': OutputDestination('misc', orm.Dict),
    'kpoints': OutputDestination('kpoints', orm.KpointsData),
    'structure': OutputDestination('structure', orm.StructureData),
    'trajectory': OutputDestination('trajectory', orm.TrajectoryData),
    'forces': OutputDestination('arrays', orm.ArrayData),
    'stress': OutputDestination('stress', orm.ArrayData),
    'bands': OutputDestination('bands', orm.BandsData),
    'dos': OutputDestination('dos', orm.ArrayData),
    'energies': OutputDestination('arrays', orm.ArrayData),
    'projectors': OutputDestination('arrays', orm.ArrayData),
    'born_charges': OutputDestination('arrays', orm.ArrayData),
    'dielectrics': OutputDestination('arrays', orm.ArrayData),
}


class VaspParser(Parser):
    """Class for parsing VASP output files and storing the results in AiiDA."""

    def __init__(self, node):
        """
        Initialize the Parser instance
        """
        super(VaspParser, self).__init__(node)

    def parse(self, **kwargs):
        """
        Parse outputs, store results in database.
        """
        error_code = self.compose_output()
        if error_code is not None:
            return error_code
        return None

    def compose_output(self):
        """
        Compose the output of a VASP calculation.
        """

        settings_config: ParserSettingsConfig = ParserSettingsConfig.model_rebuild(
            self.node.inputs.settings.get_dict().get('parser_settings', {})
        )

        quantities = {}
        # Parse the files
        with self.retrieved.open('vasprun.xml', 'rb') as handler:
            parser = VasprunParser(handler=handler, settings=settings_config)
            quantities['vasprun'] = parser.get_all_quantities()

        with self.retrieved.open('OUTCAR', 'r') as handler:
            parser = OutcarParser(handler=handler, settings=settings_config)
            quantities['outcar'] = parser.get_all_quantities()

        with self.retrieved.open('vasp_output', 'r'):
            parser = PoscarParser(handler=handler, settings=settings_config)
            quantities['vasp_output'] = parser.get_all_quantities()

        with self.retrieved.open('CONTCAR', 'r'):
            parser = PoscarParser(handler=handler, settings=settings_config)
            quantities['contcar'] = parser.get_all_quantities()

        if 'chgcar' in settings_config:
            with self.retrieved.open('CHGCAR', 'r'):
                parser = ChgcarParser(handler=handler, settings=settings_config)
                quantities['chgcar'] = parser.get_all_quantities()

        if settings_config.parse_ibzkpt is True:
            with self.retrieved.open('IBZKPT', 'r'):
                parser = KpointsParser(handler=handler, settings=settings_config)
                quantities['ibzkpt'] = parser.get_all_quantities()

        if settings_config.parse_doscar is True:
            with self.retrieved.open('DOSCAR', 'r'):
                parser = DoscarParser(handler=handler, settings=settings_config)
                quantities['doscar'] = parser.get_all_quantities()

        if settings_config.parse_eigenval is True:
            with self.retrieved.open('EIGENVAL', 'r'):
                parser = EigenvalParser(handler=handler, settings=settings_config)
                quantities['eigenval'] = parser.get_all_quantities()

        # Compose the results
        out_node = {name: {} for name in set([config.link_name for config in QUANTITIES_CONFIG.values()])}

        for quantity_name, quantity_config in QUANTITIES_CONFIG.items():
            for fname, parsed_quantities in quantities.items():
                for name in parsed_quantities.keys():
                    if name == quantity_name:
                        if name in out_node[quantity_config.link_name]:
                            self.logger.info(f'Duplicate quantity {name}. Overwritten it with that from {fname}.')
                        out_node[quantity_config.link_name][quantity_name] = parsed_quantities[name]

        # Compose the misc node

        # Compose the structure node

        # Compose the arrays node

        # Compose the kpoints node

        # Compose the chgcar node
