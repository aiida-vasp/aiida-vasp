"""
Parser module for composing output of a VASP calculation.

The simplified parser outputs the following nodes:

1. A `misc` node that stores simply summary information such as total energies,
    total run times, any warnings issues during the calculation and if the calculation
    was finished.
2. A `arrays` node that for storing the properties that are arrays by nature
    such as the dielectric function, hessian matrix as well as the energies of each
    SCF cycle.
3. A `trajectory` node for storing the trajectory of geometry optimisation and AIMD.
4. A `band` node for storing the band structure.
5. A `dos` node for storing the band structure.
6. A `chargcar` node for parsed electronic density as array.

Main difference from the previous version
1. Individual properties are not stored in the separate node in order to reduce the number of nodes created. Scalar and
    small arrays are stored in misc (final energies and forces), and larger none-standard output arrays are stored in
    the `array` output node.
2. Standard outputs each as the bands and dos are stored in dedicated nodes.
3. `pydantic` is used to validate the parser settings *at submission* time.
"""

from typing import List

import numpy as np
from aiida import orm
from aiida.parsers.parser import Parser
from pydantic import BaseModel, Field

from aiida_vasp.parsers.content_parsers import *


class ParserError(RuntimeError):
    pass


class QuantityMissingError(ParserError):
    pass


class RequiredQuantityMissingError(ParserError):
    pass


DEFAULT_QUANTITIES = (
    'total_energies',
    'maximum_stress',
    'maximum_force',
    'notifications',
    'run_status',
    'run_stats',
    'version',
    'band_properties',
)

DEFAULT_EXCLUDED_QUANTITIES = ('energies', 'trajectory', 'kpoints', 'chgcar', 'projectors')

DEFAULT_EXCLUDED_NODE = tuple()

DEFAULT_REQUIRED_QUANTITIES = ('run_status', 'run_stats', 'total_energies')


class ParserSettingsConfig(BaseModel):
    """
    Settings for the VASP parser.
    """

    include_quantity: List[str] = Field(description='Properties to include', default_factory=lambda: [])
    exclude_quantity: List[str] = Field(description='Quantities to be excluded', default_factory=lambda: [])
    required_quantity: List[str] = Field(
        description='Quantities that most be present', default_factory=lambda: list(DEFAULT_REQUIRED_QUANTITIES)
    )
    include_node: List[str] = Field(description='Output node to include', default_factory=lambda: [])
    exclude_node: List[str] = Field(description='Output node to exclude', default_factory=lambda: [])
    kpoints_from_ibzkpt: bool = False


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

        user_config: ParserSettingsConfig = ParserSettingsConfig(
            **self.node.inputs.settings.get_dict().get('parser_settings', {})
        )
        self.user_config = user_config

        quantities_each = {}

        # Apply the modifiers
        quantities_to_exclude = [key for key in DEFAULT_EXCLUDED_QUANTITIES if key not in user_config.include_quantity]
        quantities_to_exclude += user_config.exclude_quantity
        nodes_to_exclude = [key for key in DEFAULT_EXCLUDED_NODE if key not in user_config.include_node]
        nodes_to_exclude += user_config.exclude_node

        # Parse the files
        with self.retrieved.open('vasprun.xml', 'rb') as handler:
            parser = VasprunParser(handler=handler)
            quantities_each['vasprun.xml'] = parser.get_all_quantities()

        with self.retrieved.open('OUTCAR', 'r') as handler:
            parser = OutcarParser(handler=handler)
            quantities_each['outcar'] = parser.get_all_quantities()

        with self.retrieved.open('vasp_output', 'r') as handler:
            parser = StreamParser(handler=handler)
            quantities_each['vasp_output'] = parser.get_all_quantities()

        with self.retrieved.open('CONTCAR', 'r') as handler:
            parser = PoscarParser(handler=handler)
            quantities_each['contcar'] = parser.get_all_quantities()

        if 'chgcar' not in quantities_to_exclude:
            with self.retrieved.open('CHGCAR', 'r') as handler:
                parser = ChgcarParser(handler=handler)
                quantities_each['chgcar'] = parser.get_all_quantities()

        if user_config.kpoints_from_ibzkpt:
            with self.retrieved.open('IBZKPT', 'r') as handler:
                parser = KpointsParser(handler=handler)
                quantities_each['ibzkpt'] = parser.get_all_quantities()

        # Do we really need DOSCAR and EIGENVAL
        # if user_config.parse_doscar is True:
        #     with self.retrieved.open('DOSCAR', 'r'):
        #         parser = DoscarParser(handler=handler, settings=user_config)
        #         quantities_each['doscar'] = parser.get_all_quantities()

        # if user_config.parse_eigenval is True:
        #     with self.retrieved.open('EIGENVAL', 'r'):
        #         parser = EigenvalParser(handler=handler, settings=user_config)
        #         quantities_each['eigenval'] = parser.get_all_quantities()

        # Remove the quantities
        for name, parsed_quantities in quantities_each.items():
            for sub_key in list(parsed_quantities.keys()):
                if sub_key in quantities_to_exclude:
                    del parsed_quantities[sub_key]

        # Check in required quantities are present
        for name in user_config.required_quantity:
            exists = False
            for _, value in quantities_each.items():
                if value.get(name) is not None:
                    exists = True
                    break
            if exists is False:
                raise RequiredQuantityMissingError(f'Required quantity {name} is missing.')

        # Create the outputs
        for name in ['misc', 'structure', 'trajectory', 'kpoints', 'arrays', 'band', 'dos']:
            if name in nodes_to_exclude:
                continue
            node = None
            try:
                node = getattr(self, '_compose_' + name)(quantities_each)
            except QuantityMissingError:
                continue
            if node is not None:
                self.out(name, node)

    def _compose_misc(self, quantities_each):
        misc_quantities = (
            'total_energies',
            'maximum_stress',
            'maximum_force',
            'notifications',
            'run_status',
            'run_stats',
            'version',
            'final_forces',
            'final_stress' 'site_magnetization' 'band_properties',
        )

        out_dict = {}
        gather_quantities(quantities_each, 'vasprun.xml', out_dict, misc_quantities)
        gather_quantities(quantities_each, 'OUTCAR', out_dict, misc_quantities)
        return orm.Dict(dict=out_dict)

    def _compose_structure(self, quantities_each):
        """Parse structure from the vasprun.xml file, if failed, use
        that from CONTCAR
        """

        data = None
        if 'vasprun.xml' in quantities_each:
            data = quantities_each['vasprun.xml'].get('structure')
        if data is None:
            data = quantities_each['CONTCAR']['structure']
        if data is None:
            raise QuantityMissingError()

        node = orm.StructureData()
        node.set_cell(data['unitcell'])
        for site in data['sites']:
            node.append_atom(position=site['position'], symbols=site['symbol'], name=site['kind_name'])
        return node

    def _compose_arrays(self, quantities_each):
        """Generate the generic `Arrays` output"""
        array_quantities = ('projectors', 'born_charges', 'dielectrics', 'hessian', 'dynmat', 'energies')
        out_arrays = {}
        gather_quantities(quantities_each, 'vasprun.xml', out_arrays, array_quantities, flatten_dict=True)
        # Remove None values in the arrays
        out_arrays = {key: value for key, value in out_arrays.items() if value is not None}
        if out_arrays:
            return orm.ArrayData(out_arrays)
        return None

    def _compose_kpoints(self, quantities_each):
        kpoints_data = None
        if self.user_config.kpoints_from_ibzkpt is True:
            kpoints_data = quantities_each['ibzkpt']['kpoints']
        elif 'vasprun.xml' in quantities_each:
            kpoints_data = quantities_each['vasprun.xml'].get('kpoints')

        if kpoints_data is not None:
            node = orm.KpointsData()
            if kpoints_data['mode'] == 'explicit':
                node.set_kpoints(
                    kpoints_data['points'], weights=kpoints_data['weights'], cartesian=kpoints_data['cartesian']
                )
            elif kpoints_data['mode'] == 'automatic':
                node.set_kpoints_mesh(kpoints_data['divisions'], offset=kpoints_data['shifts'])
            else:
                raise ValueError(f'Unknown kpoints mode {kpoints_data["mode"]}')
            return node
        raise QuantityMissingError('No valid kpoints data to use')

    def _compose_trajectory(self, quantities_each):
        """Compose the band output"""

        node = orm.TrajectoryData()

        if 'vasprun.xml' in quantities_each:
            traj_data = quantities_each['vasprun.xml'].get('trajectory')
            if traj_data is None:
                return None
            for item in traj_data:
                for key, value in traj_data[item].items():
                    if key == 'symbols':
                        node.base.attributes.set(key, value)
                    else:
                        node.set_array(key, value)
            return node
        return None

    def _compose_band(self, quantities_each):
        """Compose the band structure node"""
        if 'vasprun.xml' in quantities_each:
            deigen = quantities_each['vasprun.xml']['eigenvalues']
            docc = quantities_each['vasprun.xml']['occupancies']
            if 'total' in deigen:
                eigenvalues = np.array(deigen['total'])
                occupancies = np.array([docc['total']])
            else:
                eigenvalues = np.array([deigen['up'], deigen['down']])
                occupancies = np.array([docc['up'], docc['down']])
            node = orm.BandsData()
            kpoints = self._compose_kpoints(quantities_each)
            node.set_kpointsdata(kpoints)
            node.set_bands(eigenvalues, occupations=occupancies)
            return node

    def _compose_dos(self, quantities_each):
        """Compose the DOS node"""
        arrays_dict = {}
        if 'vasprun.xml' in quantities_each:
            gather_quantities(quantities_each, 'dos', arrays_dict, ['dos'], flatten_dict=True)
        if arrays_dict:
            node = orm.ArrayData(arrays_dict)
            return node


def gather_quantities(quantities_each, namespace, dst, fields, flatten_dict=False):
    """
    Gather quantities and put them into the target dictionary
    """
    for key, value in quantities_each.get(namespace, {}).items():
        if key in fields:
            if isinstance(value, dict) and flatten_dict:
                # flatten the dictionary - prepend the key with the name of the quantity
                for key2, value2 in value.items():
                    dst[key + '_' + key2] = value2
            else:
                dst[key] = value
