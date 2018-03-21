"""Test that Poscar and MultiPotcar work together correctly."""
# pylint: disable=unused-import,redefined-outer-name,unused-argument,unused-wildcard-import,wildcard-import
import re

import pytest

from aiida_vasp.io.potcar import MultiPotcarIo
from aiida_vasp.utils.aiida_utils import get_data_class
from aiida_vasp.utils.fixtures import *
from aiida_vasp.utils.fixtures.data import POTCAR_MAP


@pytest.mark.parametrize(['vasp_structure'], [('str',)], indirect=True)
def test_poscar_potcar(fresh_aiida_env, vasp_structure_poscar, potcar_family):
    """Test that the order of potentials corresponds between POSCAR and POTCAR."""
    poscar_lines = vasp_structure_poscar.poscar_str().splitlines()
    poscar_counts = [int(count) for count in poscar_lines[5].split(' ')]
    index = 0
    poscar_sites = poscar_lines[7:]
    poscar_pot_order = []
    for count in poscar_counts:
        poscar_pot_order.append(poscar_sites[index].rsplit(' ', 1)[1])
        index += count
    potcar_cls = get_data_class('vasp.potcar')
    potcar_dict = potcar_cls.get_potcars_dict(elements=['As', 'In', 'In_d'], family_name=potcar_family, mapping=POTCAR_MAP)
    mulpotio = MultiPotcarIo.from_structure(structure=vasp_structure_poscar.structure, potentials_dict=potcar_dict)
    assert [re.sub('_.*', '', symbol) for symbol in vasp_structure_poscar.potentials_order] == poscar_pot_order
    assert [potio.node.full_name for potio in mulpotio.potcars] == [
        POTCAR_MAP[element] for element in vasp_structure_poscar.potentials_order
    ]
