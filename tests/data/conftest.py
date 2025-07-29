"""
Fixtures related to data.

-------------------------
Here different pytest fixtures are set up. They typically contain computers,
VASP inputs etc. which you would need to mock a VASP job in this plugin. It
also contains the set up of the mock VASP executable, which is used to test the
workchains.
"""

import pytest
from aiida import orm
from aiida.plugins import DataFactory

from aiida_vasp.data.potcar import PotcarData, PotcarFileData


@pytest.fixture
def vasp2w90_params(aiida_profile, vasp_params):
    vasp_params_data = vasp_params()
    incar_data = orm.Dict(dict=vasp_params_data.code.get_dict().update({'lwannier90': True}))
    return incar_data


@pytest.fixture
def potcar_node_pair(aiida_profile, data_path):
    """Create a POTCAR node pair."""
    potcar_path = data_path('potcar', 'As', 'POTCAR')
    potcar_file_node = PotcarFileData(file=potcar_path)
    potcar_file_node.store()
    return {'file': potcar_file_node, 'potcar': PotcarData.find_one(symbol='As')}


@pytest.fixture()
def ref_retrieved(data_path):
    """Fixture: retrieved directory from an NSCF vasp run."""
    retrieved = DataFactory('core.folder')()
    retrieved.put_object_from_tree(path=data_path('basic_run'))
    return retrieved
