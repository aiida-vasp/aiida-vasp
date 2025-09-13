from aiida_vasp.utils.temp_profile import *
print(load_temp_profile())


# Uncomment the below line to create a localhost Computer if you have not done so
comp = orm.Computer('localhost', 'localhost', transport_type='core.local', scheduler_type='core.direct')
comp.store()


# Some configuration may be needed for first-time user
from pathlib import Path
import os
from subprocess import check_output
comp.set_workdir('/tmp/aiida_run/')
comp.configure()
vasp_path = check_output(['which', 'mock-vasp'],universal_newlines=True).strip()
vasp_code = orm.InstalledCode(comp, vasp_path, default_calc_job_plugin='vasp.vasp')
vasp_code.label ='mock-vasp'
vasp_code.store()
os.environ['MOCK_VASP_REG_BASE'] = str((Path() / 'mock_registry').absolute())
os.environ['MOCK_VASP_UPLOAD_PREFIX'] = 'relax'

# Upload the POTCAR files
from aiida_vasp.data.potcar import PotcarData, PotcarFileData
from pathlib import Path

print(PotcarData.upload_potcar_family(str(Path('potcars').absolute()), "PBE.EXAMPLE", "PBE.EXAMPLE"))

# Setting up the silicon structure
from ase.build import bulk
si = bulk('Si', 'diamond', 5.4)
si_node = orm.StructureData(ase=si)

from aiida import orm
from aiida_vasp.protocols.generator import VaspRelaxInputGenerator

upd = VaspRelaxInputGenerator()
# Override the vasp.potential_family input of the builder
upd.get_builder(structure=si_node, code='mock-vasp@localhost', overrides={
   'vasp': {'potential_family': 'PBE.EXAMPLE'}}
   )
upd.builder

results = upd.run_get_node()
