from aiida import load_dbenv, is_dbenv_loaded
if not is_dbenv_loaded():
    load_dbenv()

from aiida.orm import Code, Computer
from aiida.orm.querytool import QueryTool
from aiida_vasp.calcs.maker import VaspMaker
import sys


cifname = sys.argv[1]
mkcalc = VaspMaker(structure=cifname)
mkcalc.code = Code.get_from_string('asevasp@monch')
mkcalc.kpoints.set_kpoints_mesh([8, 8, 8])
mkcalc.add_parameters(
    system=mkcalc.structure.filename,
    npar=8
)
mkcalc.recipe = 'test_sc'

mkcalc.computer = Computer.get('monch')
mkcalc.queue = 'dphys_compute'

v5 = mkcalc.new()
v5.set_resources({
    'num_machines': 8,
    'num_mpiprocs_per_machine': 2})
# ~ v5.set_max_memory_kb(8000000)

tag = sys.argv[2]
q = QueryTool()
q.set_class(mkcalc.calc_cls)
q.add_extra_filter('tag', '=', tag)
ql = map(lambda c: c.get_extra('test-nr'), q.run_query())

last_tn = ql and max(ql) or 0

v5.store_all()

v5.set_extras({'tag': tag, 'test-nr': last_tn+1})
print 'pk: {}, test-nr: {}'.format(v5.pk, last_tn+1)
v5.submit()
