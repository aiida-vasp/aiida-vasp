from aiida.orm import CalculationFactory, Code, DataFactory
from aiida.djsite.db.testbase import AiidaTestCase
from aiida.tools.codespecific.vasp.io.incar import IncarParser
from common import Common
import tempfile


class VaspCalcBaseTest(AiidaTestCase):
    def setUp(self):
        self.calc_cls = CalculationFactory('vasp.base.BasicCalculation')
        Common.import_paw()
        Paw = DataFactory('vasp.paw')
        self.code = Code()
        self.code.set_computer(self.computer)
        self.code.set_remote_computer_exec((self.computer, '/bin/foo'))
        self.paw_in = Paw.load_paw(element='In')[0]
        self.paw_as = Paw.load_paw(element='As')[0]
        self.calc = self._get_calc('c', 'm')
        self.inp = self.calc.get_inputs_dict()
        self.tmp, self.tmpf = tempfile.mkstemp()

    def tearDown(self):
        import os
        os.remove(self.tmpf)

    def _get_calc(self, stype, ktype):
        calc = self.calc_cls()
        calc.use_code(self.code)
        calc.set_computer(self.computer)
        calc.set_resources({'num_machines': 1, 'num_mpiprocs_per_machine':1})
        calc.use_settings(Common.settings())
        if stype == 's':
            calc.use_structure(Common.structure())
        else:
            calc.use_structure(Common.cif())
        if ktype == 'm':
            calc.use_kpoints(Common.kpoints_mesh())
        else:
            calc.use_kpoints(Common.kpoints_list())
        calc.use_paw(self.paw_in, kind='In')
        calc.use_paw(self.paw_as, kind='As')
        return calc

    def test_store(self):
        c_sm= self._get_calc('s', 'm')
        c_sm.store_all()
        self.assertIsNotNone(c_sm.pk)

        c_sl = self._get_calc('s', 'l')
        c_sl.store_all()
        self.assertIsNotNone(c_sl.pk)

        c_cm = self._get_calc('c', 'm')
        c_cm.store_all()
        self.assertIsNotNone(c_cm.pk)

        c_cl = self._get_calc('c', 'l')
        c_cl.store_all()
        self.assertIsNotNone(c_cm.pk)

    def test_write_incar(self):
        self.calc.write_incar(self.inp, self.tmpf)
        icp = IncarParser(self.tmpf)
        for k, v in self.calc.inp.settings.get_dict().iteritems():
            self.assertIn(str(v), icp.result[k])
