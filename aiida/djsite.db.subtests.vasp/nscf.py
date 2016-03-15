from aiida.djsite.db.testbase import AiidaTestCase
from aiida.orm import CalculationFactory, Code
from aiida.common.folders import SandboxFolder
from common import Common
import tempfile
import os


class NscfCalcTest(AiidaTestCase):
    def setUp(self):
        self.calc_cls = CalculationFactory('vasp.nscf')
        self.code = Code()
        self.code.set_computer(self.computer)
        self.code.set_remote_computer_exec((self.computer, '/bin/foo'))
        Common.import_paw()
        self.tmpd, self.tmpf = tempfile.mkstemp()

    def tearDown(self):
        os.remove(self.tmpf)

    def _get_calc(self, kp=None):
        if not kp:
            kp = Common.kpoints_mesh()
        calc = self.calc_cls()
        calc.use_code(self.code)
        calc.set_computer(self.computer)
        calc.use_settings(Common.settings())
        calc.use_structure(Common.cif())
        calc.use_kpoints(kp)
        calc.use_paw(Common.paw_in(), kind='In')
        calc.use_paw(Common.paw_as(), kind='As')
        return calc, calc.get_inputs_dict()

    def test_verify(self):
        calc, inp = self._get_calc(kp=Common.kpoints_list())
        calc.verify_inputs(inp)
        calc.use_kpoints(Common.kpoints_mesh())
        inp = calc.get_inputs_dict()
        calc.verify_inputs(inp)

    def test_prepare(self):
        calc, inp = self._get_calc()
        with SandboxFolder() as sf:
            ci = calc._prepare_for_submission(sf, inp)
        self.assertIn('EIGENVAL', ci.retrieve_list)
        self.assertIn('DOSCAR', ci.retrieve_list)
        self.assertIn('wannier90.win', ci.retrieve_list)
        self.assertIn('wannier90.mmn', ci.retrieve_list)
        self.assertIn('wannier90.amn', ci.retrieve_list)
        self.assertIn('wannier90.eig', ci.retrieve_list)
