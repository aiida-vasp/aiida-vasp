from aiida.djsite.db.testbase import AiidaTestCase
from aiida.orm import CalculationFactory, Code
from aiida.common.folders import SandboxFolder
from common import Common


class ScfCalcTest(AiidaTestCase):
    def setUp(self):
        self.calc_cls = CalculationFactory('vasp.scf')
        self.code = Code()
        self.code.set_computer(self.computer)
        self.code.set_remote_computer_exec((self.computer, '/bin/foo'))
        Common.import_paw()

    def tearDown(self):
        pass

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
        self.assertRaises(AttributeError, calc.verify_inputs, inp)
        calc.use_kpoints(Common.kpoints_mesh())
        inp = calc.get_inputs_dict()
        calc.verify_inputs(inp)

    def test_prepare(self):
        calc, inp = self._get_calc()
        with SandboxFolder() as sf:
            ci = calc._prepare_for_submission(sf, inp)
        self.assertIn('CHGCAR', ci.retrieve_list)
        self.assertIn('WAVECAR', ci.retrieve_list)
        self.assertIn('IBZKPT', ci.retrieve_list)
