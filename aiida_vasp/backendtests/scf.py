"""Test creation and preparation of ScfCalculation"""
from aiida.backends.testbase import AiidaTestCase
from aiida.orm import CalculationFactory, Code
from aiida.common.folders import SandboxFolder

from .common import Common


class ScfCalcTest(AiidaTestCase):
    """Test creation and preparation of ScfCalculation"""

    def setUp(self):
        """Set up test environment"""
        self.calc_cls = CalculationFactory('vasp.scf')
        self.code = Code()
        self.code.set_computer(self.computer)
        self.code.set_remote_computer_exec((self.computer, '/bin/foo'))
        Common.import_paw()

    def tearDown(self):
        pass

    def _get_calc(self, kpoints=None):
        """Create a calculation to test"""
        if not kpoints:
            kpoints = Common.kpoints_mesh()
        calc = self.calc_cls()
        calc.use_code(self.code)
        calc.set_computer(self.computer)
        calc.use_parameters(Common.parameters())
        calc.use_structure(Common.cif())
        calc.use_kpoints(kpoints)
        calc.use_paw(Common.paw_in(), kind='In')
        calc.use_paw(Common.paw_as(), kind='As')
        return calc, calc.get_inputs_dict()

    def test_verify(self):
        """Check that correct inputs are successfully verified"""
        calc, inp = self._get_calc(kpoints=Common.kpoints_list())
        self.assertRaises(AttributeError, calc.verify_inputs, inp)
        calc.use_kpoints(Common.kpoints_mesh())
        inp = calc.get_inputs_dict()
        calc.verify_inputs(inp)

    # pylint: disable=protected-access
    def test_prepare(self):
        """Check that preparing writes all necessary files"""
        calc, inp = self._get_calc()
        with SandboxFolder() as sandbox_f:
            calc_info = calc._prepare_for_submission(sandbox_f, inp)
        self.assertIn('CHGCAR', calc_info.retrieve_list)
        self.assertIn('WAVECAR', calc_info.retrieve_list)
        self.assertIn('IBZKPT', calc_info.retrieve_list)
