"""Test common code in VaspCalcBase"""
from aiida.backends.testbase import AiidaTestCase
from aiida.orm import CalculationFactory, Code


class VaspCalcBaseTest(AiidaTestCase):
    """Test common code in VaspCalcBase"""

    def setUp(self):
        self.calc_cls = CalculationFactory('vasp.base.VaspCalcBase')

    def test_check_inputs_fail(self):
        calc = self.calc_cls()
        self.assertRaises(ValueError, calc.check_input, calc.get_inputs_dict(), 'code')
        self.assertRaises(ValueError, calc.verify_inputs, calc.get_inputs_dict())

    def test_check_inputs_succ(self):
        calc = self.calc_cls()
        calc.use_code(Code())
        self.assertTrue(calc.check_input(calc.get_inputs_dict(), 'code'))
        self.assertTrue(calc.verify_inputs(calc.get_inputs_dict()))
