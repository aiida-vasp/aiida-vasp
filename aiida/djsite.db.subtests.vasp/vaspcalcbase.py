from aiida.orm import CalculationFactory, Code
from aiida.djsite.db.testbase import AiidaTestCase


class VaspCalcBaseTest(AiidaTestCase):
    def setUp(self):
        self.calc_cls = CalculationFactory('vasp.base.VaspCalcBase')

    def test_check_inputs_Fail(self):
        calc = self.calc_cls()
        self.assertRaises(
            ValueError,
            calc.check_input,
            calc.get_inputs_dict(),
            'code'
        )
        self.assertRaises(
            ValueError,
            calc.verify_inputs,
            calc.get_inputs_dict()
        )

    def test_check_inputs_Succ(self):
        calc = self.calc_cls()
        calc.use_code(Code())
        self.assertTrue(
            calc.check_input(calc.get_inputs_dict(), 'code')
        )
        self.assertTrue(calc.verify_inputs(calc.get_inputs_dict()))
