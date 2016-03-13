from aiida.orm import CalculationFactory
from aiida.djsite.db.testbase import AiidaTestCase
# ~ from os.path import realpath, join, dirname


class VaspCalcBaseTest(AiidaTestCase):
    def setUp(self):
        self.calc_cls = CalculationFactory('vasp.base.BasicCalculation')
