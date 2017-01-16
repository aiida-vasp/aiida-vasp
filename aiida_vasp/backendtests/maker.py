from aiida.djsite.db.testbase import AiidaTestCase
from aiida.orm.calculation.job.vasp.maker import VaspMaker


class VaspMakerTest(AiidaTestCase):
    '''
    py:class:VaspMakerTest:
        tests for VaspMaker
    '''
    def setUp(self):
        self.tm = VaspMaker()
