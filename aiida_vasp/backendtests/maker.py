from aiida.djsite.db.testbase import AiidaTestCase
from aiida_vasp.calcs.maker import VaspMaker


class VaspMakerTest(AiidaTestCase):
    '''
    py:class:VaspMakerTest:
        tests for VaspMaker
    '''

    def setUp(self):
        self.tm = VaspMaker()
