"""Test creation and preparation of WannierCalculation"""
import tempfile
import os

from aiida.backends.testbase import AiidaTestCase
from aiida.orm import CalculationFactory, Code
from aiida.common.folders import SandboxFolder

from .common import Common


# pylint: disable=protected-access
class WannierCalcTest(AiidaTestCase):
    """Test creation and preparation of WannierCalculation"""

    def setUp(self):
        """Set up test environment"""
        self.calc_cls = CalculationFactory('vasp.wannier')
        self.code = Code()
        self.code.set_computer(self.computer)
        self.code.set_remote_computer_exec((self.computer, '/bin/foo'))
        Common.import_paw()
        self.tmpd, self.tmpf = tempfile.mkstemp()
        self.wdat = Common.wdat()
        self.wdat.add_file(self.tmpf, 'test1')
        self.wdat.add_file(self.tmpf, 'test2')
        self.wdat._make_archive()

    def tearDown(self):
        os.remove(self.tmpf)

    def _get_calc(self):
        """Create a calculation for testing"""
        calc = self.calc_cls()
        calc.use_code(self.code)
        calc.set_computer(self.computer)
        calc.use_parameters(Common.win())
        calc.inp.parameters.update_dict({'hr_plot': True})
        calc.use_data(self.wdat)
        return calc, calc.get_inputs_dict()

    def test_verify(self):
        calc, inp = self._get_calc()
        calc.verify_inputs(inp)

    def test_prepare(self):
        """Check that preparing creates all necessary files"""
        calc, inp = self._get_calc()
        with SandboxFolder() as sandbox_f:
            calc_info = calc._prepare_for_submission(sandbox_f, inp)
            inputs = sandbox_f.get_content_list()
        self.assertEquals(set(inputs), {'wannier90.win', 'test1', 'test2'})
        self.assertIn(['wannier90*', '.', 0], calc_info.retrieve_list)

    def test_write_win(self):
        """Check that the .win file is written correctly"""
        calc, inp = self._get_calc()
        calc.write_win(inp, self.tmpf)
        with open(self.tmpf, 'r') as win:
            res = win.read()
        win_res = Common.win_res() + '\nhr_plot = T'
        self.assertEquals(res, win_res)

    def test_parse_with_retrieved(self):
        """Check that parsing is successful and that output links are created correctly"""
        calc, _ = self._get_calc()
        pars = calc.get_parserclass()(calc)
        success, _ = pars.parse_with_retrieved({
            'retrieved':
            Common.retrieved_nscf()
        })
        self.assertTrue(success)
