"""
Test creation and preparation of an AmnCalculation
"""
import tempfile
import os

from aiida.backends.testbase import AiidaTestCase
from aiida.orm import CalculationFactory, Code
from aiida.common.folders import SandboxFolder

from .common import Common


class AmnCalcTest(AiidaTestCase):
    """Test creation and preparation of an AmnCalculation"""

    # pylint: disable=protected-access
    def setUp(self):
        """Set up test environment"""
        self.calc_cls = CalculationFactory('vasp.amn')
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

    def _get_calc(self, kpoints=None, no_wdat=False):
        """Create a calculation for testing"""
        if not kpoints:
            kpoints = Common.kpoints_mesh()
        calc = self.calc_cls()
        calc.use_code(self.code)
        calc.set_computer(self.computer)
        calc.use_settings(Common.settings())
        calc.inp.settings.update_dict({'icharg': 11})
        calc.use_structure(Common.cif())
        calc.use_kpoints(kpoints)
        calc.use_paw(Common.paw_in(), kind='In')
        calc.use_paw(Common.paw_as(), kind='As')
        calc.use_charge_density(Common.charge_density())
        calc.use_wavefunctions(Common.wavefunctions())
        calc.use_wannier_settings(Common.win())
        if not no_wdat:
            calc.use_wannier_data(self.wdat)
        return calc, calc.get_inputs_dict()

    def test_verify(self):
        """Check that correct inputs are verified successfully"""
        calc, inp = self._get_calc(kpoints=Common.kpoints_list())
        calc.verify_inputs(inp)
        calc.use_kpoints(Common.kpoints_mesh())
        inp = calc.get_inputs_dict()
        calc.verify_inputs(inp)
        calc.use_settings(Common.settings())
        inp = calc.get_inputs_dict()
        calc.verify_inputs(inp)

    # pylint: disable=protected-access
    def test_prepare(self):
        """Check that all input files are written"""
        calc, inp = self._get_calc()
        with SandboxFolder() as sandbox_f:
            calc_info = calc._prepare_for_submission(sandbox_f, inp)
            inputs_list = sandbox_f.get_content_list()
        self.assertEquals(
            set(inputs_list), {
                'INCAR', 'KPOINTS', 'POSCAR', 'POTCAR', 'CHGCAR', 'WAVECAR',
                'wannier90.win', 'test1', 'test2'
            })
        self.assertIn(['wannier90*', '.', 0], calc_info.retrieve_list)
        calc.use_settings(Common.settings())
        inp = calc.get_inputs_dict()
        calc.verify_inputs(inp)
        with SandboxFolder() as sandbox_f:
            calc._prepare_for_submission(sandbox_f, inp)
            inputs_list = sandbox_f.get_content_list()
        self.assertEquals(
            set(inputs_list), {
                'INCAR', 'KPOINTS', 'POSCAR', 'POTCAR', 'WAVECAR',
                'wannier90.win', 'test1', 'test2'
            })
        calc, inp = self._get_calc(no_wdat=True)
        with SandboxFolder() as sandbox_f:
            calc_info = calc._prepare_for_submission(sandbox_f, inp)
            inputs_list = sandbox_f.get_content_list()
        self.assertEquals(
            set(inputs_list), {
                'INCAR', 'KPOINTS', 'POSCAR', 'POTCAR', 'CHGCAR', 'WAVECAR',
                'wannier90.win'
            })

    def test_write_chgcar(self):
        """Check that CHGCAR file is written correctly"""
        calc, inp = self._get_calc()
        calc.write_chgcar(inp, self.tmpf)
        with open(self.tmpf, 'r') as chg:
            res = chg.read()
        self.assertEquals(res, Common.charge_density_res())

    def test_write_wavecar(self):
        """Check that WAVECAR file is written correctly"""
        calc, inp = self._get_calc()
        calc.write_wavecar(inp, self.tmpf)
        with open(self.tmpf, 'r') as wav:
            res = wav.read()
        self.assertEquals(res, Common.wavefunctions_res())

    def test_parse_with_retrieved(self):
        """Check that parsing is successful and that wannier_data is among the outputs"""
        calc, _ = self._get_calc()
        pars = calc.get_parserclass()(calc)
        success, outs = pars.parse_with_retrieved({
            'retrieved':
            Common.retrieved_nscf()
        })
        outs = dict(outs)
        self.assertTrue(success)
        self.assertIn('wannier_data', outs)
