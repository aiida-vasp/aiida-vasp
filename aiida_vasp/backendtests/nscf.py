"""Test creating and preparing NscfCalculation"""
import tempfile
import os

from aiida.backends.testbase import AiidaTestCase
from aiida.orm import CalculationFactory, Code
from aiida.common.folders import SandboxFolder

from .common import Common


class NscfCalcTest(AiidaTestCase):
    """Test creating and preparing NscfCalculation"""

    def setUp(self):
        """Set up test environment"""
        self.calc_cls = CalculationFactory('vasp.nscf')
        self.code = Code()
        self.code.set_computer(self.computer)
        self.code.set_remote_computer_exec((self.computer, '/bin/foo'))
        Common.import_paw()
        self.tmpd, self.tmpf = tempfile.mkstemp()

    def tearDown(self):
        os.remove(self.tmpf)

    def _get_calc(self, kpoints=None):
        """Create a calculation for testing"""
        if not kpoints:
            kpoints = Common.kpoints_mesh()
        calc = self.calc_cls()
        calc.use_code(self.code)
        calc.set_computer(self.computer)
        calc.use_parameters(Common.parameters())
        calc.inp.parameters.update_dict({'icharg': 11})
        calc.use_structure(Common.cif())
        calc.use_kpoints(kpoints)
        calc.use_paw(Common.paw_in(), kind='In')
        calc.use_paw(Common.paw_as(), kind='As')
        calc.use_charge_density(Common.charge_density())
        calc.use_wavefunctions(Common.wavefunctions())
        return calc, calc.get_inputs_dict()

    def test_verify(self):
        """Check correct inputs get verified successfully"""
        calc, inp = self._get_calc(kpoints=Common.kpoints_list())
        calc.verify_inputs(inp)
        calc.use_kpoints(Common.kpoints_mesh())
        inp = calc.get_inputs_dict()
        calc.verify_inputs(inp)
        calc.use_parameters(Common.parameters())
        inp = calc.get_inputs_dict()
        calc.verify_inputs(inp)

    # pylint: disable=protected-access
    def test_prepare(self):
        """Check that preparing creates all necessary files"""
        calc, inp = self._get_calc()
        with SandboxFolder() as sandbox_f:
            calc_info = calc._prepare_for_submission(sandbox_f, inp)
            inputs = sandbox_f.get_content_list()
        self.assertEquals(
            set(inputs),
            {'INCAR', 'KPOINTS', 'POSCAR', 'POTCAR', 'CHGCAR', 'WAVECAR'})
        self.assertIn('EIGENVAL', calc_info.retrieve_list)
        self.assertIn('DOSCAR', calc_info.retrieve_list)
        self.assertIn(['wannier90*', '.', 0], calc_info.retrieve_list)
        calc.use_parameters(Common.parameters())
        inp = calc.get_inputs_dict()
        calc.verify_inputs(inp)
        with SandboxFolder() as sandbox_f:
            calc._prepare_for_submission(sandbox_f, inp)
            inputs = sandbox_f.get_content_list()
        self.assertEquals(
            set(inputs), {'INCAR', 'KPOINTS', 'POSCAR', 'POTCAR', 'WAVECAR'})

    def test_write_chgcar(self):
        """Test that CHGAR file is written correctly"""
        calc, inp = self._get_calc()
        calc.write_chgcar(inp, self.tmpf)
        with open(self.tmpf, 'r') as chg:
            res = chg.read()
        self.assertEquals(res, Common.charge_density_res())

    def test_write_wavecar(self):
        """Test that WAVECAR file is written correctly"""
        calc, inp = self._get_calc()
        calc.write_wavecar(inp, self.tmpf)
        with open(self.tmpf, 'r') as wav:
            res = wav.read()
        self.assertEquals(res, Common.wavefunctions_res())

    def test_parse_with_retrieved(self):
        """Check that parsing is successful and creates the right output links"""
        calc, _ = self._get_calc()
        pars = calc.get_parserclass()(calc)
        success, outs = pars.parse_with_retrieved({
            'retrieved':
            Common.retrieved_nscf()
        })
        outs = dict(outs)
        self.assertTrue(success)
        self.assertIn('bands', outs)
        self.assertIn('dos', outs)
        self.assertIn('wannier_parameters', outs)
        self.assertIn('wannier_data', outs)
        self.assertIn('results', outs)
