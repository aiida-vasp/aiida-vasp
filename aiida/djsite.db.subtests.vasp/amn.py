from aiida.djsite.db.testbase import AiidaTestCase
from aiida.orm import CalculationFactory, Code
from aiida.common.folders import SandboxFolder
from common import Common
import tempfile
import os


class AmnCalcTest(AiidaTestCase):
    def setUp(self):
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

    def _get_calc(self, kp=None):
        if not kp:
            kp = Common.kpoints_mesh()
        calc = self.calc_cls()
        calc.use_code(self.code)
        calc.set_computer(self.computer)
        calc.use_settings(Common.settings())
        calc.inp.settings.update_dict({'icharg': 11})
        calc.use_structure(Common.cif())
        calc.use_kpoints(kp)
        calc.use_paw(Common.paw_in(), kind='In')
        calc.use_paw(Common.paw_as(), kind='As')
        calc.use_charge_density(Common.charge_density())
        calc.use_wavefunctions(Common.wavefunctions())
        calc.use_wannier_settings(Common.win())
        calc.use_wannier_data(self.wdat)
        return calc, calc.get_inputs_dict()

    def test_verify(self):
        calc, inp = self._get_calc(kp=Common.kpoints_list())
        calc.verify_inputs(inp)
        calc.use_kpoints(Common.kpoints_mesh())
        inp = calc.get_inputs_dict()
        calc.verify_inputs(inp)
        calc.use_settings(Common.settings())
        inp = calc.get_inputs_dict()
        calc.verify_inputs(inp)

    def test_prepare(self):
        calc, inp = self._get_calc()
        with SandboxFolder() as sf:
            ci = calc._prepare_for_submission(sf, inp)
            il = sf.get_content_list()
        self.assertEquals(set(il),
                          {'INCAR', 'KPOINTS', 'POSCAR',
                           'POTCAR', 'CHGCAR', 'WAVECAR', 'wannier90.win',
                           'test1', 'test2'})
        self.assertIn(['wannier90*', '.', 0], ci.retrieve_list)
        calc.use_settings(Common.settings())
        inp = calc.get_inputs_dict()
        calc.verify_inputs(inp)
        with SandboxFolder() as sf:
            calc._prepare_for_submission(sf, inp)
            il = sf.get_content_list()
        self.assertEquals(set(il),
                          {'INCAR', 'KPOINTS', 'POSCAR',
                           'POTCAR', 'WAVECAR', 'wannier90.win',
                           'test1', 'test2'})

    def test_write_chgcar(self):
        calc, inp = self._get_calc()
        calc.write_chgcar(inp, self.tmpf)
        with open(self.tmpf, 'r') as chg:
            res = chg.read()
        self.assertEquals(res, Common.charge_density_res())

    def test_write_wavecar(self):
        calc, inp = self._get_calc()
        calc.write_wavecar(inp, self.tmpf)
        with open(self.tmpf, 'r') as wav:
            res = wav.read()
        self.assertEquals(res, Common.wavefunctions_res())

    def test_parse_with_retrieved(self):
        calc, inpt = self._get_calc()
        pars = calc.get_parserclass()(calc)
        ok, outs = pars.parse_with_retrieved({
            'retrieved': Common.retrieved_nscf()
        })
        outs = dict(outs)
        self.assertIn('wannier_data', outs)
