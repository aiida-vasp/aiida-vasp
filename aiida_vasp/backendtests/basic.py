from aiida.orm import CalculationFactory, Code, DataFactory
from aiida.djsite.db.testbase import AiidaTestCase
from aiida_vasp.utils.io.incar import IncarParser
from common import Common
import tempfile
import os


class BasicCalcTest(AiidaTestCase):
    def setUp(self):
        self.calc_cls = CalculationFactory('vasp.base.BasicCalculation')
        Common.import_paw()
        Paw = DataFactory('vasp.paw')
        self.code = Code()
        self.code.set_computer(self.computer)
        self.code.set_remote_computer_exec((self.computer, '/bin/foo'))
        self.paw_in = Paw.load_paw(element='In')[0]
        self.paw_as = Paw.load_paw(element='As')[0]
        self.tmp, self.tmpf = tempfile.mkstemp()

    def tearDown(self):
        import os
        os.remove(self.tmpf)

    def _get_calc(self, stype, ktype):
        calc = self.calc_cls()
        calc.use_code(self.code)
        calc.set_computer(self.computer)
        calc.set_resources({'num_machines': 1, 'num_mpiprocs_per_machine': 1})
        calc.use_parameters(Common.parameters())
        if stype == 's':
            calc.use_structure(Common.structure())
        else:
            calc.use_structure(Common.cif())
        if ktype == 'm':
            calc.use_kpoints(Common.kpoints_mesh())
        else:
            calc.use_kpoints(Common.kpoints_list())
        calc.use_paw(self.paw_in, kind='In')
        calc.use_paw(self.paw_as, kind='As')
        return calc

    def test_store(self):
        c_sm = self._get_calc('s', 'm')
        c_sm.store_all()
        self.assertIsNotNone(c_sm.pk)

        c_sl = self._get_calc('s', 'l')
        c_sl.store_all()
        self.assertIsNotNone(c_sl.pk)

        c_cm = self._get_calc('c', 'm')
        c_cm.store_all()
        self.assertIsNotNone(c_cm.pk)

        c_cl = self._get_calc('c', 'l')
        c_cl.store_all()
        self.assertIsNotNone(c_cm.pk)

    def test_write_incar(self):
        calc = self._get_calc('s', 'm')
        inp = calc.get_inputs_dict()
        calc.write_incar(inp, self.tmpf)
        icp = IncarParser(self.tmpf)
        for k, v in calc.inp.parameters.get_dict().iteritems():
            self.assertIn(str(v), icp.result[k])

    def test_write_poscar_structure(self):
        calc = self._get_calc('s', 'm')
        inp = calc.get_inputs_dict()
        from ase.io.vasp import read_vasp
        calc.write_poscar(inp, self.tmpf)
        wd = os.getcwd()
        os.chdir(os.path.dirname(self.tmpf))
        poscar = None
        poscar = read_vasp(self.tmpf)
        os.chdir(os.path.dirname(wd))
        self.assertIsNotNone(poscar)

    def test_write_poscar_cif(self):
        calc = self._get_calc('c', 'm')
        inp = calc.get_inputs_dict()
        from ase.io.vasp import read_vasp
        calc.write_poscar(inp, self.tmpf)
        wd = os.getcwd()
        os.chdir(os.path.dirname(self.tmpf))
        poscar = None
        poscar = read_vasp(self.tmpf)
        os.chdir(os.path.dirname(wd))
        self.assertIsNotNone(poscar)

    def test_write_kpoints_mesh(self):
        from aiida_vasp.utils.io.kpoints import KpParser
        calc = self._get_calc('c', 'm')
        inp = calc.get_inputs_dict()
        calc.write_kpoints(inp, self.tmpf)
        kpp = KpParser(self.tmpf)
        self.assertTrue((kpp.kpoints == Common.kpoints_mesh_res()).all())

    def test_write_kpoints_list(self):
        from aiida_vasp.utils.io.kpoints import KpParser
        calc = self._get_calc('c', 'l')
        inp = calc.get_inputs_dict()
        calc.write_kpoints(inp, self.tmpf)
        kpp = KpParser(self.tmpf)
        kres, wres = Common.kpoints_list_res()
        self.assertTrue((kpp.kpoints == kres).all())
        self.assertTrue((kpp.weights == wres).all())

    def test_write_potcar(self):
        calc = self._get_calc('c', 'm')
        inp = calc.get_inputs_dict()
        calc.write_potcar(inp, self.tmpf)
        with open(self.tmpf, 'r') as pcf:
            pcs = pcf.read()
        with open(os.path.expanduser('~/tmp/potcar'), 'w') as tpc:
            tpc.write(pcs)
        self.assertIn('In_d', pcs)
        self.assertIn('As', pcs)
        self.assertEquals(pcs.count('End of Dataset'), 2)

    def test_elements(self):
        calc = self._get_calc('c', 'm')
        self.assertRaises(AttributeError, calc.get_attr, 'elements')
        calc._prestore()
        self.assertEquals(['In', 'As'], calc.elements)

    def test_new_parameters(self):
        calc = self.calc_cls()
        calc.use_parameters(calc.new_parameters(dict={'bla': 3}))
        self.assertEquals(calc.inp.parameters.get_dict().get('bla'), 3)

    def test_new_structure(self):
        calc = self.calc_cls()
        calc.use_structure(calc.new_structure())

    def test_new_kpoints(self):
        calc = self.calc_cls()
        calc.use_kpoints(calc.new_kpoints())

    def test_load_paw(self):
        calc = self.calc_cls()
        calc.use_paw(calc.load_paw(
            symbol=self.paw_in.symbol,
            family='TEST'), kind='test')
