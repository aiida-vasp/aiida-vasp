from vasp5 import Vasp5Calculation
import unittest
import tempfile


class Vasp5CalcTest(unittest.TestCase):
    '''
    Test Case for py:class:`~aiida.orm.calculation.job.vasp.vasp5.Vasp5Calculation`
    '''
    def setUp(self):
        self.calc = Vasp5Calculation()

    def test_inputs(self):
        self.assertTrue(hasattr(self.calc, 'use_code'))
        self.assertTrue(hasattr(self.calc, 'use_settings'))
        self.assertTrue(hasattr(self.calc, 'use_structure'))
        self.assertTrue(hasattr(self.calc, 'use_paw'))
        self.assertTrue(hasattr(self.calc, 'use_kpoints'))
        self.assertTrue(hasattr(self.calc, 'use_charge_density'))
        self.assertTrue(hasattr(self.calc, 'use_wavefunctions'))

    def test_internal_params(self):
        self.assertTrue(self.caclc.get_parser_name())

    def test_settings_property(self):
        self.calc.use_settings(self.calc.new_settings(dict={'A': 0}))
        self.assertEqual(self.calc.settings, {'a': 0})

    def test_write_incar(self):
        '''
        write out an INCAR tag to a tempfile and check wether
        it was written correctly
        '''
        inc = self.calc.new_settings(dict={'sytsem': 'InAs'})
        dst = tempfile.mkstemp()[1]
        self.calc.use_settings(inc)
        self.calc.write_incar({}, dst)
        with open(dst, 'r') as incar:
            self.assertEqual(incar.read().strip(), 'SYSTEM = InAs')

    def test_write_paw(self):
        '''
        concatenate two paws into a tmp POTCAR and check wether
        each is contained in the result
        '''
        self.calc.use_paw(self.Paw.load_paw(family='LDA', 'In_d'), kind='In')
        self.calc.use_paw(self.Paw.load_paw(family='LDA', 'As'), kind='As')
        dst = tempfile.mkstemp()[1]
        self.calc.write_paw({}, dst)
        with open(dst, 'r') as potcar:
            x = potcar.read()
            with open(self.calc.inp.paw_In, 'r') as paw_In:
                a = paw_In.read()
                self.assertIn(a, x)
            with open(self.calc.inp.paw_As, 'r') as paw_As:
                b = paw_As.read()
                self.assertIn(b, x)

    def test_write_structure(self, inputdict, dst):
        '''
        feed a structure into calc and write it to a POSCAR temp file
        check for nonemptiness of the file
        '''
        larray = array([[0,.5,.5],
                        [.5,0,.5],
                        [.5,.5,0]])
        alat = 6.058
        structure = StructureData(cell=larray*alat)
        structure.append_atom(position=[0, 0, 0], symbols='In')
        structure.append_atom(position=[.25,.25,.25], symbols='As')
        self.calc.use_structure(structure)
        dst = tempfile.mkstemp()[1]
        self.calc.write_structure(self, {}, dst)
        with open(dst, 'r') as poscar:
            self.assertTrue(poscar.read())

    def test_need_kp_false(self):
        self.calc.use_settings(
            self.calc.new_settings(dict={'kspacing': 0.5, 'kgamma': True}))
        self.assertFalse(self.calc._need_kp())

    def test_need_kp_true(self):
        self.calc.use_settings(self.calc.new_settings())
        self.assertTrue(self.calc._need_kp())

    def test_need_chgd_none(self):
        self.calc.use_settings(self.calc.new_settings())
        self.assertFalse(self.calc._need_chgd())

    def test_need_chgd_icharg(self):
        for i in [0, 2, 4, 10, 12]:
            self.calc.use_settings(
                self.calc.new_settings(dict={'icharg': i}))
            self.assertFalse(self.calc._need_chgd())
        for i in [1, 11]:
            self.calc.use_settings(
                self.calc.new_settings(dict={'icharg': i}))
            self.assertTrue(self.calc._need_chgd())

    def test_need_wfn_none(self):
        self.assertFalse(self.calc.need_wfn())
        self.calc.use_wavefunctions(self.calc.new_wavefunctions())
        self.assertTrue(self.calc.need_wfn())

    def test_need_wfn_istart(self):
        self.calc.use_settings(
            self.calc.new_settings(dict={'istart': 0}))
        self.assertFalse(self.calc.need_wfn())
        for i in [1, 2, 3]:
            self.calc.use_settings(
                self.calc.new_settings(dict={'istart': i}))
            self.assertTrue(self.calc.need_wfn())

    def test_get_paw_linkname(self):
        self.assertEqual(self.calc._get_paw_linkname('In'), 'paw_In')
