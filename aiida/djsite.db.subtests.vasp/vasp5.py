from aiida.orm import CalculationFactory, DataFactory
from aiida.djsite.db.testbase import AiidaTestCase
import numpy as np
import tempfile
from os.path import dirname, realpath, join


class Vasp5CalcTest(AiidaTestCase):
    '''
    Test Case for
    py:class:`~aiida.orm.calculation.job.vasp.vasp5.Vasp5Calculation`
    '''
    def setUp(self):
        self.calc = CalculationFactory('vasp.vasp5')()
        DataFactory('vasp.paw').import_family(
            realpath(join(dirname(__file__), 'LDA')), family_override='TEST')

        larray = np.array([[0, .5, .5],
                           [.5, 0, .5],
                           [.5, .5, 0]])
        alat = 6.058
        self.structure = DataFactory('structure')(cell=larray*alat)
        self.structure.append_atom(position=[0, 0, 0], symbols='In')
        self.structure.append_atom(position=[.25, .25, .25], symbols='As')

        cifpath = realpath(join(dirname(__file__),
                                'data', 'EntryWithCollCode43360.cif'))
        self.cif = DataFactory('cif').get_or_create(cifpath)[0]

    def test_inputs(self):
        self.assertTrue(hasattr(self.calc, 'use_code'))
        self.assertTrue(hasattr(self.calc, 'use_settings'))
        self.assertTrue(hasattr(self.calc, 'use_structure'))
        self.assertTrue(hasattr(self.calc, 'use_paw'))
        self.assertTrue(hasattr(self.calc, 'use_kpoints'))
        self.assertTrue(hasattr(self.calc, 'use_charge_density'))
        self.assertTrue(hasattr(self.calc, 'use_wavefunctions'))

    def test_internal_params(self):
        self.assertTrue(self.calc.get_parser_name())

    def test_settings_property(self):
        self.calc.use_settings(self.calc.new_settings(dict={'A': 0}))
        self.assertEqual(self.calc._settings, {'a': 0})

    def test_write_incar(self):
        '''
        write out an INCAR tag to a tempfile and check wether
        it was written correctly
        '''
        inc = self.calc.new_settings(dict={'system': 'InAs'})
        dst = tempfile.mkstemp()[1]
        self.calc.use_settings(inc)
        self.calc.write_incar({}, dst)
        with open(dst, 'r') as incar:
            self.assertEqual(incar.read().strip(), 'SYSTEM = InAs')

    def test_write_potcar(self):
        '''
        concatenate two paws into a tmp POTCAR and check wether
        each is contained in the result
        '''
        self.calc.use_structure(self.structure)
        self.calc.use_paw(
            self.calc.load_paw(family='TEST', symbol='In_d'), kind='In')
        self.calc.use_paw(
            self.calc.load_paw(family='TEST', symbol='As'), kind='As')
        dst = tempfile.mkstemp()[1]
        self.calc.write_potcar(self.calc.get_inputs_dict(), dst)
        with open(dst, 'r') as potcar:
            x = potcar.read()
            with open(self.calc.inp.paw_In.potcar, 'r') as paw_In:
                a = paw_In.read()
                self.assertIn(a, x)
            with open(self.calc.inp.paw_As.potcar, 'r') as paw_As:
                b = paw_As.read()
                self.assertIn(b, x)

    def test_write_poscar(self):
        '''
        feed a structure into calc and write it to a POSCAR temp file
        check for nonemptiness of the file
        '''
        self.calc.use_structure(self.structure)
        dst = tempfile.mkstemp()[1]
        self.calc.write_poscar({}, dst)
        with open(dst, 'r') as poscar:
            self.assertTrue(poscar.read())

    def test_write_poscar_cif(self):
        '''
        feed a cif file into calc and write it to a POSCAR temp file
        make sure the file is not empty
        '''
        self.calc.use_structure(self.cif)
        dst = tempfile.mkstemp()[1]
        self.calc.write_poscar({}, dst)
        with open(dst, 'r') as poscar:
            self.assertTrue(poscar.read())

    def test_write_kpoints(self):
        '''
        feed kpoints into calc and write to KPOINTS temp file
        verify the file is not empty
        '''
        kp = self.calc.new_kpoints()
        kp.set_kpoints_mesh([4, 4, 4])
        self.calc.use_kpoints(kp)
        self.calc.use_settings(self.calc.new_settings())
        dst = tempfile.mkstemp()[1]
        self.calc.write_kpoints(self.calc.get_inputs_dict(), dst)
        with open(dst, 'r') as kpoints:
            self.assertTrue(kpoints.read())

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
        self.calc.use_settings(self.calc.new_settings())
        self.assertFalse(self.calc._need_wfn())
        self.calc.use_wavefunctions(self.calc.new_wavefunctions())
        self.assertTrue(self.calc._need_wfn())

    def test_need_wfn_istart(self):
        self.calc.use_settings(
            self.calc.new_settings(dict={'istart': 0}))
        self.assertFalse(self.calc._need_wfn())
        for i in [1, 2, 3]:
            self.calc.use_settings(
                self.calc.new_settings(dict={'istart': i}))
            self.assertTrue(self.calc._need_wfn(),
                            msg='_need_wfn not True for istart=%s' % i)

    def test_get_paw_linkname(self):
        self.assertEqual(self.calc._get_paw_linkname('In'), 'paw_In')

    def test_Paw(self):
        self.assertIs(self.calc.Paw, DataFactory('vasp.paw'))

    def test_load_paw(self):
        paw_A = self.calc.load_paw(family='TEST', symbol='As')
        paw_B = self.calc.Paw.load_paw(family='TEST', symbol='As')[0]
        self.assertEqual(paw_A.pk, paw_B.pk)

    def test_new_setting(self):
        self.assertIsInstance(self.calc.new_settings(),
                              DataFactory('parameter'))

    def test_new_structure(self):
        self.assertIsInstance(self.calc.new_structure(),
                              DataFactory('structure'))

    def test_new_kpoints(self):
        self.assertIsInstance(self.calc.new_kpoints(),
                              DataFactory('array.kpoints'))

    def test_new_charge_density(self):
        self.assertIsInstance(self.calc.new_charge_density(),
                              DataFactory('vasp.chargedensity'))

    def test_new_wavefunctions(self):
        self.assertIsInstance(self.calc.new_wavefunctions(),
                              DataFactory('vasp.wavefun'))
