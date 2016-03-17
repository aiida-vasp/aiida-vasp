from aiida.orm import CalculationFactory, DataFactory
from aiida.tools.codespecific.vasp.default_paws import lda, gw
import os


class VaspMaker(object):
    '''
    py:class:: VaspMaker(*args, **kwargs)
        simplifies the task of creating a Vasp5Calculation from scratch
        further simplifies creating certain often used types of calculations

        :keyword structure: A StructureData node or a
        (relative) path to either a .cif file or a POSCAR file. Defaults to
        :type structure: str or StructureData
        a new empty structure node recieved from calc_cls.

        :keyword calc_cls: the class that VaspMaker will use when creating Calculation nodes.
        defaults to 'vasp.vasp5'.
        if a string is given, it will be passed to aiida's CalculationFactory
        :type calc_cls: str or vasp.BasicCalculation subclass

        :keyword continue_from: A vasp calculation node with charge_density and
        wavefunction output links. VaspMaker will create calculations that
        start with those as inputs.
        :type continue_from: vasp calculation node

        :keyword copy_from: A vasp calculation. It's inputs will be used as defaults
        for the created calculations.
        :type copy_from: vasp calculation node

        :keyword charge_density: chargedensity node from a previously run calculation
        :type charge_density: ChargedensityData
        :keyword wavefunctions: wavefunctions node from a previously run calculation
        :type wavefunctions: WavefunData
        :keyword array.KpointsData kpoints: kpoints node to use for input
        :keyword str paw_family: The name of a PAW family stored in the db
        :keyword str paw_defaults: A dictionary mapping element symbols -> PAW symbols
        :keyword str label: value for the calculation label
        :keyword str computer: computer name, defaults to code's if code is given
        :keyword str code: code name, if any Calculations are given, defaults to their code
        :keyword str resources: defaults to copy_from.get_resources() or None
        :keyword str queue: defaults to queue from given calculation, if any, or None

        py:method:: new()
        :returns: an instance of :py:attr:calc_cls, initialized with the data held by the VaspMaker

        py:method:: add_settings(**kwargs)
        Adds keys to the settings (INCAR keywords), if settings is already stored, makes a copy.
        Does not overwrite previously set keywords.

        py:method:: rewrite_settings(**kwargs)
        Same as :py:meth:add_settings, but also overwrites keywords.

        py:method:: set_kpoints_mesh

        py:method:: set_kpoints_path

        py:attribute:: structure
        Used to initialize the created calculations as well as other nodes (like kpoints).
        When changed, can trigger changes in other data nodes.

        py:attribute:: calc_cls
        Vasp calculation class to be used in :py:meth:new

        py:attribute:: computer

        py:attribute:: code

        py:attribute:: queue

        py:attribute:: settings
        A readonly shortcut to the contents of the settings node

        py:attribute:: kpoints
        The kpoints node to be used, may be copied to have py:func:set_cell called.

        py:attribute:: wavefunction

        py:attribute:: charge_density

        py:attribute:: elements
        Chemical symbols of the elements contained in py:attr:structure
    '''
    def __init__(self, *args, **kwargs):
        self._init_defaults(*args, **kwargs)
        self._calcname = kwargs.get('calc_cls')
        if 'continue_from' in kwargs:
            self._init_from(kwargs['continue_from'])
        if 'copy_from' in kwargs:
            self._copy_from(kwargs['copy_from'])

    def _init_defaults(self, *args, **kwargs):
        calcname = kwargs.get('calc_cls', 'vasp.vasp5')
        if isinstance(calcname, (str, unicode)):
            self.calc_cls = CalculationFactory(calcname)
        else:
            self.calc_cls = calcname
        self.label = kwargs.get('label', 'unlabeled')
        self._computer = kwargs.get('computer')
        self._code = kwargs.get('code')
        self._settings = kwargs.get('settings', self.calc_cls.new_settings())
        self._set_default_structure(kwargs.get('structure'))
        self._paw_fam = kwargs.get('paw_family', 'PBE')
        self._paw_def = kwargs.get('paw_defaults')
        self._paws = {}
        self._set_default_paws()
        self._kpoints = kwargs.get('kpoints', self.calc_cls.new_kpoints())
        self.kpoints = self._kpoints
        self._charge_density = kwargs.get('charge_density', None)
        self._wavefunctions = kwargs.get('wavefunctions', None)
        self._recipe = None
        self._queue = kwargs.get('queue')
        self._resources = kwargs.get('resources', {})

    def _copy_from(self, calc):
        ins = calc.get_inputs_dict()
        if not self._calcname:
            self.calc_cls = calc.__class__
        self.label = calc.label + '_copy'
        self._computer = calc.get_computer()
        self._code = calc.get_code()
        self._settings = ins.get('settings')
        self._structure = ins.get('structure')
        self._paws = {}
        for paw in filter(lambda i: 'paw' in i[0], ins.iteritems()):
            self._paws[paw[0].replace('paw_', '')] = paw[1]
        self._kpoints = ins.get('kpoints')
        self._charge_density = ins.get('charge_density')
        self._wavefunctions = ins.get('wavefunctions')
        self._queue = calc.get_queue_name()
        self._resources = calc.get_resources()

    def _set_default_structure(self, structure):
        if not structure:
            self._structure = self.calc_cls.new_structure()
        elif isinstance(structure, (str, unicode)):
            structure = os.path.abspath(structure)
            if os.path.splitext(structure)[1] == '.cif':
                self._structure = DataFactory('cif').get_or_create(structure)[0]
            elif os.path.basename(structure) == 'POSCAR':
                from ase.io.vasp import read_vasp
                pwd = os.path.abspath(os.curdir)
                os.chdir(os.path.dirname(structure))
                atoms = read_vasp('POSCAR')
                os.chdir(pwd)
                self._structure = self.calc_cls.new_structure()
                self._structure.set_ase(atoms)
        else:
            self._structure = structure

    def _init_from(self, prev):
        self._copy_from(prev)
        if 'structure' in prev.get_outputs_dict():
            self.structure = prev.out.structure
        self.rewrite_settings(istart=1, icharg=11)
        self.wavefunctions = prev.out.wavefunctions
        self.charge_density = prev.out.charge_density

    def new(self):
        calc = self.calc_cls()
        calc.use_code(self._code)
        calc.use_structure(self._structure)
        for k in self.elements:
            calc.use_paw(self._paws[k], kind=k)
        calc.use_settings(self._settings)
        calc.use_kpoints(self._kpoints)
        calc.set_computer(self._computer)
        calc.set_queue_name(self._queue)
        if self._charge_density:
            calc.use_charge_density(self._charge_density)
        if self._wavefunctions:
            calc.use_wavefunctions(self._wavefunctions)
        calc.label = self.label
        calc.set_resources(self._resources)
        return calc

    @property
    def structure(self):
        return self._structure

    @structure.setter
    def structure(self, val):
        self._structure = val
        self._set_default_lda_paws()
        if self._kpoints.pk:
            self._kpoints = self._kpoints.copy()
        self._kpoints.set_cell(self._structure.get_ase().get_cell())

    @property
    def settings(self):
        return self._settings.get_dict()

    @property
    def kpoints(self):
        return self._kpoints

    @kpoints.setter
    def kpoints(self, kp):
        self._kpoints = kp
        self._kpoints.set_cell(self._structure.get_ase().get_cell())

    def set_kpoints_path(self, value=None, weights=None, **kwargs):
        '''
        py:method:: set_kpoints_path([value=None][, weights=None], **kwargs)

        Calls kpoints' set_kpoints_path method with value, automatically adds weights.
        Copies the kpoints node if it's already stored.
        '''
        if self._kpoints.pk:
            self.kpoints = self._kpoints.copy()
        self._kpoints.set_kpoints_path(value=value, **kwargs)
        if 'weights' not in kwargs:
            kpl = self._kpoints.get_kpoints()
            wl = [1. for i in kpl]
            self._kpoints.set_kpoints(kpl, weights=wl)

    def set_kpoints_mesh(self, *args, **kwargs):
        '''
        py:method:: set_kpoints_mesh(*args, **kwargs)

        Passes arguments on to kpoints.set_kpoints_mesh, copies if it was already stored.
        '''
        if self._kpoints.pk:
            self.kpoints = self.calc_cls.new_kpoints()
        self._kpoints.set_kpoints_mesh(*args, **kwargs)

    def set_kpoints_list(self, *args, **kwargs):
        '''
        py:method:: set_kpoints_list(*args, **kwargs)

        Passes arguments on to kpoints.set_kpoints, copies if it was already stored.
        '''
        if self._kpoints.pk:
            self.kpoints = self.calc_cls.new_kpoints()
        self._kpoints.set_kpoints(*args, **kwargs)

    @property
    def wavefunctions(self):
        return self._wavefunctions

    @wavefunctions.setter
    def wavefunctions(self, val):
        self._wavefunctions = val
        self.add_settings(istart=1)

    @property
    def charge_density(self):
        return self._charge_density

    @charge_density.setter
    def charge_density(self, val):
        self._charge_density = val
        self.add_settings(icharg=11)

    @property
    def code(self):
        return self._code

    @code.setter
    def code(self, val):
        self._code = val
        self._computer = val.get_computer()

    @property
    def computer(self):
        return self._computer

    @computer.setter
    def computer(self, val):
        self._computer = val

    @property
    def queue(self):
        return self._queue

    @queue.setter
    def queue(self, val):
        self._queue = val

    @property
    def resources(self):
        return self._resources

    def add_settings(self, **kwargs):
        if self._settings.pk:
            self._settings = self._settings.copy()
        for k, v in kwargs.iteritems():
            if k not in self.settings:
                self._settings.update_dict({k: v})

    def rewrite_settings(self, **kwargs):
        if self._settings.pk:
            self._settings = self._settings.copy()
        self._settings.update_dict(kwargs)

    def _set_default_paws(self, overwrite=False):
        if self._paw_fam.lower() == 'lda':
            defaults = self._paw_def or lda
        elif self._paw_fam.lower() in ['pbe', 'gw']:
            defaults = self._paw_def or gw
        else:
            if not self._paw_def:
                print 'keyword paw_family was not LDA or PBE'
                print 'and no paw_defaults keyword was given!'
                print 'manual paw initialization required'
                return None
            else:
                defaults = self._paw_def
        for k in self.elements:
            if k not in self._paws or overwrite:
                paw = self.calc_cls.Paw.load_paw(
                    family=self._paw_fam, symbol=defaults[k])[0]
                self._paws[k] = paw

    @property
    def elements(self):
        return set(self._structure.get_ase().get_chemical_symbols())

    def pkcmp(self, nodeA, nodeB):
        if nodeA.pk < nodeB.pk:
            return -1
        elif nodeA.pk > nodeB.pk:
            return 1
        else:
            return 0

    def verify_incar(self):
        if not self.struct:
            raise ValueError('need structure,')
        magmom = self.incar.get('magmom', [])
        lsorb = self.incar.get('lsorbit', False)
        lnonc = self.incar.get('lnoncollinear', False)
        ok = True
        nmag = len(magmom)
        nsit = len(self.struct.sites)
        if lsorb:
            if lnonc:
                if magmom and not nmag == 3*nsit:
                    ok = False
            else:
                if not nmag == nsit:
                    ok = False
        else:
            if not nmag == nsit:
                ok = False
        return ok

    def check_magmom(self):
        magmom = self.settings.get('magmom', [])
        st_magmom = self._structure.get_ase().get_initial_magnetic_moments()
        lsf = self.noncol and 3 or 1
        nio = self.n_ions
        s_mm = nio * lsf
        mm = len(magmom)
        if magmom and st_magmom:
            return s_mm == mm

    def set_magmom_1(self, val):
        magmom = [val]
        magmom *= self.n_ions
        magmom *= self.noncol and 3 or 1
        self.rewrite_settings(magmom=magmom)

    @property
    def nbands(self):
        return self.n_ions * 3 * self.noncol and 3 or 1

    @property
    def n_ions(self):
        return self.structure.get_ase().get_number_of_atoms()

    @property
    def n_elec(self):
        res = 0
        for k in self._structure.get_ase().get_chemical_symbols():
            res += self._paws[k].valence
        return res

    @property
    def noncol(self):
        lsorb = self.settings.get('lsorbit', False)
        lnonc = self.settings.get('lnoncollinear', False)
        return lsorb or lnonc

    @property
    def icharg(self):
        return self.settings['icharg']

    @icharg.setter
    def icharg(self, value):
        if value not in [0, 1, 2, 4, 10, 11, 12]:
            raise ValueError('invalid ICHARG value for vasp 5.3.5')
        else:
            self.settings['icharg'] = value

    @property
    def recipe(self):
        return self._recipe

    @recipe.setter
    def recipe(self, val):
        if self._recipe and self._recipe != val:
            raise ValueError('recipe is already set to something else')
        self._init_recipe(val)
        self._recipe = val

    def _init_recipe(self, recipe):
        if recipe == 'test_sc':
            self._init_recipe_test_sc()
        else:
            raise ValueError('recipe not recognized')

    def _init_recipe_test_sc(self):
        self.add_settings(
            gga='PE',
            gga_compat=False,
            ismear=0,
            lorbit=11,
            lsorbit=True,
            sigma=0.05,
        )
