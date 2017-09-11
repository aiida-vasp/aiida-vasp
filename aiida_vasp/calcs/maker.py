"""
A utility class to simplify creating VASP-Calculations, loosely follows the builder design pattern
"""
import os

from aiida.orm import CalculationFactory, DataFactory

# pylint: disable=too-many-instance-attributes,too-many-public-methods
from aiida_vasp.calcs.vasp import ordered_unique_list


class VaspMaker(object):
    """
    simplifies creating a Scf, Nscf or AmnCalculation from scratch interactively or
    as a copy or continuation of a previous calculation
    further simplifies creating certain often used types of calculations

    Most of the required information can be given as keyword arguments to
    the constructor or set via properties later on.

    The input information is stored in the instance and the calculation
    is only built in the :py:meth:`new` method. This also makes it possible
    to create a set of similar calculations in an interactive setting very
    quickly.

    :param structure: A StructureData node or a
        (relative) path to either a .cif file or a POSCAR file. Defaults to
        a new empty structure node recieved from calc_cls.
    :type structure: str or StructureData

    :keyword calc_cls: the class that VaspMaker will use when creating
        Calculation nodes.
        defaults to 'vasp.vasp'.
        if a string is given, it will be passed to aiida's CalculationFactory
    :type calc_cls: str or vasp.BasicCalculation subclass

    :keyword continue_from: A vasp calculation node with charge_density and
        wavefunction output links. VaspMaker will create calculations that
        start with those as inputs.
    :type continue_from: vasp calculation node

    :keyword copy_from: A vasp calculation. It's inputs will be used as
        defaults for the created calculations.
    :type copy_from: vasp calculation node

    :keyword charge_density: chargedensity node from a previously run
        calculation
    :type charge_density: ChargedensityData
    :keyword wavefunctions: wavefunctions node from a previously run
        calculation
    :type wavefunctions: WavefunData
    :keyword array.KpointsData kpoints: kpoints node to use for input
    :keyword str paw_family: The name of a PAW family stored in the db
    :keyword str paw_map: A dictionary mapping element symbols -> PAW
        symbols
    :keyword str label: value for the calculation label
    :keyword str computer: computer name, defaults to code's if code is
        given
    :keyword str code: code name, if any Calculations are given, defaults
        to their code
    :keyword str resources: defaults to copy_from.get_resources() or None
    :keyword str queue: defaults to queue from given calculation, if any,
        or None

    .. py:method:: new()

        :returns: an instance of :py:attr:`calc_cls`, initialized with the data
        held by the VaspMaker

    .. py:method:: add_parameters(**kwargs)

        Adds keys to the parameters (INCAR keywords), if parameters is already
        stored, makes a copy.
        Does not overwrite previously set keywords.

    .. py:method:: rewrite_parameters(**kwargs)

        Same as :py:meth:`add_parameters`, but also overwrites keywords.

    .. py:attribute:: structure

        Used to initialize the created calculations as well as other nodes
        (like kpoints).
        When changed, can trigger changes in other data nodes.

    .. py:attribute:: calc_cls

        Vasp calculation class to be used in :py:meth:`new`

    .. py:attribute:: computer

    .. py:attribute:: code

    .. py:attribute:: queue

    .. py:attribute:: parameters

        A readonly shortcut to the contents of the parameters node

    .. py:attribute:: kpoints

        The kpoints node to be used, may be copied to have py:func:set_cell
        called.

    .. py:attribute:: wavefunction

    .. py:attribute:: charge_density

    .. py:attribute:: elements

        Chemical symbols of the elements contained in py:attr:structure
    """

    def __init__(self, *args, **kwargs):
        self._recipe = None
        self._parameters = None
        self._queue = None
        self._computer = None
        self._code = None
        self._wannier_data = None
        self._charge_density = None
        self._wavefunctions = None
        self._wannier_parameters = None
        self._kpoints = None
        self._structure = None

        self._init_defaults(*args, **kwargs)
        self._calcname = kwargs.get('calc_cls')
        if 'continue_from' in kwargs:
            self._init_from(kwargs['continue_from'])
        if 'copy_from' in kwargs:
            self._copy_from(kwargs['copy_from'])

    def _init_defaults(self, *args, **kwargs):  # pylint: disable=unused-argument
        """Set default values"""
        calcname = kwargs.get('calc_cls', 'vasp.vasp')
        if isinstance(calcname, (str, unicode)):
            self.calc_cls = CalculationFactory(calcname)
        else:
            self.calc_cls = calcname
        self.label = kwargs.get('label', 'unlabeled')
        self._computer = kwargs.get('computer')
        self._code = kwargs.get('code')
        self._parameters = kwargs.get('parameters',
                                      self.calc_cls.new_parameters())
        self._set_default_structure(kwargs.get('structure'))
        self._paw_fam = kwargs.get('paw_family', 'PBE')
        self._paw_def = kwargs.get('paw_map')
        self._paws = {}
        self._set_default_paws()
        self._kpoints = kwargs.get('kpoints', self.calc_cls.new_kpoints())
        self.kpoints = self._kpoints
        self._charge_density = kwargs.get('charge_density', None)
        self._wavefunctions = kwargs.get('wavefunctions', None)
        self._wannier_parameters = kwargs.get('wannier_parameters', None)
        self._wannier_data = kwargs.get('wannier_data', None)
        self._recipe = None
        self._queue = kwargs.get('queue')
        self._resources = kwargs.get('resources', {})

    def _copy_from(self, calc):
        """Copy data links from another calculation"""
        ins = calc.get_inputs_dict()
        if not self._calcname:
            self.calc_cls = calc.__class__
        self.label = calc.label + '_copy'
        self._computer = calc.get_computer()
        self._code = calc.get_code()
        self._parameters = ins.get('parameters')
        self._structure = ins.get('structure')
        self._paws = {}
        for paw in [i for i in ins.iteritems() if 'paw' in i[0]]:
            self._paws[paw[0].replace('paw_', '')] = paw[1]
        self._kpoints = ins.get('kpoints')
        self._charge_density = ins.get('charge_density')
        self._wavefunctions = ins.get('wavefunctions')
        self._wannier_parameters = ins.get('wannier_parameters')
        self._wannier_data = ins.get('wannier_data')
        self._queue = calc.get_queue_name()
        self._resources = calc.get_resources()

    def _set_default_structure(self, structure):
        """Set a structure depending on what was given, empty if nothing was given"""
        if not structure:
            self._structure = self.calc_cls.new_structure()
        elif isinstance(structure, (str, unicode)):
            structure = os.path.abspath(os.path.expanduser(structure))
            if os.path.splitext(structure)[1] == '.cif':
                self._structure = DataFactory('cif').get_or_create(structure)[
                    0]
            elif os.path.basename(structure) == 'POSCAR':
                from ase.io.vasp import read_vasp
                atoms = read_vasp(os.path.abspath(structure))
                self._structure = self.calc_cls.new_structure()
                self._structure.set_ase(atoms)
        else:
            self._structure = structure

    def _init_from(self, prev):
        """Initialize from an already run calculation"""
        out = prev.get_outputs_dict()
        self._copy_from(prev)
        if 'structure' in out:
            self.structure = prev.out.structure
        self.rewrite_parameters(istart=1, icharg=11)
        self.wavefunctions = prev.out.wavefunctions
        self.charge_density = prev.out.charge_density
        self._wannier_parameters = out.get('wannier_parameters',
                                           self._wannier_parameters)
        self._wannier_data = out.get('wannier_data', self.wannier_data)

    def new(self):
        """Create a new (unstored) Calculation node from previously set properties"""
        calc = self.calc_cls()
        calc.use_code(self._code)
        calc.use_structure(self._structure)
        for k in self.elements:
            calc.use_paw(self._paws[k], kind=k)
        calc.use_parameters(self._parameters)
        calc.use_kpoints(self._kpoints)
        calc.set_computer(self._computer)
        calc.set_queue_name(self._queue)
        if self._charge_density:
            calc.use_charge_density(self._charge_density)
        if self._wavefunctions:
            calc.use_wavefunctions(self._wavefunctions)
        if self._wannier_parameters:
            calc.use_wannier_parameters(self._wannier_parameters)
        if self._wannier_data:
            calc.use_wannier_data(self._wannier_data)
        calc.label = self.label
        calc.set_resources(self._resources)
        return calc

    @property
    def structure(self):
        return self._structure

    @structure.setter
    def structure(self, val):
        self._set_default_structure(val)
        self._set_default_paws()
        if self._kpoints.pk:
            self._kpoints = self._kpoints.copy()
        self._kpoints.set_cell(self._structure.get_ase().get_cell())

    @property
    def parameters(self):
        return self._parameters.get_dict()

    @property
    def kpoints(self):
        return self._kpoints

    @kpoints.setter
    def kpoints(self, kpoints):
        self._kpoints = kpoints
        self._kpoints.set_cell(self._structure.get_ase().get_cell())

    def set_kpoints_path(self, value=None, weights=None, **kwargs):
        '''
        Calls kpoints' set_kpoints_path method with value, automatically adds
        weights.
        Copies the kpoints node if it's already stored.
        '''
        if self._kpoints.is_stored:
            self.kpoints = self.calc_cls.new_kpoints()
        self._kpoints.set_kpoints_path(value=value, **kwargs)
        if not weights:
            kp_list = self._kpoints.get_kpoints()
            weights = [1. for _ in kp_list]
            self._kpoints.set_kpoints(kp_list, weights=weights)

    def set_kpoints_mesh(self, *args, **kwargs):
        '''
        Passes arguments on to kpoints.set_kpoints_mesh, copies if it was
        already stored.
        '''
        if self._kpoints.pk:
            self.kpoints = self.calc_cls.new_kpoints()
        self._kpoints.set_kpoints_mesh(*args, **kwargs)

    def set_kpoints_list(self, kpoints, weights=None, **kwargs):
        '''
        Passes arguments on to kpoints.set_kpoints, copies if it was already
        stored.
        '''
        import numpy as np
        if self._kpoints.pk:
            self.kpoints = self.calc_cls.new_kpoints()
        if not weights:
            weights = np.ones(len(kpoints), dtype=float)
        self._kpoints.set_kpoints(kpoints, weights=weights, **kwargs)

    @property
    def wavefunctions(self):
        return self._wavefunctions

    @wavefunctions.setter
    def wavefunctions(self, val):
        self._wavefunctions = val
        self.add_parameters(istart=1)

    @property
    def charge_density(self):
        return self._charge_density

    @charge_density.setter
    def charge_density(self, val):
        self._charge_density = val
        self.add_parameters(icharg=11)

    @property
    def wannier_parameters(self):
        return self._wannier_parameters

    @wannier_parameters.setter
    def wannier_parameters(self, val):
        self._wannier_parameters = val
        if 'lwannier90' not in self.parameters:
            self.add_parameters(lwannier90=True)

    @property
    def wannier_data(self):
        return self._wannier_data

    @wannier_data.setter
    def wannier_data(self, val):
        self._wannier_data = val

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

    @resources.setter
    def resources(self, val):
        if isinstance(val, dict):
            self._resources.update(val)
        else:
            self._resources['num_machines'] = val[0]
            self._resources['num_mpiprocs_per_machine'] = val[1]

    def add_parameters(self, **kwargs):
        if self._parameters.pk:
            self._parameters = self._parameters.copy()
        for key, value in kwargs.iteritems():
            if key not in self.parameters:
                self._parameters.update_dict({key: value})

    def rewrite_parameters(self, **kwargs):
        if self._parameters_conflict(kwargs):
            if self._parameters.pk:
                self._parameters = self._parameters.copy()
            self._parameters.update_dict(kwargs)

    def _parameters_conflict(self, parameters):
        conflict = False
        for key, value in parameters.iteritems():
            conflict |= (self.parameters.get(key) != value)
        return conflict

    def _set_default_paws(self):
        for key in self.elements:
            if key not in self._paws:
                if self._paw_def is None:
                    raise ValueError(
                        "The 'paw_map' keyword is required. Pre-defined potential "
                        "mappings are defined in 'aiida.tools.codespecific.vasp.default_paws'."
                    )
                try:
                    paw = self.calc_cls.Paw.load_paw(
                        family=self._paw_fam, symbol=self._paw_def[key])[0]
                except KeyError:
                    raise ValueError(
                        "The given 'paw_map' does not contain a mapping for element '{}'".
                        format(key))
                self._paws[key] = paw

    @property
    def elements(self):
        return ordered_unique_list(
            self._structure.get_ase().get_chemical_symbols())

    @staticmethod
    def compare_pk(node_a, node_b):
        """Compare two nodes by primary key"""
        if node_a.pk < node_b.pk:
            return -1
        elif node_a.pk > node_b.pk:
            return 1
        return 0

    def verify_parameters(self):
        """
        Verify input parameters

        :return: (bool, string) success, message
        """
        if not self._structure:
            raise ValueError('need structure,')
        magmom = self.parameters.get('magmom', [])
        lsorb = self.parameters.get('lsorbit', False)
        lnonc = self.parameters.get('lnoncollinear', False)
        parameters_ok = True
        msg = 'Everything ok'
        nmag = len(magmom)
        nsit = self.n_ions
        if lsorb:
            if lnonc:
                if magmom and not nmag == 3 * nsit:
                    parameters_ok = False
                    msg = 'magmom has wrong dimension'
            else:
                if magmom and not nmag == nsit:
                    parameters_ok = False
                    msg = 'magmom has wrong dimension'
        else:
            if magmom and not nmag == nsit:
                parameters_ok = False
                msg = 'magmom has wrong dimension'
        return parameters_ok, msg

    def check_magmom(self):
        """Check that the magnetic moment given in parameters matches the one of the structure"""
        magnetic_moment = self.parameters.get('magmom', [])
        structure_mm = self._structure.get_ase().get_initial_magnetic_moments()
        lsf = 3 if self.noncol else 1
        nio = self.n_ions
        structure_mm_dim = nio * lsf
        magnetic_moment_dim = len(magnetic_moment)
        if magnetic_moment and structure_mm:
            return structure_mm_dim == magnetic_moment_dim
        return True

    def set_magmom_1(self, val):
        magmom = [val]
        magmom *= self.n_ions
        magmom *= self.noncol and 3 or 1
        self.rewrite_parameters(magmom=magmom)

    @property
    def nbands(self):
        return self.n_ions * 3 * (self.noncol and 3 or 1)

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
        lsorb = self.parameters.get('lsorbit', False)
        lnonc = self.parameters.get('lnoncollinear', False)
        return lsorb or lnonc

    @property
    def icharg(self):
        return self.parameters.get('icharg', 'default')

    @icharg.setter
    def icharg(self, value):
        if value not in [0, 1, 2, 4, 10, 11, 12]:
            raise ValueError('invalid ICHARG value for vasp 5.3.5')
        else:
            self.parameters['icharg'] = value

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
        self.add_parameters(
            gga='PE',
            gga_compat=False,
            ismear=0,
            lorbit=11,
            lsorbit=True,
            sigma=0.05, )
