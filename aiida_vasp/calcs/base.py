try:
    from collections import ChainMap
except ImportError:
    from chainmap import ChainMap

from aiida.orm import JobCalculation, DataFactory
from aiida.common.utils import classproperty
from aiida.common.datastructures import CalcInfo, CodeInfo


def make_use_methods(inputs, bases):
    '''
    creates the _use_methods classproperty from a list of
    base classes as well as :py:class:`Input` instances. For use
    in :py:class:`CalcMeta` during class creation.
    '''
    @classproperty
    def _use_methods(cls):
        retdict = JobCalculation._use_methods
        for b in bases:
            if hasattr(b, '_use_methods'):
                retdict.update(b._use_methods)
        for inp, dct in inputs.iteritems():
            ln = dct['linkname']
            if isinstance(ln, classmethod):
                dct['linkname'] = getattr(cls, ln.__func__.__name__)
        retdict.update(inputs)
        return retdict
    return _use_methods


def make_init_internal(cls, **kwargs):
    '''
    returns a _update internal_params method, required in class creation
    by :py:class:`CalcMeta`
    '''
    def _update_internal_params(self):
        for k, v in kwargs.iteritems():
            setattr(self, k, v)
    return _update_internal_params


def seqify(seq_arg):
    if not isinstance(seq_arg, list) and not isinstance(seq_arg, tuple):
        seq_arg = [seq_arg]
    return tuple(seq_arg)


def datify(type_or_str):
    if isinstance(type_or_str, str):
        return DataFactory(type_or_str)
    else:
        return type_or_str


class Input(object):
    '''
    Utility class that handles mapping between CalcMeta's and Calculation's
    interfaces for defining input nodes
    '''
    def __init__(self, types, param=None, ln=None, doc=''):
        self.types = tuple(map(datify, seqify(types)))
        self.param = param
        self._ln = ln
        self.doc = doc

    def set_ln(self, value):
        if not self._ln:
            self._ln = value

    def get_dict(self):
        ret = {
            'valid_types': self.types,
            'additional_parameter': self.param,
            'linkname': self._ln,
            'docstring': self.doc
        }
        return ret

    @classmethod
    def it_filter(cls, item):
        return isinstance(item[1], Input)

    @classmethod
    def filter_classdict(cls, classdict):
        return filter(cls.it_filter, classdict.iteritems())


class IntParam(object):
    '''
    Utility class that handles mapping between CalcMeta's and Calculation's
    internal parameter interfaces
    '''

    pmap = {
        'default_parser': '_default_parser',
        'input_file_name': '_INPUT_FILE_NAME',
        'output_file_name': '_OUTPUT_FILE_NAME'
    }

    @classmethod
    def it_filter(cls, item):
        return item[0] in cls.pmap

    @classmethod
    def filter_classdict(cls, classdict):
        return filter(cls.it_filter, classdict.iteritems())

    @classmethod
    def map_param(cls, item):
        k, v = item
        return (cls.pmap[k], v)

    @classmethod
    def map_params(cls, classdict):
        return map(cls.map_param, cls.filter_classdict(classdict))

    @classmethod
    def k_filter(cls, key):
        return key in cls.pmap

    @classmethod
    def get_keylist(cls, classdict):
        return filter(cls.k_filter, classdict)


class CalcMeta(JobCalculation.__metaclass__):
    '''
    Metaclass that allows simpler and clearer Calculation class writing.
    Takes :py:class:`Input` instances from the class and converts it to
    the correct entries in the finished class's _use_methods classproperty.
    finds class attributes corresponding to 'internal parameters'
    and creates the finished class's _update_internal_params method.
    '''
    def __new__(cls, name, bases, classdict):
        inputs = {}
        delete = []
        inputobj = Input.filter_classdict(classdict)
        for k, v in inputobj:
            if not v.param:
                v.set_ln(k)
            else:
                v.set_ln(classdict['_get_{}_linkname'.format(k)])
            inputs[k] = v.get_dict()
            delete.append(k)
        internals = dict(IntParam.map_params(classdict))
        delete.extend(IntParam.get_keylist(classdict))
        classdict['_use_methods'] = make_use_methods(inputs, bases)
        classdict['_update_internal_params'] = make_init_internal(
            cls, **internals)
        for k in delete:
            classdict.pop(k)
        Calc = super(CalcMeta, cls).__new__(cls, name, bases, classdict)
        return Calc


class VaspCalcBase(JobCalculation):
    '''
    Base class of all calculations utilizing VASP.

    * sets :py:class:`CalcMeta` as it's __metaclass__
    * Defines internal parameters common to all vasp calculations.
    * provides a basic, extendable implementation of _prepare_for_submission
    * provides hooks, so subclasses can extend the behaviour without

    having to reimplement common functionality
    '''
    __metaclass__ = CalcMeta
    input_file_name = 'INCAR'
    output_file_name = 'OUTCAR'

    @classmethod
    def max_retrieve_list(cls):
        '''returns a list of all possible output files from a VASP run'''
        retrieve_list = [
            'CHG',
            'CHGCAR',
            'CONTCAR',
            'DOSCAR',
            'EIGENVAL',
            'ELFCAR',
            'IBZKPT',
            'LOCPOT',
            'OSZICAR',
            'OUTCAR',
            'PCDAT',
            'PROCAR',
            'PROOUT',
            'STOPCAR',
            'TMPCAR',
            'WAVECAR',
            'XDATCAR',
            ('wannier90*', '.', 0),
            'vasprun.xml'
        ]
        return retrieve_list

    def _prepare_for_submission(self, tempfolder, inputdict):
        '''
        Writes the four minimum output files,
        INCAR, POSCAR, POTCAR, KPOINTS. Delegates the
        construction and writing / copying to write_<file> methods.
        That way, subclasses can use any form of input nodes and just
        have to implement the write_xxx method accordingly.
        Subclasses can extend by calling the super method and if neccessary
        modifying it's output CalcInfo before returning it.
        '''
        # write input files
        incar = tempfolder.get_abs_path('INCAR')
        structure = tempfolder.get_abs_path('POSCAR')
        paw = tempfolder.get_abs_path('POTCAR')
        kpoints = tempfolder.get_abs_path('KPOINTS')

        self.verify_inputs(inputdict)
        self.write_incar(inputdict, incar)
        self.write_poscar(inputdict, structure)
        self.write_potcar(inputdict, paw)
        self.write_kpoints(inputdict, kpoints)
        self.write_additional(tempfolder, inputdict)

        # calcinfo
        calcinfo = CalcInfo()
        calcinfo.uuid = self.uuid
        calcinfo.retrieve_list = self.max_retrieve_list()
        codeinfo = CodeInfo()
        codeinfo.code_uuid = self.get_code().uuid
        codeinfo.code_pk = self.get_code().pk
        calcinfo.codes_info = [codeinfo]

        return calcinfo

    def _init_internal_params(self):
        '''must be present on all JobCalculation subclasses, that
        set internal parameters'''
        super(VaspCalcBase, self)._init_internal_params()
        self._update_internal_params()

    def verify_inputs(self, inputdict, *args, **kwargs):
        '''
        Hook to be extended by subclasses with checks for input nodes.
        Is called once before submission.
        '''
        self.check_input(inputdict, 'code')
        return True

    def check_input(self, inputdict, linkname, check_fn=lambda: True):
        '''
        :py:classmethod:: check_input(inputdict, linkname[, check_fn])

        check wether the given linkname is in the inputdict.

        This is meant to be called prior to submission and will raise an
        Exception if no input node is found for the linkname.

        :param function check_fn: A function that returns True if the
            check should be performed or False if not.

        '''
        notset_msg = 'input not set: %s'
        if check_fn():
            if linkname not in inputdict:
                raise ValueError(notset_msg % linkname)
        return True

    def store(self, *args, **kwargs):
        '''adds a _prestore subclass hook for operations that
        should be done just before storing'''
        self._prestore()
        super(VaspCalcBase, self).store(*args, **kwargs)

    def _prestore(self):
        '''Subclass hook for updating attributes etc, just before storing'''
        pass

    def write_additional(self, tempfolder, inputdict):
        """
        Subclass hook to write additional input files.
        """
        pass
