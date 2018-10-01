# pylint: disable=abstract-method,invalid-metaclass,ungrouped-imports
# explanation: pylint wrongly complains about Node not implementing query
"""Base and meta classes for VASP calculations"""
import os
import six
from py import path as py_path  # pylint: disable=no-name-in-module,no-member

from aiida.orm import JobCalculation, DataFactory
from aiida.common.utils import classproperty
from aiida.common.datastructures import CalcInfo, CodeInfo
from aiida.common.exceptions import ValidationError
from aiida.common.folders import SandboxFolder

from aiida_vasp.utils.aiida_utils import get_data_node, cmp_get_transport

try:
    from aiida.orm.implementation.general.node import _AbstractNodeMeta as __absnode__
except ImportError:
    __absnode__ = JobCalculation.__metaclass__


def make_use_methods(inputs, bases):
    """
    Create the _use_methods classproperty from a list of base classes as well as :py:class:`Input` instances.

    For use in :py:class:`CalcMeta` during class creation.
    """

    @classproperty
    # pylint: disable=protected-access,no-member
    def _use_methods(cls):
        """Automatically generated _use_methods function."""
        retdict = JobCalculation._use_methods
        for base_class in bases:
            if hasattr(base_class, '_use_methods'):
                retdict.update(base_class._use_methods)
        for _, dct in inputs.iteritems():
            link_name = dct['linkname']
            if isinstance(link_name, classmethod):
                dct['linkname'] = getattr(cls, link_name.__func__.__name__)
        retdict.update(inputs)
        return retdict

    return _use_methods


def make_init_internal(cls, **kwargs):  # pylint: disable=unused-argument
    """Returns a _update internal_params method, required in class creation by :py:class:`CalcMeta`."""

    def _update_internal_params(self):
        for key, value in kwargs.iteritems():
            setattr(self, key, value)

    return _update_internal_params


def seqify(seq_arg):
    if not isinstance(seq_arg, list) and not isinstance(seq_arg, tuple):
        seq_arg = [seq_arg]
    return tuple(seq_arg)


def datify(type_or_str):
    if isinstance(type_or_str, str):
        return DataFactory(type_or_str)
    return type_or_str


class Input(object):
    """
    Utility class that handles mapping between CalcMeta's and Calculation's interfaces for defining input nodes.

    Examples::

        parameters = Input(types='parameter', doc='input parameters.')
        structure = Input(types=['structure', 'cif'], doc='input structure')
        potential = Input(types='vasp.potcar', param='kind')

    Usage::

        @six.add_metaclass(CalcMeta)
        MyCalculation(JobCalculation):
            potential = Input(types='vasp.potcar', param='kind')

        my_calc = MyCalculation
        potential = load_node(...)
        my_calc.use_potential(potential, kind='In')

    """

    def __init__(self, types, param=None, ln=None, doc=''):
        self.types = tuple(map(datify, seqify(types)))
        self.param = param
        self._ln = ln
        self.doc = doc

    def set_ln(self, value):
        if not self._ln:
            self._ln = value

    def get_dict(self):
        """Create a dictionary representation of this input"""
        ret = {'valid_types': self.types, 'additional_parameter': self.param, 'linkname': self._ln, 'docstring': self.doc}
        return ret

    @classmethod
    def it_filter(cls, item):
        return isinstance(item[1], Input)

    @classmethod
    def filter_classdict(cls, classdict):
        return filter(cls.it_filter, classdict.iteritems())


class IntParam(object):
    """Utility class that handles mapping between CalcMeta's and Calculation's internal parameter interfaces."""

    pmap = {'default_parser': '_default_parser', 'input_file_name': '_INPUT_FILE_NAME', 'output_file_name': '_OUTPUT_FILE_NAME'}

    @classmethod
    def it_filter(cls, item):
        return item[0] in cls.pmap

    @classmethod
    def filter_classdict(cls, classdict):
        return filter(cls.it_filter, classdict.iteritems())

    @classmethod
    def map_param(cls, item):
        key, value = item
        return (cls.pmap[key], value)

    @classmethod
    def map_params(cls, classdict):
        return map(cls.map_param, cls.filter_classdict(classdict))

    @classmethod
    def k_filter(cls, key):
        return key in cls.pmap

    @classmethod
    def get_keylist(cls, classdict):
        return filter(cls.k_filter, classdict)


class CalcMeta(__absnode__):
    """
    Metaclass that allows simpler and clearer Calculation class writing.

    Takes :py:class:`Input` instances from the class and converts it to
    the correct entries in the finished class's _use_methods classproperty.
    finds class attributes corresponding to 'internal parameters'
    and creates the finished class's _update_internal_params method.
    """

    def __new__(mcs, name, bases, classdict):
        inputs = {}
        delete = []
        inputobj = Input.filter_classdict(classdict)
        for key, value in inputobj:
            if not value.param:
                value.set_ln(key)
            else:
                value.set_ln(classdict['_get_{}_linkname'.format(key)])
            inputs[key] = value.get_dict()
            delete.append(key)
        internals = dict(IntParam.map_params(classdict))
        delete.extend(IntParam.get_keylist(classdict))
        classdict['_use_methods'] = make_use_methods(inputs, bases)
        classdict['_update_internal_params'] = make_init_internal(mcs, **internals)
        for key in delete:
            classdict.pop(key)
        calc_cls = super(CalcMeta, mcs).__new__(mcs, name, bases, classdict)
        return calc_cls


@six.add_metaclass(CalcMeta)
class VaspCalcBase(JobCalculation):
    """
    Base class of all calculations utilizing VASP.

    * sets :py:class:`CalcMeta` as it's metaclass
    * Defines internal parameters common to all vasp calculations.
    * provides a basic, extendable implementation of _prepare_for_submission
    * provides hooks, so subclasses can extend the behaviour without

    having to reimplement common functionality
    """

    input_file_name = 'INCAR'
    output_file_name = 'OUTCAR'

    restart_folder = Input(types='remote', doc='A remote folder to restart from after crashing')

    @classmethod
    def max_retrieve_list(cls):
        """Return a list of all possible output files from a VASP run."""
        retrieve_list = [
            'CHG', 'CHGCAR', 'CONTCAR', 'DOSCAR', 'EIGENVAL', 'ELFCAR', 'IBZKPT', 'LOCPOT', 'OSZICAR', 'OUTCAR', 'PCDAT', 'PROCAR',
            'PROOUT', 'STOPCAR', 'TMPCAR', 'WAVECAR', 'XDATCAR', ['wannier90*', '.', 0], 'vasprun.xml'
        ]
        return retrieve_list

    def _prepare_for_submission(self, tempfolder, inputdict):
        """
        Writes the four minimum output files, INCAR, POSCAR, POTCAR, KPOINTS.

        Delegates the construction and writing / copying to write_<file> methods.
        That way, subclasses can use any form of input nodes and just
        have to implement the write_xxx method accordingly.
        Subclasses can extend by calling the super method and if neccessary
        modifying it's output CalcInfo before returning it.
        """
        # write input files
        incar = tempfolder.get_abs_path('INCAR')
        structure = tempfolder.get_abs_path('POSCAR')
        potentials = tempfolder.get_abs_path('POTCAR')
        kpoints = tempfolder.get_abs_path('KPOINTS')

        remote_copy_list = []

        self.verify_inputs(inputdict)
        if self._is_restart(inputdict):
            remote_copy_list.extend(self.remote_copy_restart_folder(inputdict))
        self.write_incar(inputdict, incar)
        self.write_poscar(inputdict, structure)
        self.write_potcar(inputdict, potentials)
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
        calcinfo.remote_copy_list = remote_copy_list

        return calcinfo

    def _init_internal_params(self):
        """Must be present on all JobCalculation subclasses, that set internal parameters."""
        super(VaspCalcBase, self)._init_internal_params()
        self._update_internal_params()

    def remote_copy_restart_folder(self, inputdict):
        """Add all files required for restart to the list of files to be copied from the previous calculation."""
        restart_folder = inputdict['restart_folder']
        computer = self.get_computer()
        excluded = ['INCAR', '_aiidasubmit.sh', '.aiida']
        copy_list = [(computer.uuid, os.path.join(restart_folder.get_remote_path(), name), '.')
                     for name in restart_folder.listdir()
                     if name not in excluded]
        return copy_list

    def verify_inputs(self, inputdict, *args, **kwargs):  # pylint: disable=unused-argument
        """
        Hook to be extended by subclasses with checks for input nodes.

        Is called once before submission.
        """
        self.check_input(inputdict, 'code')
        self.check_restart_folder(inputdict)
        return True

    @staticmethod
    def check_input(inputdict, linkname, check_fn=lambda: True):
        """
        Check wether the given linkname is in the inputdict.

        This is meant to be called prior to submission and will raise an
        Exception if no input node is found for the linkname.

        :param function check_fn: A function that returns True if the
            check should be performed or False if not.

        """
        notset_msg = 'input not set: %s'
        if check_fn():
            if linkname not in inputdict:
                raise ValidationError(notset_msg % linkname)
        return True

    def check_restart_folder(self, inputdict):
        restart_folder = inputdict.get('restart_folder', None)
        if restart_folder:
            previous_calc = restart_folder.get_inputs(node_type=JobCalculation)[0]
            if not self.get_computer().pk == previous_calc.get_computer().pk:
                raise ValidationError('Calculation can not be restarted on another computer')

    # pylint: disable=no-self-use
    def _is_restart(self, inputdict):
        restart_folder = inputdict.get('restart_folder', None)
        is_restart = False
        if restart_folder:
            is_restart = True
        return is_restart

    def store(self, *args, **kwargs):
        """Adds a _prestore subclass hook for operations that should be done just before storing."""
        self._prestore()
        super(VaspCalcBase, self).store(*args, **kwargs)

    def _prestore(self):
        """Subclass hook for updating attributes etc, just before storing"""
        pass

    def write_additional(self, tempfolder, inputdict):
        """Subclass hook to write additional input files."""
        pass

    @classmethod
    def immigrant(cls, code, remote_path, **kwargs):
        """
        Create an immigrant with appropriate inputs from a code and a remote path on the associated computer.

        More inputs are required to pass resources information, if the POTCAR file is missing from the folder
        or if additional settings need to be passed, e.g. parser instructions.

        :param code: a Code instance for the code originally used.
        :param remote_path: The directory on the code's computer in which the simulation was run.
        :param resources: dict. The resources used during the run (defaults to 1 machine, 1 process).
        :param potcar_spec: dict. If the POTCAR file is not present anymore, this allows to pass a family and
            mapping to find the right POTCARs.
        :param settings: dict. Used for non-default parsing instructions, etc.
        """
        from aiida_vasp.calcs import immigrant as imgr
        remote_path = py_path.local(remote_path)
        proc_cls = imgr.VaspImmigrantJobProcess.build(cls)
        builder = proc_cls.get_builder()
        builder.code = code
        builder.options.max_wallclock_seconds = 1  # pylint: disable=no-member
        builder.options.resources = kwargs.get('resources', {'num_machines': 1, 'num_mpiprocs_per_machine': 1})  # pylint: disable=no-member
        settings = kwargs.get('settings', {})
        settings.update({'import_from_path': remote_path.strpath})
        builder.settings = get_data_node('parameter', dict=settings)
        with cmp_get_transport(code.get_computer()) as transport:
            with SandboxFolder() as sandbox:
                sandbox_path = py_path.local(sandbox.abspath)
                transport.get(remote_path.join('INCAR').strpath, sandbox_path.strpath)
                transport.get(remote_path.join('POSCAR').strpath, sandbox_path.strpath)
                transport.get(remote_path.join('POTCAR').strpath, sandbox_path.strpath, ignore_nonexisting=True)
                transport.get(remote_path.join('KPOINTS').strpath, sandbox_path.strpath)
                builder.parameters = imgr.get_incar_input(sandbox_path)
                builder.structure = imgr.get_poscar_input(sandbox_path)
                builder.potential = imgr.get_potcar_input(sandbox_path, potcar_spec=kwargs.get('potcar_spec', None))
                builder.kpoints = imgr.get_kpoints_input(sandbox_path)
                cls._immigrant_add_inputs(transport, remote_path=remote_path, sandbox_path=sandbox_path, builder=builder, **kwargs)
        return proc_cls, builder

    @classmethod
    def _immigrant_add_inputs(cls, transport, remote_path, sandbox_path, builder, **kwargs):
        pass
