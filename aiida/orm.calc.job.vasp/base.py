from aiida.orm import JobCalculation
from aiida.common.utils import classproperty
from aiida.common.datastructures import CalcInfo, CodeInfo


def make_use_methods(inputs):
    @classproperty
    def _use_methods(cls):
        retdict = JobCalculation._use_methods
        for inp, dct in inputs.iteritems():
            ln = dct['linkname']
            if isinstance(ln, classmethod):
                dct['linkname'] = getattr(cls, ln.__func__.__name__)
        retdict.update(inputs)
        return retdict
    return _use_methods


def make_init_internal(cls, **kwargs):
    def _update_internal_params(self):
        for k, v in kwargs.iteritems():
            setattr(self, k, v)
    return _update_internal_params


def seqify(seq_arg):
    if not isinstance(seq_arg, list) and not isinstance(seq_arg, tuple):
        seq_arg = (seq_arg)
    return seq_arg


class Input(object):
    '''
    Utility class that handles mapping between CalcMeta's and Calculation's
    interfaces for defining input nodes
    '''
    def __init__(self, types, param=None, ln=None, doc=''):
        self.types = seqify(types)
        self.param = param
        self._ln = ln
        self.doc = ''

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
    Metaclass responsible to simplify routine Calculation setup
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
        for k in delete:
            classdict.pop(k)
        classdict['_use_methods'] = make_use_methods(inputs)
        classdict['_update_internal_params'] = make_init_internal(
            cls, **internals)
        Calc = super(CalcMeta, cls).__new__(cls, name, bases, classdict)
        return Calc


class VaspCalcBase(JobCalculation):
    __metaclass__ = CalcMeta
    input_file_name = 'INCAR'
    output_file_name = 'OUTCAR'

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
        self.write_structure(inputdict, structure)
        self.write_paw(inputdict, paw)
        self.write_kpoints(inputdict, kpoints)

        # calcinfo
        calcinfo = CalcInfo()
        calcinfo.uuid = self.uuid
        calcinfo.retrieve_list = [
            'CHG',
            'CHGCAR',
            'CONTCAR',
            'DOSCAR',
            'EIGENVAL',
            'OSZICAR',
            'OUTCAR',
            'PCDAT',
            'PROCAR',
            'WAVECAR',
            'XDATCAR',
            'vasprun.xml'
        ]
        codeinfo = CodeInfo()
        codeinfo.code_uuid = self.get_code().uuid
        codeinfo.code_pk = self.get_code().pk
        calcinfo.codes_info = [codeinfo]

        return calcinfo

    def _init_internal_params(self):
        super(VaspCalcBase, self)._init_internal_params()
        self._update_internal_params()

    def verify_inputs(self):
        pass


# only implements gamma centered monkhorst
# no expert mode with manual basis vectors
# no fully automatic mode
# no original monkhorst (equivalent to 0.5 shift though)
kpmtemp = '''\
Automatic mesh
0
Gamma
{N[0]} {N[1]} {N[2]}
{s[0]} {s[1]} {s[2]}\
'''

# only implements direct explicit list with weights
# no cartesian
# no tetrahedron method
# no bandstructure lines
kpltemp = '''\
Explicit list
{N}
Direct
{klist}\
'''
kplitemp = '''\
k[0] k[1] k[2] w\
'''


class TentativeVaspCalc(VaspCalcBase):
    '''vasp 3.5.3 calc with nice hopefully interface'''
    incar = Input(types='parameter', doc='input parameters')
    structure = Input(types='structure', doc='aiida structure')
    paw = Input(types='vasp.potpaw', doc='paw symbols or files', param='kind')
    kpoints = Input(types='array.kpoints', doc='kpoints array or mesh')
    default_parser = 'vasp.tentative'

    def _init_internal_params(self):
        super(TentativeVaspCalc, self)._init_internal_params()
        self._update_internal_params()

    def verify_inputs(self, inputdict):
        incar = inputdict['incar'].get_dict()
        keys = map(lambda k: k.lower(), incar.keys())
        need_kp = True
        if 'kspacing' in keys and 'kgamma' in keys:
            need_kp = False
        self.use_kp = True
        kp = inputdict.get('kpoints')
        if not need_kp:
            msg = 'INCAR contains KSPACING and KGAMMA: '
            if not kp:
                msg += 'KPOINTS omitted'
                self.use_kp = False
            else:
                msg += 'KPOINTS still used'
            self.logger.info(msg)

    def write_incar(self, inputdict, dst):
        from incar import dict_to_incar
        with open(dst, 'w') as incar:
            incar.write(dict_to_incar(inputdict['incar'].get_dict()))

    @classmethod
    def _get_paw_linkname(cls, kind):
        return 'paw_{}'.format(kind)

    def write_paw(self, inputdict, dst):
        import subprocess32 as sp
        structure = inputdict['structure']
        catcom = ['cat']
        # order the symbols according to order given in structure
        for kind in structure.get_kind_names():
            paw = inputdict[self._get_paw_linkname(kind)]
            catcom.append(paw.get_abs_path('POTCAR'))
        # cat the pawdata nodes into the file
        with open(dst, 'w') as pc:
            sp.check_call(catcom, stdout=pc)

    def write_structure(self, inputdict, dst):
        from ase.io.vasp import write_vasp
        structure = inputdict['structure']
        self.elements = structure.get_kind_names()
        with open(dst, 'w') as structure:
            write_vasp(structure, structure.get_ase(), vasp5=True)

    def write_kpoints(self, inputdict, dst):
        if self.use_kp:
            kp = inputdict['kpoints']
            try:
                mesh, offset = kp.get_kpoints_mesh()
                with open(dst, 'w') as kpoints:
                    kps = kpmtemp.format(N=mesh, s=offset)
                    kpoints.write(kps)
            except AttributeError:
                kpl, weights = kp.get_kpoints(also_weights=True)
                kw = zip(kpl, weights)
                with open(dst, 'w') as kpoints:
                    kpls = '\n'.join([kplitemp.format(k[0], k[1]) for k in kw])
                    kps = kpltemp.format(N=len(kw), klist=kpls)
                    kpoints.write(kps)

    @property
    def elements(self):
        return self.get_attr('elements')

    @elements.setter
    def elements(self, value):
        self._set_attr('elements', value)


class TestVC(VaspCalcBase):
    test_input = Input(types='parameter', doc='bla')
    default_parser = 'vasp.vasp'
