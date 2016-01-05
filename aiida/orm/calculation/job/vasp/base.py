from aiida.orm import JobCalculation, DataFactory, Node
from aiida.common.utils import classproperty

def make_use_methods(inputs):
    @classproperty
    def _use_methods(cls):
        retdict = JobCalculation._use_methods
        retdict.update(inputs)
        return retdict
    return _use_methods

def make_init_internal(cls, **kwargs):
    def _init_internal_params(self):
        self._super()._init_internal_params()
        for k, v in kwargs.iteritems():
            setattr(self, k, v)
    return _init_internal_params

def sequencify(argname):
    def wrapper(func):
        def new_func(*args, **kwargs):
            tlf = kwargs[argname]
            if isinstance(tlf, list) or isinstance(tlf, tuple):
                return func(*args, **kwargs)
            else:
                kwargs[argname] = [tlf]
                return func(*args, **kwargs)
        return new_func
    return wrapper

@sequencify('valid_types')
def input(valid_types = [], additional_parameter = None,
          linkname = None, docstring = ''):
    vt = []
    for t in valid_types:
        if isinstance(t, Node):
            vt.append(t)
        else:
            vt.append(DataFactory(t))
    inp = {
        'valid_types': vt,
        'additional_parameter': additional_parameter,
        'linkname': linkname,
        'docstring': docstring
    }
    return inp

def _super(self):
    return super(self.__class__, self)

class CalcMeta(JobCalculation.__metaclass__):
    def __new__(cls, name, bases, classdict):
        classdict['_super'] = _super
        inputs = {}
        internals = {}
        delete = []
        for k, v in classdict.iteritems():
            if isinstance(v, dict):
                if 'valid_types' in v:
                    if not v['linkname'] and not v['additional_parameter']:
                        v['linkname'] = k
                    elif not v['linkname'] and v['additional_parameter']:
                        v['linkname'] = classdict['_get_{}_linkname'.format(k)]
                    inputs.update({k: v})
                    delete.append(k)
            elif k == 'default_parser':
                internals['_default_parser'] = v
                delete.append(k)
            elif k == 'input_file_name':
                internals['_INPUT_FILE_NAME'] = v
                delete.append(k)
            elif k == 'output_file_name':
                internals['_OUTPUT_FILE_NAME'] = v
                delete.append(k)
        for k in delete:
            classdict.pop(k)
        classdict['_use_methods'] = make_use_methods(inputs)
        classdict['_init_internal_params'] = make_init_internal(cls, **internals)
        Calc = super(CalcMeta, cls).__new__(cls, name, bases, classdict)
        #setattr(Calc, '_init_internale_params', make_init_internal(**internals))
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
        poscar = tempfolder.get_abs_path('POSCAR')
        potcar = tempfolder.get_abs_path('POTCAR')
        kpoints = tempfolder.get_abs_path('KPOINTS')

        self.write_incar(inputdict, incar)
        self.write_poscar(inputdict, poscar)
        self.write_potcar(inputdict, potcar)
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
        calcinfo.codes_info = [codeinfo]

        return calcinfo

class TentativeVaspCalc(VaspCalcBase):
    '''vasp 3.5.3 calc with nice hopefully interface'''
    incar = input(valid_types='parameter', docstring='input parameters')
    poscar = input(valid_types='structure', docstring='aiida structure')
    potcar = input(valid_types='parameter', docstring='paw symbols or files')
    kpcar = input(valid_types='array.kpoints', docstring='kpoints array or mesh')
    default_parser = 'vasp.tentative'

    def write_incar(self, inputdict, dst):
        from incar import dict_to_incar
        with open(dst) as incar:
            incar.write(dict_to_incar(inputdict['incar']))

    def write_potcar(self, inputdict, dst):
        structure = inputdict['poscar']
        symbols = inputdict['potcar']
        sym = []
        # order the symbols according to order given in structure
        for kind in structure.get_kind_names():
            sym.append(symbols[kind])
        # find or create singlefile nodes for each symbol
        TODO()
        # concatenate
        TODO()
        ##
        ## Or actually it might be better to have a PAW node for each
        ## Kind and provide an external function to retrieve a PAW
        ## for a given symbol from the data bank
        ## That'd be back to additional paramter for potcar inp
        # then find all the linknames
        # and cat them in the right order
        ##

class TestVC(VaspCalcBase):
    test_input = input(valid_types='parameter', docstring='bla')
    default_parser = 'vasp.vasp'
