from aiida.orm import JobCalculation, DataFactory
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

def input(valid_types = [], additional_parameter = None,
          linkname = None, docstring = ''):
    if not isinstance(valid_types, str):
    inp = {
        'valid_types': DataFactory(valid_types),
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
        print name
        inputs = {}
        internals = {}
        delete = []
        for k, v in classdict.iteritems():
            if isinstance(v, dict):
                if 'valid_types' in v:
                    if not v['linkname'] and not v['additional_parameter']:
                        v['linkname'] = k
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
    potcar = input(valid_types=('parameter, singlefile')
    default_parser = 'vasp.tentative'

    def write_incar(self, inputdict, dst):
        from incar import dict_to_incar
        with open(dst) as incar:
            incar.write(dict_to_incar(inputdict['incar']))

    def write_poscar

class TestVC(VaspCalcBase):
    test_input = input(valid_types='parameter', docstring='bla')
    default_parser = 'vasp.vasp'
