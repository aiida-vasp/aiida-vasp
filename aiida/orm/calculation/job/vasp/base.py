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
    inp = {
        'valid_types': DataFactory(valid_types),
        'additional_parameter': additional_parameter,
        'linkname': linkname,
        'docstring': docstring
    }
    return inp

class CalcMeta(JobCalculation.__metaclass__):
    def __new__(cls, name, bases, classdict):
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

    def _super(self):
        return super(VaspCalcBase, self)

class TestVC(VaspCalcBase):
    test_input = input(valid_types = 'parameter', docstring = 'bla')
    default_parser = 'vasp.vasp'
