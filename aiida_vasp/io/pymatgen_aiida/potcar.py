"""Find, import, compose and write POTCAR files."""
from functools import update_wrapper

from pymatgen.io.vasp import PotcarSingle
from aiida.common.utils import md5_file

from aiida_vasp.utils.aiida_utils import get_data_class


def delegate_method_kwargs(prefix='_init_with_'):
    """
    Get a kwargs delegating decorator.

    :params prefix: (str) common prefix of delegate functions
    """

    def decorator(meth):
        """Decorate a class method to delegate kwargs."""

        def wrapper(*args, **kwargs):
            for kwarg, value in kwargs.items():
                getattr(args[0], prefix + kwarg)(value)
            meth(*args, **kwargs)

        update_wrapper(wrapper, meth)
        return wrapper

    return decorator


class PotcarIo(object):
    """
    Use pymatgen.io.vasp.Potcar to deal with VASP pseudopotential IO.

    Instanciate with one of the following kwargs:

    :param pymatgen: a pymatgen.io.vasp.PotcarSingle instance
    :param path: (string) absolute path to the POTCAR file
    :param potcar_node: a PotcarData node
    :param potcar_file_node: a PotcarFileNode
    """

    def __init__(self, **kwargs):
        """Init from Potcar object or delegate to kwargs initializers."""
        self.potcar_obj = None
        self.md5 = None
        self.init_with_kwargs(**kwargs)

    @delegate_method_kwargs(prefix='_init_with_')
    def init_with_kwargs(self, **kwargs):
        """Delegate initialization to _init_with - methods."""

    def _init_with_path(self, filepath):
        self.potcar_obj = PotcarSingle(filepath)
        self.md5 = md5_file(filepath)
        get_data_class('vasp.potcar').get_or_create(file=filepath)

    def _init_with_potcar_file_node(self, node):
        with node.get_file_obj() as potcar_fo:
            self.potcar_obj = PotcarSingle(potcar_fo.read())
        self.md5 = node.md5

    def _init_with_potcar_node(self, node):
        self._init_with_potcar_file_node(node.find_file_node())

    @property
    def pymatgen(self):
        return self.potcar_obj

    @property
    def file_node(self):
        return get_data_class('vasp.potcar').find(md5=self.md5).find_file_node()

    @property
    def node(self):
        return get_data_class('vasp.potcar').find(md5=self.md5)
