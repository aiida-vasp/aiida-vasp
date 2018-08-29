"""Find, import, compose and write POTCAR files."""
from functools import update_wrapper
from itertools import groupby
import re
import six

from py import path as py_path  # pylint: disable=no-name-in-module,no-member

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

    def _init_with_path(self, file_path):
        node, _ = get_data_class('vasp.potcar').get_or_create_from_file(file_path=file_path)
        self.md5 = node.md5

    def _init_with_potcar_file_node(self, node):
        self.md5 = node.md5

    def _init_with_potcar_node(self, node):
        self._init_with_potcar_file_node(node.find_file_node())

    def _init_with_contents(self, contents):
        node, _ = get_data_class('vasp.potcar').get_or_create_from_contents(contents)
        self.md5 = node.md5

    @property
    def pymatgen(self):
        if not self.potcar_obj:
            self.potcar_obj = self.file_node.get_pymatgen()
        return self.potcar_obj

    @property
    def file_node(self):
        return get_data_class('vasp.potcar').find_one(md5=self.md5).find_file_node()

    @property
    def node(self):
        return get_data_class('vasp.potcar').find_one(md5=self.md5)

    @property
    def content(self):
        return self.file_node.get_content()

    @classmethod
    def from_(cls, potcar):
        """Determine the best guess at how the input represents a POTCAR file and construct a PotcarIo instance based on that."""
        if isinstance(potcar, (six.string_types)):
            if py_path.local(potcar).exists():
                potcar = cls(path=potcar)
            else:
                potcar = cls(contents=potcar)
        elif isinstance(potcar, get_data_class('vasp.potcar')):
            potcar = cls(potcar_node=potcar)
        elif isinstance(potcar, get_data_class('vasp.potcar_file')):
            potcar = cls(potcar_file_node=potcar)
        elif isinstance(potcar, PotcarIo):
            pass
        else:
            potcar = cls.from_(str(potcar))
        return potcar

    def __eq__(self, other):
        return self.md5 == other.md5


class MultiPotcarIo(object):
    """Handle file i/o for POTCAR files with one or more potentials."""

    def __init__(self, potcars=None):
        self._potcars = []
        if potcars:
            for potcar in potcars:
                self.append(PotcarIo.from_(potcar))

    def append(self, potcar):
        self._potcars.append(PotcarIo.from_(potcar))

    def write(self, path):
        path = py_path.local(path)
        with path.open('w') as dest_fo:
            for potcar in self._potcars:
                dest_fo.write(potcar.content)

    @classmethod
    def read(cls, path):
        """Read a POTCAR file that may contain one or more potentials into a list of PotcarIo objects."""
        potcars = cls()
        path = py_path.local(path)
        with path.open('r') as potcar_fo:
            potcar_strings = re.compile(r"\n?(\s*.*?End of Dataset\n)", re.S).findall(potcar_fo.read())

        for potcar_contents in potcar_strings:
            potcars.append(PotcarIo.from_(potcar_contents))
        return potcars

    @property
    def potcars(self):
        return self._potcars

    @classmethod
    def from_structure(cls, structure, potentials_dict):
        """Create a MultiPotcarIo from an AiiDA `StructureData` object and a dictionary with a potential for each kind in the structure."""
        symbol_order = cls.potentials_order(structure)
        return cls(potcars=[potentials_dict[symbol] for symbol in symbol_order])

    def get_potentials_dict(self, structure):
        """
        Get a dictionary {kind_name: PotcarData} that would fit the structure.

        If the PotcarData contained in MultiPotcarIo do not match the structure, an exception is raised.
        """
        structure_elements = structure.get_symbols_set()
        if structure_elements != self.element_symbols:
            raise ValueError('structure elements do not match POTCAR elements')
        if len(structure.get_kind_names()) != len(structure_elements):
            raise ValueError('structure has more kind names than elements')

        element_potcars = {potcario.node.element: potcario.node for potcario in self.potcars}
        return {kind.name: element_potcars[kind.symbol] for kind in structure.kinds}

    @property
    def element_symbols(self):
        return {potcario.node.element for potcario in self.potcars}

    @classmethod
    def potentials_order(cls, structure):
        return [kind[0] for kind in cls.count_kinds(structure)]

    @classmethod
    def count_kinds(cls, structure):
        """
        Count consecutive kinds that compose the different sites.

        :return: [(kind_name, num), ... ]
        """
        kind_name_order = [site.kind_name for site in structure.sites]
        groups = groupby(kind_name_order)
        counts = [(label, sum(1 for _ in group)) for label, group in groups]
        return counts

    @property
    def max_enmax(self):
        return max([potcario.pymatgen.enmax for potcario in self.potcars])
