"""
POTCAR parser.

--------------
The file parser that handles the parsing of POTCAR files. Also contains methods to
find, import, compose and write POTCAR files.
"""

import os
import re
from itertools import groupby
from pathlib import Path

from parsevasp.potcar import Potcar

from aiida_vasp.parsers.content_parsers.base import BaseFileParser
from aiida_vasp.utils.aiida_utils import get_data_class
from aiida_vasp.utils.delegates import delegate_method_kwargs


class PotcarParser(BaseFileParser):
    """A lightweight interface that provides access to POTCAR metadata parsing.

    Similar to the other content parser for VASP in structure, but only used directly in the POTCAR
    handling logic.

    """

    DEFAULT_SETTINGS = {}

    PARSABLE_QUANTITIES = {}

    def _init_from_handler(self, handler):
        """Initialize using a file like handler.

        Parameters
        ----------
        handler : object
            A file like object that provides the necessary content to be parsed.

        """

        try:
            self._content_parser = Potcar(file_handler=handler, logger=self._logger)
        except SystemExit:
            self._logger.warning('Parsevasp exited abnormally.')

    @property
    def metadata(self):
        """Return the metadata Potcar instance."""
        return self._content_parser

    def _init_from_data(self, data):
        """No need to init from an AiiDA data structure."""
        raise NotImplementedError('PotcarParser does not implement a _init_from_data() method.')

    def _content_data_to_content_parser(self):
        """Since no need to accept AiiDA data structure, no need to convert it."""
        raise NotImplementedError('PotcarParser does not implement a _content_data_to_content_parser() method.')


class PotcarIo:  # pylint: disable=useless-object-inheritance
    """
    Deals with VASP input output of POTCAR files.

    Instanciate with one of the following kwargs:

    :param path: (string) absolute path to the POTCAR file
    :param potcar_node: a PotcarData node
    :param potcar_file_node: a PotcarFileNode
    :param contents: a string with the POTCAR content
    """

    def __init__(self, **kwargs):
        """Init from Potcar object or delegate to kwargs initializers."""
        self.potcar_obj = None
        self.sha512 = None
        self.init_with_kwargs(**kwargs)

    @delegate_method_kwargs(prefix='_init_with_')
    def init_with_kwargs(self, **kwargs):
        """Delegate initialization to _init_with - methods."""

    def _init_with_path(self, file_path):
        """Initialize with a path."""
        node, _ = get_data_class('vasp.potcar').get_or_create_from_file(file_path=file_path)
        self.sha512 = node.sha512

    def _init_with_potcar_file_node(self, node):
        """Initialize with an existing potential file node."""
        self.sha512 = node.sha512

    def _init_with_potcar_node(self, node):
        """Initialize with an existing potential node."""
        self._init_with_potcar_file_node(node.find_file_node())

    def _init_with_contents(self, contents):
        """Initialize with a string."""
        try:
            contents = contents.encode('utf-8')
        except AttributeError:
            pass
        node, _ = get_data_class('vasp.potcar').get_or_create_from_contents(contents)
        self.sha512 = node.sha512

    @property
    def file_node(self):
        return get_data_class('vasp.potcar').find_one(sha512=self.sha512).find_file_node()

    @property
    def node(self):
        return get_data_class('vasp.potcar').find_one(sha512=self.sha512)

    @property
    def content(self):
        return self.file_node.get_content()

    @classmethod
    def from_(cls, potcar):
        """Determine the best guess at how the input represents a POTCAR file and construct
        a PotcarIo instance based on that."""
        if isinstance(potcar, str):
            try:
                path_exists = Path(potcar).exists()
            except OSError:
                # We failed possibly due to a too long filename or that the potcar content is in fact
                # potcar content, revert to the os module to check if it exists
                path_exists = os.path.exists(potcar)
            if path_exists:
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
            potcar = cls.from_(potcar)
        return potcar

    def __eq__(self, other):
        return self.sha512 == other.sha512


class MultiPotcarIo:  # pylint: disable=useless-object-inheritance
    """Handle file i/o for POTCAR files with one or more potentials."""

    def __init__(self, potcars=None):
        self._potcars = []
        if potcars:
            for potcar in potcars:
                self.append(PotcarIo.from_(potcar))

    def append(self, potcar):
        self._potcars.append(PotcarIo.from_(potcar))

    def write(self, path):
        path = Path(path)
        with path.open('wb') as dest_fo:
            for potcar in self._potcars:
                dest_fo.write(potcar.content)

    @classmethod
    def read(cls, path):
        """Read a POTCAR file that may contain one or more potentials into a list of PotcarIo objects."""
        potcars = cls()
        path = Path(path)
        with path.open('r', encoding='utf8') as potcar_fo:
            potcar_strings = re.compile(r'\n?(\s*.*?End of Dataset\n)', re.S).findall(potcar_fo.read())

        for potcar_contents in potcar_strings:
            potcars.append(PotcarIo.from_(potcar_contents))
        return potcars

    @property
    def potcars(self):
        return self._potcars

    @classmethod
    def from_structure(cls, structure, potentials_dict):
        """Create a MultiPotcarIo from an AiiDA `StructureData` object and a dictionary with a
        potential for each kind in the structure."""
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
