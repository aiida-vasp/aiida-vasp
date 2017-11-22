# pylint: disable=abstract-method
"""
Attempt to create a convenient but licence-respecting storage system that also guarantees provenience.

Consists of two classes, PotcarData and PotcarFileData. Between the two data node classes exists a
one to one mapping but never a DbLink of any kind. The mapping must be defined in terms of a POTCAR
file hash sum.

Reasons for not using a file system based solution in general:

    * simplicity -> no necessity to define an fs based storage / retrieval schema
    * storage schema can be updated without manual user interaction
    * with fs based it is possible to lose enhanced provenance locally by deleting a file
    * This is easier to share between machines for same user / group members

Reasons for not using fs paths:

    * migrating to a new machine involves reinstating file hierarchy, might be non-trivial
    * corner cases with links, recursion etc

Reasons for not using pymatgen system:

    * changing an environment variable would invalidate provenance / disable reusing potentials
    * would block upgrading to newer pymatgen versions if they decide to change


Note::

    An fs based solution can be advantageous but should be 'expert mode' and not
    default solution due to provenance tradeoffs.

The following requirements have to be met:

    * The file hash attribute of PotcarFileData is unique in the Db
    * The file hash attribute of PotcarData is unique in the Db
    * Both classes can easily and quickly be found via the hash attribute
    * A PotcarData node can be exported without exporting the PotcarFileData node
    * The corresponding PotcarData node can be created at any time from the PotcarFileData node
    * PotcarFileData nodes are expected to be grouped in DbGroups called 'families'
    * The PotcarFileData nodes can be found according to their 'functional type' (pymatgen term)

The following would be nice to also allow optionally:

    * To pre-upload the files to a remote computer from a db and concat them right on there (to save traffic)
    * To use files directly on the remote computer (disclaimer: will never be as secure / tested)
    * To use existing pymatgen-style potentials library (disclaimer: support might break)

It is not to be expected for hundreds of distinct Potcar families to be present in the same database.

The mechanism for reading a POTCAR file into the Db::

    +-----------------------+
    [ parsing a POTCAR file ]
    +-----------------------+
            |
            v
            pmg_potcar = Potcar.from_file()
            |
            v
     _-----------------------------------------------_
    ( exists for PotcarFileData with pmg_potcar.hash? )-----> no
     ^-----------------------------------------------^        |
            |                                                 v
            v                                                 create
            yes                                               |
            |                                                 |
            v                                                 v
     _-------------------------_                             _-------------------------_
    ( Family given to parse to? ) -------------> no -+      ( Family given to parse to? )
     ^-------------------------^                     |       ^-------------------------^
            |                                        |        |         |
            v                                        |        |         no
            yes<------------------------------------]|[-------+         |
            |                                        |                  choose family according to functional type (with fallback?)
            v                                        |                  |
            add existing PotcarFileData to family<--]|[-----------------+
            |                                        |
            |                     +------------------+
            v                     v
     _--------------------------------_
    ( exists corresponding PotcarData? )-----> no -----> create
     ^--------------------------------^ <------------------+
            |
            v
            return corresponding PotcarData

The mechanism for writing one or more PotcarData to file (from a calculation)::

    +-----------------------+
    [ Writing a POTCAR file ]
    +-----------------------+
            |
            v
            for each PotcarData node:
                get corresponding PotcarFileData <-> query for same symbol, family, hash, do not use links
            |
            v
            for each PotcarFileData:
                create a pymatgen PotcarSingle object
            |
            v
            create a pymatgen Potcar object from all the PotcarSingle objects
            (maybe need to take care to order same as in POSCAR)
            |
            v
            use Potcar.write_file

"""
from pymatgen.io.vasp import PotcarSingle
from aiida.common.utils import md5_file
from aiida.orm.querybuilder import QueryBuilder
from aiida.orm.data import Data

from aiida_vasp.data.archive import ArchiveData


class PotcarFileData(ArchiveData):
    """Store a POTCAR file in the db."""

    def set_file(self, filepath):
        self.add_file(filepath)

    def add_file(self, src_abs, dst_filename=None):
        """Add the POTCAR file to the archive and set attributes."""
        if self._filelist:
            raise AttributeError('Can only hold one POTCAR file')
        super(PotcarFileData, self).add_file(src_abs, dst_filename)
        self.set_attr('md5', md5_file(src_abs))
        potcar = PotcarSingle.from_file(src_abs)
        self.set_attr('title', potcar.keywords['TITEL'])
        self.set_attr('functional', potcar.functional)
        self.set_attr('element', potcar.element)
        self.set_attr('symbol', potcar.symbol)

    def get_file_obj(self):
        return self.archive.extractfile(self.archive.members[0])

    @property
    def md5(self):
        return self.get_attr('md5')

    @property
    def title(self):
        return self.get_attr('title')

    @property
    def functional(self):
        return self.get_attr('functional')

    @property
    def element(self):
        return self.get_attr('element')

    @property
    def symbol(self):
        return self.get_attr('symbol')


class PotcarData(Data):
    """Store enough metadata about a POTCAR file to identify it."""
    _meta_attrs = ['md5', 'title', 'functional', 'element', 'symbol']

    def set_potcar_file_node(self, potcar_file_node):
        for attr_name in self._meta_attrs:
            self.set_attr(attr_name, potcar_file_node.get_attr(attr_name))

    def find_file_node(self):
        """Find and return the matching PotcarFileData node."""
        query_builder = QueryBuilder()
        query_builder.append(PotcarFileData, tag='potcar_file')
        filters = {}
        for attr_name in self._meta_attrs:
            filters['attributes.{}'.format(attr_name)] = {'==': self.get_attr(attr_name)}
        query_builder.add_filter('potcar_file', filters)
        return query_builder.one()

    @property
    def md5(self):
        return self.get_attr('md5')

    @property
    def title(self):
        return self.get_attr('title')

    @property
    def functional(self):
        return self.get_attr('functional')

    @property
    def element(self):
        return self.get_attr('element')

    @property
    def symbol(self):
        return self.get_attr('symbol')
