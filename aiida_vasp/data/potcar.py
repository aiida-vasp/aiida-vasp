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
            pmg_potcar = PotcarData.get_or_create_from_file()
            |
            v
     _----------------------------------------------_
    ( exists for PotcarFileData with pmg_potcar.md5? )-----> no
     ^----------------------------------------------^         |
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
import tarfile
import tempfile
import shutil
from contextlib import contextmanager

import six
import py
from pymatgen.io.vasp import PotcarSingle
from aiida.backends.utils import get_automatic_user
from aiida.common import aiidalogger
from aiida.common.utils import md5_file, classproperty
from aiida.common.exceptions import UniquenessError, MultipleObjectsError, NotExistent
from aiida.orm import Group
from aiida.orm.data import Data
from aiida.orm.querybuilder import QueryBuilder

from aiida_vasp.data.archive import ArchiveData


@contextmanager
def temp_dir():
    try:
        tempdir = tempfile.mkdtemp()
        yield tempdir
    finally:
        shutil.rmtree(tempdir)


class PotcarMetadataMixin(object):
    """Provide common Potcar metadata access and querying functionality."""

    @classmethod
    def query_by_attrs(cls, **kwargs):
        """Find a Data node by attributes."""
        label = 'label'
        query_builder = cls.querybuild(label=label)
        filters = {}
        for attr_name, attr_val in kwargs.items():
            filters['attributes.{}'.format(attr_name)] = {'==': attr_val}
        query_builder.add_filter(label, filters)
        return query_builder

    @classmethod
    def find(cls, **kwargs):
        """Find a node by POTCAR metadata attributes given in kwargs."""
        query_builder = cls.query_by_attrs(**kwargs)
        return query_builder.one()[0]

    @classmethod
    def exists(cls, **kwargs):
        """Answers the question wether a node with attributes given in kwargs exists."""
        return bool(cls.query_by_attrs(**kwargs).count() >= 1)

    @property
    def md5(self):
        """Md5 hash of the POTCAR file (readonly)."""
        return self.get_attr('md5')

    @property
    def title(self):
        """Title of the POTCAR file (readonly)."""
        return self.get_attr('title')

    @property
    def functional(self):
        """Functional class of the POTCAR potential (readonly)."""
        return self.get_attr('functional')

    @property
    def element(self):
        """Chemical element described by the POTCAR (readonly)."""
        return self.get_attr('element')

    @property
    def symbol(self):
        """Element symbol property (VASP term) of the POTCAR potential (readonly)."""
        return self.get_attr('symbol')

    @property
    def original_file_name(self):
        """The name of the original file uploaded into AiiDA"""
        return self.get_attr('original_filename')

    def verify_unique(self):
        """Raise a UniquenessError if an equivalent node exists."""
        if self.exists(md5=self.md5):
            raise UniquenessError('A {} node already exists for this file.'.format(str(self.__class__)))

        other_attrs = self.get_attrs()
        other_attrs.pop('md5')
        if self.exists(**other_attrs):
            raise UniquenessError('A {} node with these attributes but a different file exists.'.format(str(self.__class__)))


class PotcarFileData(ArchiveData, PotcarMetadataMixin):
    """
    Store a POTCAR file in the db, never use as input to a calculation or workflow.

    .. warning:: Warning! Sharing nodes of this type may be illegal!

    In general POTCAR files may underly licence agreements, such as the ones distributed
    by the VASP group to VASP licence holders. Take care to not share such licenced data
    with non-licence holders.

    When writing a calculation plugin or workflow, do not use this as an input type, use :class:`aiida_vasp.data.potcar.PotcarData` instead!
    """

    _query_type_string = 'data.vasp.potcar_file.'
    _plugin_type_string = 'data.vasp.potcar_file.PotcarFileData'

    def set_file(self, filepath):
        """Initialize from a file path."""
        self.add_file(filepath)

    def add_file(self, src_abs, dst_filename=None):
        """Add the POTCAR file to the archive and set attributes."""
        if self._filelist:
            raise AttributeError('Can only hold one POTCAR file')
        super(PotcarFileData, self).add_file(src_abs, dst_filename)
        self._set_attr('md5', md5_file(src_abs))
        potcar = PotcarSingle.from_file(src_abs)
        self._set_attr('title', potcar.keywords['TITEL'])
        self._set_attr('functional', potcar.functional)
        self._set_attr('element', potcar.element)
        self._set_attr('symbol', potcar.symbol)
        self._set_attr('original_filename', src_abs)

    def store(self, with_transaction=True):
        """Ensure uniqueness and existence of a matching PotcarData node before storing."""
        _ = PotcarData.get_or_create(self)
        self.verify_unique()
        return super(PotcarFileData, self).store(with_transaction=with_transaction)

    def verify_unique(self):
        """Raise a UniquenessError if an equivalent node exists."""
        if self.exists(md5=self.md5):
            raise UniquenessError('A PotcarFileData already exists for this file.')

        other_attrs = self.get_attrs()
        other_attrs.pop('md5')
        if self.exists(**other_attrs):
            raise UniquenessError('A PotcarFileData with these attributes but a different file exists.')

    def get_file_obj(self):
        """Open a readonly file object to read the stored POTCAR file."""
        return self.archive.extractfile(self.archive.members[0])

    @classmethod
    def get_or_create(cls, filepath):
        """Get or create (store) a PotcarFileData node."""
        md5 = md5_file(filepath)
        if cls.exists(md5=md5):
            created = False
            node = cls.find(md5=md5)
        else:
            created = True
            node = cls(file=filepath)
            node.store()
        return node, created


class PotcarData(Data, PotcarMetadataMixin):
    """
    Store enough metadata about a POTCAR file to identify and find it.

    Meant to be used as an input to calculations. This node type holds no
    licenced data and can be freely shared without legal repercussions.
    """

    _meta_attrs = ['md5', 'title', 'functional', 'element', 'symbol', 'original_filename']
    _query_type_string = 'data.vasp.potcar.'
    _plugin_type_string = 'data.vasp.potcar.PotcarData'

    GROUP_TYPE = 'data.vasp.potcar.family'

    def set_potcar_file_node(self, potcar_file_node):
        """Initialize from a PotcarFileData node."""
        for attr_name in self._meta_attrs:
            self._set_attr(attr_name, potcar_file_node.get_attr(attr_name))

    def find_file_node(self):
        """Find and return the matching PotcarFileData node."""
        return PotcarFileData.find(**self.get_attrs())

    def store(self, with_transaction=True):
        """Ensure uniqueness before storing."""
        self.verify_unique()
        return super(PotcarData, self).store(with_transaction=with_transaction)

    @classmethod
    def get_or_create(cls, file_node):
        """Get or create (store) a PotcarData node."""
        if cls.exists(md5=file_node.md5):
            created = False
            node = cls.find(md5=file_node.md5)
        else:
            created = True
            node = cls(potcar_file_node=file_node)
            node.store()
        return node, created

    @classmethod
    def get_or_create_from_file(cls, file_path):
        """Get or create (store) a PotcarData node from a POTCAR file."""
        md5 = md5_file(file_path)
        file_node = PotcarFileData.find(md5=md5) if PotcarFileData.exists(md5=md5) else PotcarFileData(file=file_path)
        node, created = cls.get_or_create(file_node)
        if not file_node.is_stored:
            file_node.store()
        return node, created

    def get_family_names(self):
        """List potcar families to which this instance belongs."""
        return [group.name for group in Group.query(nodes=self, type_string=self.potcar_family_type_string)]

    @classproperty
    def potcar_family_type_string(cls):  # pylint: disable=no-self-argument
        return cls.GROUP_TYPE

    @classmethod
    def get_potcar_group(cls, group_name):
        """
        Return the PotcarFamily group with the given name.
        """
        return Group.get(name=group_name, type_string=cls.potcar_family_type_string)

    @classmethod
    def get_potcar_groups(cls, filter_elements=None, user=None):
        """
        List all names of groups of type PotcarFamily, possibly with some filters.

        :param filter_elements: A string or a list of strings.
               If present, returns only the groups that contains one POTCAR for
               every element present in the list. Default=None, meaning that
               all families are returned. A single element can be passed as a string.
        :param user: if None (default), return the groups for all users.
               If defined, it should be either a DbUser instance, or a string
               for the username (that is, the user email).
        """
        group_query_params = {"type_string": cls.potcar_family_type_string}

        if user is not None:
            group_query_params['user'] = user

        if isinstance(filter_elements, six.string_types):
            filter_elements = [filter_elements]

        if filter_elements is not None:
            normalized_elements_set = {element.capitalize() for element in filter_elements}

            group_query_params['node_attributes'] = {'element': normalized_elements_set}

        all_upf_groups = Group.query(**group_query_params)

        groups = [(group.name, group) for group in all_upf_groups]
        # Sort by name
        groups.sort()
        # Return the groups, without name
        return [item[1] for item in groups]

    @classmethod
    def get_potcars_dict(cls, structure, family_name):
        """
        Get a dictionary {kind: POTCAR} for all elements in a structure.

        :param structure: The structure to find POTCARs for
        :param family_name: The POTCAR family to be used
        """
        group_filters = {'name': {'==': family_name}, 'type': {'==': cls.potcar_family_type_string}}
        element_filters = {'attributes.element': {'in': [kind.symbol for kind in structure.kinds]}}
        query = QueryBuilder()
        query.append(Group, tag='family', filters=group_filters)
        query.append(cls, tag='potcar', member_of='family', filters=element_filters)

        result_potcars = {}
        for kind in structure.kinds:
            potcars_of_kind = [potcar[0] for potcar in query.all() if potcar[0].element == kind.symbol]
            if not potcars_of_kind:
                raise NotExistent('No POTCAR found for element {} in family {}'.format(kind.symbol, family_name))
            elif len(potcars_of_kind) > 1:
                raise MultipleObjectsError('More than one POTCAR for element {} found in family {}'.format(kind.symbol, family_name))
            result_potcars[kind] = potcars_of_kind[0]

        return result_potcars

    @classmethod
    def get_potcars_from_structure(cls, structure, family_name):
        """
        Given a POTCAR family name and a AiiDA
        structure, return a dictionary associating each kind name with its
        UpfData object.

        :raise MultipleObjectsError: if more than one UPF for the same element is
        found in the group.
        :raise NotExistent: if no UPF for an element in the group is
        found in the group.
        """
        return {kind.name: potcar for kind, potcar in cls.get_potcars_dict(structure, family_name).items()}

    @classmethod
    def upload_potcar_family(cls, folder, group_name, group_description=None, stop_if_existing=True):
        """
        Upload a set of POTCAR potentials as a family.

        :param folder: a path containing all POTCAR files to be added.
        :param group_name: the name of the group to create. If it exists and is
            non-empty, a UniquenessError is raised.
        :param group_description: a string to be set as the group description.
            Overwrites previous descriptions, if the group was existing.
        :param stop_if_existing: if True, check for the md5 of the files and,
            if the file already exists in the DB, raises a MultipleObjectsError.
            If False, simply adds the existing UPFData node to the group.
        """
        group, _ = Group.get_or_create(name=group_name, type_string=cls.GROUP_TYPE)

        if group.user != get_automatic_user():
            raise UniquenessError(
                'There is already a UpfFamily group with name {}, but it belongs to user {}, therefore you cannot modify it'.format(
                    group_name, group.user.email))

        if group_description:
            group.description = group_description

        potcars_found = cls.recursive_upload_potcar(folder, stop_if_existing=stop_if_existing)
        num_files = len(potcars_found)
        potcars_found = [(potcar, created, file_path) for potcar, created, file_path in potcars_found if potcar not in group.nodes]

        for potcar, created, file_path in potcars_found:
            if created:
                aiidalogger.debug('New PotcarData node %s created while uploading file %s for family %s', potcar.uuid, file_path,
                                  group_name)
            else:
                aiidalogger.debug('PotcarData node %s used instead of uploading file %s to family %s', potcar.uuid, file_path, group_name)

        group.add_nodes([potcar for potcar, created, file_path in potcars_found])

        num_uploaded = len(potcars_found)

        return num_files, num_uploaded

    @classmethod
    def recursive_upload_potcar(cls, folder, stop_if_existing=True):
        """Recursively search and upload POTCAR files in a folder."""
        list_created = []
        folder = py.path.local(folder)  # pylint: disable=no-name-in-module,no-member
        for subpath in folder.listdir():
            if subpath.isdir():
                list_created.extend(cls.recursive_upload_potcar(subpath, stop_if_existing))
            elif subpath.isfile() and 'POTCAR' in subpath.basename:
                if tarfile.is_tarfile(str(subpath)):
                    with temp_dir() as staging_dir:
                        with tarfile.TarFile(str(subpath)) as potcar_archive:
                            potcar_archive.extractall(staging_dir)
                        list_created.extend(cls.recursive_upload_potcar(staging_dir, stop_if_existing))
                else:
                    potcar, created = cls.get_or_create_from_file(str(subpath))
                    if stop_if_existing and not created:
                        raise ValueError('A POTCAR with identical MD5 to {} cannot be added with the stop_if_existing kwarg.'.format(
                            str(subpath)))
                    list_created.append((potcar, created, str(subpath)))

        return list_created
