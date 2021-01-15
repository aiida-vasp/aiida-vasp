# pylint: disable=abstract-method
"""
Representation of the POTCAR files.

-----------------------------------
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
    ( exists for PotcarFileData with pmg_potcar.sha512? )-----> no
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
# pylint: disable=import-outside-toplevel
from __future__ import print_function

import re
import os
import hashlib
import tarfile
import tempfile
import shutil
from contextlib import contextmanager
from collections import namedtuple

from pathlib import Path
from pymatgen.io.vasp import PotcarSingle
from aiida.common import AIIDA_LOGGER as aiidalogger
from aiida.common.utils import classproperty
from aiida.common.exceptions import UniquenessError, NotExistent
from aiida.orm import Group
from aiida.orm import Data
from aiida.orm import QueryBuilder

from aiida_vasp.data.archive import ArchiveData
from aiida_vasp.utils.aiida_utils import get_current_user, querybuild
from aiida_vasp.utils.delegates import delegate_method_kwargs


def normalize_potcar_contents(potcar_contents):
    """Normalize whitespace in a POTCAR given as a string."""
    try:
        potcar_contents = potcar_contents.decode()
    except AttributeError:
        pass
    normalized = re.sub(r'[ \t]+', r' ', potcar_contents)  # multiple spaces
    normalized = re.sub(r'[\n\r]\s*', r'\n', normalized)  # line breaks and spaces afterwards / empty lines
    normalized = re.sub(r'^\s*', r'', normalized)  # spaces / empty lines at the very beginning
    normalized = re.sub(r'\s*$', r'\n', normalized)  # trailing endline
    return normalized


def sha512_potcar(potcar_contents):
    """Hash the contents of a POTCAR file (given as str)."""
    sha512_hash = hashlib.sha512()
    sha512_hash.update(normalize_potcar_contents(potcar_contents).encode('utf-8'))
    return sha512_hash.hexdigest()


@contextmanager
def temp_dir():
    """Temporary directory context manager that deletes the tempdir after use."""
    try:
        tempdir = tempfile.mkdtemp()
        yield Path(tempdir)
    finally:
        shutil.rmtree(tempdir)


@contextmanager
def temp_potcar(contents):
    """Temporary POTCAR file from contents."""
    with temp_dir() as tempdir:
        potcar_file = tempdir / 'POTCAR'
        with potcar_file.open('wb') as potcar_fo:
            potcar_fo.write(contents)
        yield potcar_file


def extract_tarfile(file_path):
    """Extract a .tar archive into an appropriately named folder, return the path of the folder, avoid extracting if folder exists."""
    new_path = None
    with tarfile.open(str(file_path)) as archive:
        new_dir = file_path.name.split('.tar')[0]
        new_path = file_path.parent / new_dir
        if not new_path.exists():
            archive.extractall(str(new_path))

    return new_path


def by_older(left, right):
    if left.ctime < right.ctime:
        return -1
    if left.ctime > right.ctime:
        return 1
    return 0


def by_user(left, right):
    if left.user.is_active and not right.user.is_active:
        return -1
    if not left.user.is_active and right.user.is_active:
        return 1
    return 0


class PotcarWalker(object):  # pylint: disable=useless-object-inheritance
    """
    Walk the file system and find POTCAR files under a given directory.

    Build a list of potcars including their full path and wether they are archived inside a tar archive.
    """

    def __init__(self, path):
        # Only accept a Path object or a string
        if isinstance(path, Path):
            self.path = path
        elif isinstance(path, str):
            self.path = Path(path)
        else:
            raise ValueError('The supplied path is not a Path object or a string.')
        self.potcars = set()

    def walk(self):
        """Walk the folder tree to find POTCAR, extracting any tar archives along the way."""
        if self.path.is_file():
            extracted = self.file_dispatch(self.path.parent, [], self.path.name)
            if extracted:
                self.path = extracted
                self.walk()
        else:
            for root, dirs, files in os.walk(str(self.path)):
                for file_name in files:
                    self.file_dispatch(root, dirs, file_name)

    def file_dispatch(self, root, dirs, file_name):
        """Add POTCAR files to the list and dispatch handling of different kinds of files to other methods."""
        file_path = Path(root) / file_name
        if tarfile.is_tarfile(str(file_path)):
            return self.handle_tarfile(dirs, file_path)
        if 'POTCAR' in file_name:
            self.potcars.add(file_path)
        return None

    @classmethod
    def handle_tarfile(cls, dirs, file_path):
        """Handle .tar archives: extract and add the extracted folder to be searched."""
        new_dir = extract_tarfile(file_path)
        if new_dir not in dirs:
            dirs.append(str(new_dir))
        return new_dir


class PotcarMetadataMixin(object):  # pylint: disable=useless-object-inheritance
    """Provide common Potcar metadata access and querying functionality."""
    _query_label = 'label'

    @classmethod
    def query_by_attrs(cls, query=None, **kwargs):
        """Find a Data node by attributes."""
        label = cls._query_label
        if not query:
            query = querybuild(cls, tag=label)
        filters = {}
        for attr_name, attr_val in kwargs.items():
            filters['attributes.{}'.format(attr_name)] = {'==': attr_val}
        if cls._HAS_MODEL_VERSIONING:
            filters['attributes._MODEL_VERSION'] = {'==': kwargs.get('model_version', cls._VERSION)}
        query.add_filter(label, filters)
        return query

    @classmethod
    def find(cls, **kwargs):
        """Find nodes by POTCAR metadata attributes given in kwargs."""
        query_builder = cls.query_by_attrs(**kwargs)
        if not query_builder.count():
            raise NotExistent()
        results = [result[0] for result in query_builder.all()]
        from functools import cmp_to_key
        results.sort(key=cmp_to_key(by_older))
        return results

    @classmethod
    def find_one(cls, **kwargs):
        """
        Find one single node.

        Raise an exception if there are multiple.
        """
        res = cls.find(**kwargs)
        if len(res) > 1:
            if not all([True for node in res if node.sha512 == res[0].sha512]):
                raise UniquenessError('Multiple nodes found satisfying {}'.format(kwargs))
        return res[0]

    @classmethod
    def exists(cls, **kwargs):
        """Answers the question wether a node with attributes given in kwargs exists."""
        return bool(cls.query_by_attrs(**kwargs).count() >= 1)

    @property
    def sha512(self):
        """Sha512 hash of the POTCAR file (readonly)."""
        return self.get_attribute('sha512')

    @property
    def title(self):
        """Title of the POTCAR file (readonly)."""
        return self.get_attribute('title')

    @property
    def functional(self):
        """Functional class of the POTCAR potential (readonly)."""
        return self.get_attribute('functional')

    @property
    def element(self):
        """Chemical element described by the POTCAR (readonly)."""
        return self.get_attribute('element')

    @property
    def symbol(self):
        """Element symbol property (VASP term) of the POTCAR potential (readonly)."""
        return self.get_attribute('symbol')

    @property
    def original_file_name(self):
        """The name of the original file uploaded into AiiDA."""
        return self.get_attribute('original_filename')

    @property
    def full_name(self):
        """The name of the original file uploaded into AiiDA."""
        return self.get_attribute('full_name')

    @property
    def potential_set(self):
        """The name of the original file uploaded into AiiDA."""
        return self.get_attribute('potential_set')

    def verify_unique(self):
        """Raise a UniquenessError if an equivalent node exists."""
        from copy import deepcopy
        if self.exists(sha512=self.sha512):
            raise UniquenessError('A {} node already exists for this file.'.format(str(self.__class__)))

        other_attrs = deepcopy(self.attributes)

        other_attrs.pop('sha512')
        if self.exists(**other_attrs):
            raise UniquenessError('A {} node with these attributes but a different file exists:\n{}'.format(
                str(self.__class__), str(other_attrs)))


class VersioningMixin(object):  # pylint: disable=useless-object-inheritance
    """Minimalistic Node versioning."""
    _HAS_MODEL_VERSIONING = True
    _VERSION = None

    def set_version(self):
        self.set_attribute('_MODEL_VERSION', self._VERSION)

    @property
    def model_version(self):
        return self.get_attribute('_MODEL_VERSION')

    @classmethod
    def old_versions_in_db(cls):
        """Determine whether there are Nodes created with an older version of the model."""
        label = 'versioned'
        query = querybuild(cls, tag=label)
        filters = {'attributes._MODEL_VERSION': {'<': cls._VERSION}}
        query.add_filter(label, filters)
        return bool(query.count() >= 1)


class PotcarFileData(ArchiveData, PotcarMetadataMixin, VersioningMixin):
    """
    Store a POTCAR file in the db, never use as input to a calculation or workchain.

    .. warning:: Warning! Sharing nodes of this type may be illegal!

    In general POTCAR files may underly licence agreements, such as the ones distributed
    by the VASP group to VASP licence holders. Take care to not share such licenced data
    with non-licence holders.

    When writing a calculation plugin or workchain, do not use this as an input type,
    use :class:`aiida_vasp.data.potcar.PotcarData` instead!
    """

    _query_label = 'potcar_file'
    _query_type_string = 'data.vasp.potcar_file.'
    _plugin_type_string = 'data.vasp.potcar_file.PotcarFileData.'
    _VERSION = 1

    def __init__(self, *args, **kwargs):
        # remove file in kwargs as this is not accepted in the subsequent inits
        path = kwargs.pop('file', None)
        super(PotcarFileData, self).__init__(*args, **kwargs)
        if path is not None:
            # Only allow a Path object or a string
            if isinstance(path, Path):
                self.init_with_kwargs(file=path)
            elif isinstance(path, str):
                self.init_with_kwargs(file=Path(path))
            else:
                raise ValueError('The supplied argument for file is not a Path object or a string.')

    @delegate_method_kwargs(prefix='_init_with_')
    def init_with_kwargs(self, **kwargs):
        """Delegate initialization to _init_with - methods."""

    def _init_with_file(self, filepath):
        """Initiqalize from a file path."""
        self.add_file(filepath)

    def add_file(self, src_abs, dst_filename=None):
        """Add the POTCAR file to the archive and set attributes."""
        self.set_version()
        if self._filelist:
            raise AttributeError('Can only hold one POTCAR file')
        super(PotcarFileData, self).add_file(src_abs, dst_filename)
        self.set_attribute('sha512', self.get_file_sha512(src_abs))
        # PotcarSingle needs a string for path
        potcar = PotcarSingle.from_file(str(src_abs))
        self.set_attribute('title', potcar.keywords['TITEL'])
        self.set_attribute('functional', potcar.functional)
        self.set_attribute('element', potcar.element)
        self.set_attribute('symbol', potcar.symbol)
        src_path = src_abs.resolve()
        src_rel = src_path.relative_to(src_path.parents[2])  # familyfolder/Element/POTCAR
        # Make sure we store string elements of Path in the attributes
        self.set_attribute('original_filename', str(src_rel))
        dir_name = src_path.parent
        dir_name = dir_name.name
        self.set_attribute('full_name', str(dir_name))
        self.set_attribute('potential_set', str(src_path.parts[-3]))

    @classmethod
    def get_file_sha512(cls, path):
        """Get the sha512 sum for a POTCAR file (after whitespace normalization)."""
        path = Path(path)
        with path.open('r') as potcar_fo:
            sha512 = sha512_potcar(potcar_fo.read())
        return sha512

    @classmethod
    def get_contents_sha512(cls, contents):
        """Get the sha512 sum for the contents of a POTCAR file (after normalization)."""
        return sha512_potcar(contents)

    # pylint: disable=arguments-differ
    def store(self, *args, **kwargs):
        """Ensure uniqueness and existence of a matching PotcarData node before storing."""
        self.set_version()
        _ = PotcarData.get_or_create(self)
        self.verify_unique()
        return super(PotcarFileData, self).store(*args, **kwargs)

    @contextmanager
    def get_file_obj(self):
        """Open a readonly file object to read the stored POTCAR file."""
        file_obj = None
        with self.get_archive() as archive:
            try:
                file_obj = archive.extractfile(archive.members[0])
                yield file_obj
            finally:
                if file_obj:
                    file_obj.close()

    @contextmanager
    def get_file_obj_and_tar_obj(self):
        """Return both decompressed file object and the archive object"""
        file_obj = None
        with self.get_archive() as archive:
            try:
                file_obj = archive.extractfile(archive.members[0])
                yield file_obj, archive
            finally:
                if file_obj:
                    file_obj.close()

    def export_archive(self, archive, dry_run=False):
        """Add the stored POTCAR file to an archive for export."""
        with self.get_file_obj_and_tar_obj() as objects:
            potcar_fo, tar_fo = objects
            arcname = '{}/POTCAR'.format(self.symbol)
            tarinfo = tar_fo.members[0]
            tarinfo.name = arcname
            if not dry_run:
                archive.addfile(tarinfo, fileobj=potcar_fo)
        return tarinfo.name

    def export_file(self, path, dry_run=False):
        """
        Write the contents of the stored POTCAR file to a destination on the local file system.

        :param path: path to the destination file or folder as a Path or string object

        When given a folder, the destination file will be created in a subdirectory with the name of the symbol.
        This is for conveniently exporting multiple files into the same folder structure as the POTCARs are
        distributed in.

        Examples::

            potcar_file = PotcarFileData.get_or_create(<file>)
            assert potcar_file.symbol == 'Si_d'

            potcar_file.export('./POTCAR.Si')
            ## writes to ./POTCAR.Si

            potcar_file.export('./potcars/')
            ## writes to
            ## ./
            ##  |-potcars/
            ##           |-Si_d/
            ##                 |-POTCAR
        """
        path = Path(path)
        if path.is_dir():
            path = path / self.symbol / 'POTCAR'
        if not dry_run:
            # Make sure the directory exists
            path_dir = path.parent
            path_dir.mkdir(parents=True, exist_ok=True)
            with path.open(mode='wb') as dest_fo:
                dest_fo.write(self.get_content())
        return path

    def get_content(self):
        with self.get_file_obj() as potcar_fo:
            return potcar_fo.read()

    def get_pymatgen(self):
        """Create a corresponding pymatgen ``PotcarSingle`` instance."""
        return PotcarSingle(self.get_content())

    @classmethod
    def get_or_create(cls, filepath):
        """Get or create (store) a PotcarFileData node."""
        sha512 = cls.get_file_sha512(filepath)
        if cls.exists(sha512=sha512):
            created = False
            node = cls.find_one(sha512=sha512)
        else:
            created = True
            node = cls(file=filepath)
            node.store()
        return node, created

    @classmethod
    def get_or_create_from_contents(cls, contents):
        """Get or create (store) a PotcarFileData node from a string containing the POTCAR contents."""
        with temp_potcar(contents) as potcar_file:
            return cls.get_or_create(potcar_file)


class PotcarData(Data, PotcarMetadataMixin, VersioningMixin):
    """
    Store enough metadata about a POTCAR file to identify and find it.

    Meant to be used as an input to calculations. This node type holds no
    licenced data and can be freely shared without legal repercussions.
    """

    _query_label = 'potcar'
    _query_type_string = 'data.vasp.potcar.'
    _plugin_type_string = 'data.vasp.potcar.PotcarData.'
    _VERSION = 1

    GROUP_TYPE = 'data.vasp.potcar.family'

    def __init__(self, **kwargs):
        potcar_file_node = kwargs.pop('potcar_file_node', None)
        super(PotcarData, self).__init__(**kwargs)
        if potcar_file_node is not None:
            self.set_potcar_file_node(potcar_file_node)

    def set_potcar_file_node(self, potcar_file_node):
        """Initialize from a PotcarFileData node."""
        self.set_version()
        for attr_name in potcar_file_node.attributes.keys():
            self.set_attribute(attr_name, potcar_file_node.get_attribute(attr_name))

    def find_file_node(self):
        """Find and return the matching PotcarFileData node."""
        return PotcarFileData.find_one(**self.attributes)

    # pylint: disable=arguments-differ,signature-differs
    def store(self, *args, **kwargs):
        """Ensure uniqueness before storing."""
        self.set_version()
        self.verify_unique()
        return super(PotcarData, self).store(*args, **kwargs)

    @classmethod
    def get_or_create(cls, file_node):
        """Get or create (store) a PotcarData node."""

        if cls.exists(sha512=file_node.sha512):
            created = False
            node = cls.find_one(sha512=file_node.sha512)
        else:
            created = True
            node = cls(potcar_file_node=file_node)
            node.store()
        return node, created

    @classmethod
    def get_or_create_from_file(cls, file_path):
        """Get or create (store) a PotcarData node from a POTCAR file."""
        sha512 = PotcarFileData.get_file_sha512(file_path)
        file_node = PotcarFileData.find_one(sha512=sha512) if PotcarFileData.exists(sha512=sha512) else PotcarFileData(file=file_path)
        node, created = cls.get_or_create(file_node)
        if not file_node.is_stored:
            file_node.store()
        return node, created

    @classmethod
    def get_or_create_from_contents(cls, contents):
        """Get or create (store) a PotcarData node from a string containing the POTCAR contents."""
        with temp_potcar(contents) as potcar_file:
            return cls.get_or_create_from_file(str(potcar_file))

    @classmethod
    def file_not_uploaded(cls, file_path):
        sha512 = PotcarFileData.get_file_sha512(file_path)
        return PotcarFileData.find_one(sha512=sha512) if PotcarFileData.exists(sha512=sha512) else namedtuple('potcar', ('uuid'))('-1')

    def get_family_names(self):
        """List potcar families to which this instance belongs."""
        return [group.label for group in Group.query(nodes=self, type_string=self.potcar_family_type_string)]

    @classproperty
    def potcar_family_type_string(cls):  # pylint: disable=no-self-argument
        return cls.GROUP_TYPE

    @classmethod
    def get_potcar_group(cls, group_name):
        """Return the PotcarFamily group with the given name."""
        try:
            group = Group.get(label=group_name, type_string=cls.potcar_family_type_string)
        except NotExistent:
            group = None
        return group

    @classmethod
    def get_potcar_groups(cls, filter_elements=None, filter_symbols=None):
        """
        List all names of groups of type PotcarFamily, possibly with some filters.

        :param filter_elements: list of strings.
               If present, returns only the groups that contains one POTCAR for
               every element present in the list. Default=None, meaning that
               all families are returned. A single element can be passed as a string.
        :param filter_symbols: list of strings with symbols to filter for.
        """
        group_query = QueryBuilder()
        group_query.append(Group,
                           with_node='potcar_data',
                           tag='potcar_data',
                           filters={'type_string': {
                               '==': cls.potcar_family_type_string
                           }},
                           project='*')

        groups = [group_list[0] for group_list in group_query.all()]

        if filter_elements:
            for element in filter_elements:
                idx_has_element = []
                for i, group in enumerate(groups):
                    group_filters = {'label': {'==': group.label}, 'type_string': {'==': cls.potcar_family_type_string}}
                    element_filters = {'attributes.element': {'==': element}}
                    elem_query = QueryBuilder()
                    elem_query.append(Group, tag='family', filters=group_filters)
                    elem_query.append(cls, tag='potcar', with_group='family', filters=element_filters)
                    if elem_query.count() > 0:
                        idx_has_element.append(i)
                groups = [groups[i] for i in range(len(groups)) if i in idx_has_element]

        if filter_symbols:
            for symbol in filter_symbols:
                idx_has_symbol = []
                for i, group in enumerate(groups):
                    group_filters = {'label': {'==': group.label}, 'type_string': {'==': cls.potcar_family_type_string}}
                    symbol_filters = {'attributes.symbol': {'==': symbol}}
                    symbol_query = QueryBuilder()
                    symbol_query.append(Group, tag='family', filters=group_filters)
                    symbol_query.append(cls, tag='potcar', with_group='family', filters=symbol_filters)
                    if symbol_query.count() > 0:
                        idx_has_symbol.append(i)
                groups = [groups[i] for i in range(len(groups)) if i in idx_has_symbol]

        return groups

    @classmethod
    def get_potcars_dict(cls, elements, family_name, mapping=None):
        """
        Get a dictionary {element: ``PotcarData.full_name``} for all given symbols.

        :param elements: The list of symbols to find POTCARs for
        :param family_name: The POTCAR family to be used
        :param mapping: A mapping[element] -> ``full_name``, for example: mapping={'In': 'In', 'As': 'As_d'}

        Exceptions:

         *If the mapping does not contain an item for a given element name, raise a ``ValueError``.
         *If no POTCAR is found for a given element, a ``NotExistent`` error is raised.

        If there are multiple POTCAR with the same ``full_name``, the first one
        returned by ``PotcarData.find()`` will be used.
        """
        if not mapping:
            mapping = {element: element for element in elements}
        group_filters = {'label': {'==': family_name}, 'type_string': {'==': cls.potcar_family_type_string}}
        element_filters = {'attributes.full_name': {'in': [mapping[element] for element in elements]}}
        query = QueryBuilder()
        query.append(Group, tag='family', filters=group_filters)
        query.append(cls, tag='potcar', with_group='family', filters=element_filters)

        result_potcars = {}
        for element in elements:
            if element not in mapping:
                raise ValueError('Potcar mapping must contain an item for each element in the structure, '
                                 'with the full name of the POTCAR file (i.e. "In_d", "As_h").')
            full_name = mapping[element]
            potcars_of_kind = [potcar[0] for potcar in query.all() if potcar[0].full_name == full_name]
            if not potcars_of_kind:
                raise NotExistent('No POTCAR found for full name {} in family {}'.format(full_name, family_name))
            if len(potcars_of_kind) > 1:
                result_potcars[element] = cls.find(family=family_name, full_name=full_name)[0]
            else:
                result_potcars[element] = potcars_of_kind[0]

        return result_potcars

    @classmethod
    def query_by_attrs(cls, query=None, **kwargs):
        family_name = kwargs.pop('family_name', None)
        if family_name:
            group_filters = {'label': {'==': family_name}, 'type_string': {'==': cls.potcar_family_type_string}}
            query = QueryBuilder()
            query.append(Group, tag='family', filters=group_filters)
            query.append(cls, tag=cls._query_label, with_group='family')
        return super(PotcarData, cls).query_by_attrs(query=query, **kwargs)

    @classmethod
    def get_full_names(cls, family_name=None, element=None):
        """
        Gives a set of symbols provided by this family.

        Not every symbol may be supported for every element.
        """
        query = cls.query_by_attrs(family_name=family_name, element=element)
        query.add_projection(cls._query_label, 'attributes.full_name')
        return [name[0] for name in query.all()]

    @classmethod
    def get_potcars_from_structure(cls, structure, family_name, mapping=None):
        """
        Given a POTCAR family name and a AiiDA structure, return a dictionary associating each kind name with its PotcarData object.

        :param structure: An AiiDA structure
        :param family_name: The POTCAR family to be used
        :param mapping: A mapping[kind name] -> ``full_name``, for example: mapping={'In1': 'In', 'In2': 'In_d', 'As': 'As_d'}

        The Dictionary looks as follows::

            {
                kind1.name: PotcarData_for_kind1,
                kind2.name: ...
            }

        This is to make the output of this function suitable for giving directly as input to VaspCalculation.process() instances.

        :raise MultipleObjectsError: if more than one UPF for the same element is found in the group.
        :raise NotExistent: if no UPF for an element in the group is found in the group.


        Example::

            ## using VASP recommended POTCARs
            from aiida_vasp.utils.default_paws import DEFAULT_LDA, DEFAULT_GW
            vasp_process = CalculationFactory('vasp.vasp').process()
            inputs = vasp_process.get_inputs_template()
            inputs.structure = load_node(123)
            inputs.potential = PotcarData.get_potcars_from_structure(
                structure=inputs.structure,
                family_name='PBE',
                mapping=DEFAULT_GW
            )

            ## using custom POTCAR map
            custom_mapping = {
                'In1': 'In',
                'In2': 'In_d',
                'As': 'As_d'
            }
            inputs.potential = PotcarData.get_potcars_from_structure(
                structure=inputs.structure,
                family_name='PBE',
                mapping=custom_mapping
            )
        """
        kind_names = structure.get_kind_names()
        potcar_dict = {kind_name: value
                       for kind_name, value in cls.get_potcars_dict(kind_names, # pylint: disable=unnecessary-comprehension
                                                                    family_name,
                                                                    mapping=mapping).items()}  # yapf: disable
        return potcar_dict

    @classmethod
    def _prepare_group_for_upload(cls, group_name, group_description=None, dry_run=False):
        """Prepare a (possibly new) group to upload a POTCAR family to."""
        if not dry_run:
            group, group_created = Group.objects.get_or_create(label=group_name, type_string=cls.potcar_family_type_string)
        else:
            group = cls.get_potcar_group(group_name)
            group_created = bool(not group)
            if not group:
                group = Group(label=group_name)

        if group.user.pk != get_current_user().pk:
            raise UniquenessError(
                'There is already a POTCAR family group with name {}, but it belongs to user {}, therefore you cannot modify it'.format(
                    group_name, group.user.email))

        if group_description:
            group.description = group_description
        elif group_created:
            raise ValueError('A new POTCAR family {} should be created but no description was given!'.format(group_name))

        return group

    @classmethod
    def upload_potcar_family(cls, source, group_name, group_description=None, stop_if_existing=True, dry_run=False):
        """
        Upload a set of POTCAR potentials as a family.

        :param source: a path containing all POTCAR files to be added.
        :param group_name: the name of the group to create. If it exists and is
            non-empty, a UniquenessError is raised.
        :param group_description: a string to be set as the group description.
            Overwrites previous descriptions, if the group was existing.
        :param stop_if_existing: if True, check for the sha512 of the files and,
            if the file already exists in the DB, raises a MultipleObjectsError.
            If False, simply adds the existing UPFData node to the group.
        :param dry_run: If True, do not change the database.
        """
        group = cls._prepare_group_for_upload(group_name, group_description, dry_run=dry_run)

        potcar_finder = PotcarWalker(source)
        potcar_finder.walk()
        num_files = len(potcar_finder.potcars)
        family_nodes_uuid = [node.uuid for node in group.nodes] if not dry_run else []
        potcars_tried_upload = cls._try_upload_potcars(potcar_finder.potcars, stop_if_existing=stop_if_existing, dry_run=dry_run)
        new_potcars_added = [
            (potcar, created, file_path) for potcar, created, file_path in potcars_tried_upload if potcar.uuid not in family_nodes_uuid
        ]

        for potcar, created, file_path in new_potcars_added:
            if created:
                aiidalogger.debug('New PotcarData node %s created while uploading file %s for family %s', potcar.uuid, file_path,
                                  group_name)
            else:
                aiidalogger.debug('PotcarData node %s used instead of uploading file %s to family %s', potcar.uuid, file_path, group_name)

        if not dry_run:
            group.add_nodes([potcar for potcar, created, file_path in new_potcars_added])

        num_added = len(new_potcars_added)
        num_uploaded = len([item for item in new_potcars_added if item[1]])  # item[1] refers to 'created'

        return num_files, num_added, num_uploaded

    @classmethod
    def _try_upload_potcars(cls, file_paths, stop_if_existing=True, dry_run=False):
        """Given a list of absolute paths to potcar files, try to upload them (or pretend to if dry_run=True)."""
        list_created = []
        for file_path_obj in file_paths:
            file_path = str(file_path_obj)
            try:
                if not dry_run:
                    potcar, created = cls.get_or_create_from_file(file_path)
                else:
                    potcar = cls.file_not_uploaded(file_path)
                    created = bool(potcar.uuid == -1)
                if stop_if_existing and not created:
                    raise ValueError(('A POTCAR with identical SHA512 to {} is already in the DB,'
                                      'therefore it cannot be added with the stop_if_existing kwarg.').format(file_path))
                list_created.append((potcar, created, file_path))
            except KeyError as err:
                print('skipping file {} - uploading raised {}{}'.format(file_path, str(err.__class__), str(err)))
            except AttributeError as err:
                print('skipping file {} - uploading raised {}{}'.format(file_path, str(err.__class__), str(err)))
            except IndexError as err:
                print('skipping file {} - uploading raised {}{}'.format(file_path, str(err.__class__), str(err)))

        return list_created

    @classmethod
    def export_family_folder(cls, family_name, path=None, dry_run=False):
        """
        Export a family of POTCAR nodes into a file hierarchy similar to the one POTCARs are distributed in.

        :param family_name: name of the POTCAR family
        :param path: path to a local directory, either a string or Path object, default to current directory
        :param dry_run: bool, if True, only collect the names of files that would otherwise be written.

        If ``path`` already exists, everything will be written into a subdirectory with the name of the family.
        """
        # Only allow Path or string
        if path is not None:
            if isinstance(path, (Path, str)):
                path = Path(path)
            else:
                raise ValueError('The supplied path is not a Path object or a string.')
        else:
            path = Path()

        if path.exists():
            path = path / family_name
        group = cls.get_potcar_group(family_name)
        all_file_nodes = [potcar.find_file_node() for potcar in group.nodes]
        files_written = []

        with temp_dir() as staging_dir:
            for file_node in all_file_nodes:
                new_file = file_node.export_file(staging_dir, dry_run=dry_run)
                files_written.append(path / new_file.relative_to(staging_dir))
            if not dry_run:
                # copytree uses copy2 which conserves all metadata as well
                shutil.copytree(staging_dir, path)

        return files_written

    @classmethod
    def export_family_archive(cls, family_name, path=None, dry_run=False):
        """Export a family of POTCAR nodes into a compressed archive."""
        # Only allow Path or string
        if path is not None:
            if isinstance(path, (Path, str)):
                path = Path(path)
            else:
                raise ValueError('The supplied path is not a Path object or a string.')
        else:
            path = Path()

        if path.is_dir():
            path = path / family_name

        if not path.suffix:
            name = path.name + '.tar.gz'
            path = path.parent / name

        archive = tarfile.open(str(path), 'w:gz') if not dry_run else None
        group = cls.get_potcar_group(family_name)
        all_file_nodes = [potcar.find_file_node() for potcar in group.nodes]
        files_added = []

        for file_node in all_file_nodes:
            files_added.append(file_node.export_archive(archive, dry_run=dry_run))
        if not dry_run:
            archive.close()
        return path, files_added

    def get_content(self):
        return self.find_file_node().get_content()

    def get_pymatgen(self):
        return self.find_file_node().get_pymatgen()

    @classmethod
    def find(cls, **kwargs):
        """
        Extend :py:meth:`PotcarMetadataMixin.find` with filtering by POTCAR family.

        If no POTCAR is found, raise a ``NotExistent`` exception.

        If multiple POTCAR are found, sort them by:

            * POTCARS belonging to the active user first
            * oldest first
        """
        family = kwargs.pop('family', None)
        if not family:
            return super(PotcarData, cls).find(**kwargs)
        query = cls.query_by_attrs(**kwargs)
        group_filters = {'label': {'==': family}, 'type_string': {'==': cls.potcar_family_type_string}}
        query.append(Group, tag='family', filters=group_filters, with_node=cls._query_label)
        query.add_projection(cls._query_label, '*')
        if not query.count():
            raise NotExistent()
        results = [result[0] for result in query.all()]
        from functools import cmp_to_key
        results.sort(key=cmp_to_key(by_older))
        return results
