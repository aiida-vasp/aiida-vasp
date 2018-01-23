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
import re
import hashlib
import tarfile
import tempfile
import shutil
from contextlib import contextmanager
from collections import namedtuple

from py import path as py_path  # pylint: disable=no-name-in-module,no-member
from pymatgen.io.vasp import PotcarSingle
from aiida.backends.utils import get_automatic_user
from aiida.common import aiidalogger
from aiida.common.utils import classproperty
from aiida.common.exceptions import UniquenessError, MultipleObjectsError, NotExistent
from aiida.orm import Group
from aiida.orm.data import Data
from aiida.orm.querybuilder import QueryBuilder

from aiida_vasp.data.archive import ArchiveData


def normalize_potcar_contents(potcar_contents):
    """Normalize whitespace in a POTCAR given as a string."""
    normalized = re.sub(r'[ \t]+', r' ', potcar_contents)  # multiple spaces
    normalized = re.sub(r'[\n\r]\s*', r'\n', normalized)  # line breaks and spaces afterwards / empty lines
    normalized = re.sub(r'^\s*', r'', normalized)  # spaces / empty lines at the very beginning
    normalized = re.sub(r'\s*$', r'\n', normalized)  # trailing endline
    return normalized


def md5_potcar(potcar_contents):
    """Hash the contents of a POTCAR file (given as str)."""
    md5_hash = hashlib.md5()
    md5_hash.update(normalize_potcar_contents(potcar_contents))
    return md5_hash.hexdigest()


@contextmanager
def temp_dir():
    """Temporary directory context manager that deletes the tempdir after use."""
    try:
        tempdir = tempfile.mkdtemp()
        yield py_path.local(tempdir)
    finally:
        shutil.rmtree(tempdir)


@contextmanager
def temp_potcar(contents):
    """Temporary POTCAR file from contents."""
    with temp_dir() as tempdir:
        potcar_file = tempdir.join('POTCAR')
        with potcar_file.open('w') as potcar_fo:
            potcar_fo.write(contents)
        yield potcar_file


class PotcarMetadataMixin(object):
    """Provide common Potcar metadata access and querying functionality."""
    _query_label = 'label'

    @classmethod
    def query_by_attrs(cls, **kwargs):
        """Find a Data node by attributes."""
        label = cls._query_label
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
        other_attrs.pop('original_filename')
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

    _query_label = 'potcar_file'
    _query_type_string = 'data.vasp.potcar_file.'
    _plugin_type_string = 'data.vasp.potcar_file.PotcarFileData'

    def set_file(self, filepath):
        """Initialize from a file path."""
        self.add_file(filepath)

    def set_contents(self, contents):
        """Initialize from a string."""
        with temp_potcar(contents) as potcar_file:
            self.add_file(str(potcar_file))

    def add_file(self, src_abs, dst_filename=None):
        """Add the POTCAR file to the archive and set attributes."""
        if self._filelist:
            raise AttributeError('Can only hold one POTCAR file')
        super(PotcarFileData, self).add_file(src_abs, dst_filename)
        self._set_attr('md5', self.get_file_md5(src_abs))
        potcar = PotcarSingle.from_file(src_abs)
        self._set_attr('title', potcar.keywords['TITEL'])
        self._set_attr('functional', potcar.functional)
        self._set_attr('element', potcar.element)
        self._set_attr('symbol', potcar.symbol)
        self._set_attr('original_filename', src_abs)

    @classmethod
    def get_file_md5(cls, path):
        """Get the md5 sum for a POTCAR file (after whitespace normalization)."""
        path = py_path.local(path)
        with path.open('r') as potcar_fo:
            md5 = md5_potcar(potcar_fo.read())
        return md5

    @classmethod
    def get_contents_md5(cls, contents):
        """Get the md5 sum for the contents of a POTCAR file (after normalization)."""
        return md5_potcar(contents)

    def store(self, with_transaction=True):
        """Ensure uniqueness and existence of a matching PotcarData node before storing."""
        _ = PotcarData.get_or_create(self)
        self.verify_unique()
        return super(PotcarFileData, self).store(with_transaction=with_transaction)

    @contextmanager
    def get_file_obj(self):
        """Open a readonly file object to read the stored POTCAR file."""
        file_obj = None
        try:
            file_obj = self.archive.extractfile(self.archive.members[0])
            yield file_obj
        finally:
            if file_obj:
                file_obj.close()

    def export_archive(self, archive, dry_run=False):
        """Add the stored POTCAR file to an archive for export."""
        with self.get_file_obj() as potcar_fo:
            arcname = '{}/POTCAR'.format(self.symbol)
            tarinfo = self.archive.members[0]
            tarinfo.name = arcname
            if not dry_run:
                archive.addfile(tarinfo, fileobj=potcar_fo)
        return tarinfo.name

    def export_file(self, path, dry_run=False):
        """
        Write the contents of the stored POTCAR file to a destination on the local file system.

        :param path: path to the destination file or folder

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
        path = py_path.local(path)
        if path.isdir():
            path = path.join(self.symbol, 'POTCAR')
        if not dry_run:
            with path.open(mode='w', ensure=True) as dest_fo:
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
        md5 = cls.get_file_md5(filepath)
        if cls.exists(md5=md5):
            created = False
            node = cls.find(md5=md5)
        else:
            created = True
            node = cls(file=filepath)
            node.store()
        return node, created

    @classmethod
    def get_or_create_from_contents(cls, contents):
        """Get or create (store) a PotcarFileData node from a string containing the POTCAR contents."""
        with temp_potcar(contents) as potcar_file:
            return cls.get_or_create(str(potcar_file))


class PotcarData(Data, PotcarMetadataMixin):
    """
    Store enough metadata about a POTCAR file to identify and find it.

    Meant to be used as an input to calculations. This node type holds no
    licenced data and can be freely shared without legal repercussions.
    """

    _query_label = 'potcar'
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
        md5 = PotcarFileData.get_file_md5(file_path)
        file_node = PotcarFileData.find(md5=md5) if PotcarFileData.exists(md5=md5) else PotcarFileData(file=file_path)
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
        md5 = PotcarFileData.get_file_md5(file_path)
        return PotcarFileData.find(md5=md5) if PotcarFileData.exists(md5=md5) else namedtuple('potcar', ('uuid'))('-1')

    def get_family_names(self):
        """List potcar families to which this instance belongs."""
        return [group.name for group in Group.query(nodes=self, type_string=self.potcar_family_type_string)]

    @classproperty
    def potcar_family_type_string(cls):  # pylint: disable=no-self-argument
        return cls.GROUP_TYPE

    @classmethod
    def get_potcar_group(cls, group_name):
        """Return the PotcarFamily group with the given name."""
        try:
            group = Group.get(name=group_name, type_string=cls.potcar_family_type_string)
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
        group_query.append(cls, tag='potcar_data')
        group_query.append(
            Group, group_of='potcar_data', tag='potcar_family', filters={'type': {
                '==': cls.potcar_family_type_string
            }}, project='*')

        groups = [group_list[0] for group_list in group_query.all()]

        if filter_elements:
            for element in filter_elements:
                idx_has_element = []
                for i, group in enumerate(groups):
                    group_filters = {'name': {'==': group.name}, 'type': {'==': cls.potcar_family_type_string}}
                    element_filters = {'attributes.element': {'==': element}}
                    elem_query = QueryBuilder()
                    elem_query.append(Group, tag='family', filters=group_filters)
                    elem_query.append(cls, tag='potcar', member_of='family', filters=element_filters)
                    if elem_query.count() > 0:
                        idx_has_element.append(i)
                groups = [groups[i] for i in range(len(groups)) if i in idx_has_element]

        if filter_symbols:
            for symbol in filter_symbols:
                idx_has_symbol = []
                for i, group in enumerate(groups):
                    group_filters = {'name': {'==': group.name}, 'type': {'==': cls.potcar_family_type_string}}
                    symbol_filters = {'attributes.symbol': {'==': symbol}}
                    symbol_query = QueryBuilder()
                    symbol_query.append(Group, tag='family', filters=group_filters)
                    symbol_query.append(cls, tag='potcar', member_of='family', filters=symbol_filters)
                    if symbol_query.count() > 0:
                        idx_has_symbol.append(i)
                groups = [groups[i] for i in range(len(groups)) if i in idx_has_symbol]

        return groups

    @classmethod
    def get_potcars_dict(cls, elements, family_name):
        """
        Get a dictionary {element: POTCAR} for all given symbols.

        :param symbols: The list of symbols to find POTCARs for
        :param family_name: The POTCAR family to be used
        """
        group_filters = {'name': {'==': family_name}, 'type': {'==': cls.potcar_family_type_string}}
        element_filters = {'attributes.element': {'in': elements}}
        query = QueryBuilder()
        query.append(Group, tag='family', filters=group_filters)
        query.append(cls, tag='potcar', member_of='family', filters=element_filters)

        result_potcars = {}
        for element in elements:
            potcars_of_kind = [potcar[0] for potcar in query.all() if potcar[0].element == element]
            if not potcars_of_kind:
                raise NotExistent('No POTCAR found for symbol {} in family {}'.format(element, family_name))
            elif len(potcars_of_kind) > 1:
                raise MultipleObjectsError('More than one POTCAR for symbol {} found in family {}'.format(element, family_name))
            result_potcars[potcars_of_kind[0].element] = potcars_of_kind[0]

        return result_potcars

    @classmethod
    def get_potcars_from_structure(cls, structure, family_name):
        """
        Given a POTCAR family name and a AiiDA structure, return a dictionary associating each kind name with its UpfData object.

        The Dictionary looks as follows::

            {
                (element1, ): PotcarData_for_kind1,
                (element2, ): ...
            }

        This is to make the output of this function suitable for giving directly as input to VaspCalculation.process() instances.

        :raise MultipleObjectsError: if more than one UPF for the same element is found in the group.
        :raise NotExistent: if no UPF for an element in the group is found in the group.
        """
        elements_to_name = {kind.symbol: kind.name for kind in structure.kinds}
        return {(elements_to_name[element],): potcar
                for element, potcar in cls.get_potcars_dict(elements_to_name.keys(), family_name).items()}

    @classmethod
    def _prepare_group_for_upload(cls, group_name, group_description=None, dry_run=False):
        """Prepare a (possibly new) group to upload a POTCAR family to."""
        if not dry_run:
            group, group_created = Group.get_or_create(name=group_name, type_string=cls.potcar_family_type_string)
        else:
            group = cls.get_potcar_group(group_name)
            group_created = bool(not group)
            if not group:
                group = Group(name=group_name)

        if group.user != get_automatic_user():
            raise UniquenessError(
                'There is already a PotcarFamily group with name {}, but it belongs to user {}, therefore you cannot modify it'.format(
                    group_name, group.user.email))

        if group_description:
            group.description = group_description
        elif group_created:
            raise ValueError('A new PotcarGroup {} should be created but no description was given!'.format(group_name))

        return group

    @classmethod
    def upload_potcar_family(cls, source, group_name, group_description=None, stop_if_existing=True, dry_run=False):
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
        :param dry_run: If True, do not change the database.
        """
        group = cls._prepare_group_for_upload(group_name, group_description, dry_run=dry_run)

        source = py_path.local(source)
        potcars_found = cls.recursive_upload_potcar([source], stop_if_existing=stop_if_existing, dry_run=dry_run)
        num_files = len(potcars_found)
        family_nodes_uuid = [node.uuid for node in group.nodes] if not dry_run else []
        potcars_found = [
            (potcar, created, file_path) for potcar, created, file_path in potcars_found if potcar.uuid not in family_nodes_uuid
        ]

        for potcar, created, file_path in potcars_found:
            if created:
                aiidalogger.debug('New PotcarData node %s created while uploading file %s for family %s', potcar.uuid, file_path,
                                  group_name)
            else:
                aiidalogger.debug('PotcarData node %s used instead of uploading file %s to family %s', potcar.uuid, file_path, group_name)

        if not dry_run:
            group.add_nodes([potcar for potcar, created, file_path in potcars_found])

        num_uploaded = len(potcars_found)

        return num_files, num_uploaded

    @classmethod
    def recursive_upload_potcar(cls, folders, stop_if_existing=True, dry_run=False, depth=0):
        """Recursively search and upload POTCAR files in a folder or archive."""
        if depth >= 3:
            return []
        list_created = []
        for subpath in folders:
            if subpath.isdir():
                list_created.extend(cls.recursive_upload_potcar(subpath.listdir(), stop_if_existing, dry_run=dry_run, depth=depth + 1))
            elif subpath.isfile():
                if tarfile.is_tarfile(str(subpath)):
                    with temp_dir() as staging_dir:
                        with tarfile.TarFile(str(subpath)) as potcar_archive:
                            potcar_archive.extractall(str(staging_dir))
                        list_created.extend(
                            cls.recursive_upload_potcar(staging_dir.listdir(), stop_if_existing, dry_run=dry_run, depth=depth))
                elif 'POTCAR' in subpath.basename:
                    if not dry_run:
                        potcar, created = cls.get_or_create_from_file(str(subpath))
                    else:
                        potcar = cls.file_not_uploaded(str(subpath))
                        created = bool(potcar.uuid == -1)
                    if stop_if_existing and not created:
                        raise ValueError('A POTCAR with identical MD5 to {} cannot be added with the stop_if_existing kwarg.'.format(
                            str(subpath)))
                    list_created.append((potcar, created, str(subpath)))

        return list_created

    @classmethod
    def export_family_folder(cls, family_name, path='.', dry_run=False):
        """
        Export a family of POTCAR nodes into a file hierarchy similar to the one POTCARs are distributed in.

        :param family_name: name of the POTCAR family
        :param path: path to a local directory
        :param dry_run: bool, if True, only collect the names of files that would otherwise be written.

        If ``path`` already exists, everything will be written into a subdirectory with the name of the family.
        """
        path = py_path.local(path)
        if path.exists():
            path = path.join(family_name)
        group = cls.get_potcar_group(family_name)
        all_file_nodes = [potcar.find_file_node() for potcar in group.nodes]
        files_written = []

        with temp_dir() as staging_dir:
            for file_node in all_file_nodes:
                new_file = file_node.export_file(staging_dir, dry_run=dry_run)
                files_written.append(path.join(new_file.relto(staging_dir)))
            if not dry_run:
                staging_dir.copy(path, stat=True)

        return files_written

    @classmethod
    def export_family_archive(cls, family_name, path='.', dry_run=False):
        """Export a family of POTCAR nodes into a compressed archive."""
        path = py_path.local(path)

        if path.isdir():
            path = path.join(family_name)

        if not path.ext:
            path = path.dirpath().join(path.basename + '.tar.gz')

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
        """Extend :py:meth:`PotcarMetadataMixin.find` with filtering by POTCAR family."""
        family = kwargs.pop('family', None)
        if not family:
            return super(PotcarData, cls).find(**kwargs)
        query = cls.query_by_attrs(**kwargs)
        group_filters = {'name': {'==': family}, 'type': {'==': cls.potcar_family_type_string}}
        query.append(Group, tag='family', filters=group_filters, group_of=cls._query_label)
        query.add_projection(cls._query_label, '*')
        return query.one()[0]
