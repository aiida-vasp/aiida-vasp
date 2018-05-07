# pylint: disable=abstract-method
# explanation: pylint wrongly complains about (aiida) Node not implementing query
"""
Legacy PAW Pseudopotential data node.

THIS MODULE IS DEPRECATED AND MAY BE REMOVED IN FUTURE VERSIONS.
"""
import os

from aiida.common.exceptions import NotExistent, UniquenessError
from aiida.common.utils import md5_file
from aiida.orm import Data
from aiida.orm.querybuilder import QueryBuilder


class LegacyPawData(Data):
    """Holds the files and metadata that make up a VASP PAW format pseudopotential"""
    group_type = 'data.vasp.paw.family'

    @property
    def symbol(self):
        return self.get_attr('symbol')

    @property
    def element(self):
        return self.get_attr('element')

    @property
    def potcar(self):
        return self.get_abs_path('POTCAR')

    @potcar.setter
    def potcar(self, value):  # pylint: disable=no-self-use,unused-argument
        raise DeprecationWarning('PawData is deprecated and read-only')

    @property
    def psctr(self):
        return self.get_abs_path('PSCTR')

    @psctr.setter
    def psctr(self, value):
        name = 'PSCTR'
        self.folder.insert_path(value, 'path/' + name)

    @property
    def family(self):
        return self.get_attr('family')

    @property
    def valence(self):
        return self.get_attr('valence')

    @property
    def mass(self):
        return self.get_attr('mass')

    @property
    def paw_date(self):
        return self.get_attr('paw_date')

    @classmethod
    def get_famgroup(cls, famname):
        """Returns a PAW family group if it exists, otherwise raises an exception."""
        from aiida.orm import Group
        return Group.get(name=famname, type_string=cls.group_type)

    @classmethod
    def check_family(cls, name):
        """
        :py:method: checks wether a PAW family exists.

        :returns: True if exists, False otherwise.
        """
        exists = False
        try:
            group = cls.get_famgroup(name)
            exists = bool(group)
        except NotExistent:
            exists = False
        return exists

    @classmethod
    def get_or_create_famgroup(cls, famname):
        """Returns a PAW family group, creates it if it didn't exists"""
        from aiida.orm import Group
        from aiida_vasp.utils.aiida_utils import get_current_user

        group, group_created = Group.get_or_create(name=famname, type_string=cls.group_type)

        if group.user.pk != get_current_user().pk:
            raise UniquenessError("There is already a UpfFamily group "
                                  "with name {}, but it belongs to user {},"
                                  " therefore you cannot modify it".format(famname, group.user.email))
        return group, group_created

    @classmethod
    def get_paw_groups(cls, elements=None, symbols=None, user=None):
        """Find all paw groups containing potentials with the given attributes"""
        from aiida.orm import Group
        from aiida_vasp.utils.aiida_utils import get_current_user
        params = {'type_string': cls.group_type, 'node_attributes': {'element': elements, 'symbol': symbols}}
        if user:
            params['user'] = user
        else:
            params['user'] = get_current_user()

        res = Group.query(**params)
        groups = [(g.name, g) for g in res]
        # Sort by name
        groups.sort()
        # Return the groups, without name
        return [i[1] for i in groups]

    @classmethod
    def import_family(cls, folder, familyname=None, family_desc=None, store=True, stop_if_existing=False):  # pylint: disable=unused-argument
        """Import a family from a folder like the ones distributed with VASP, usually named potpaw_XXX."""
        raise DeprecationWarning('PawData is deprecated and read-only')

    # pylint: disable=arguments-differ,unused-argument
    def store(self, *args, **kwargs):
        raise DeprecationWarning('PawData is deprecated and read-only')

    @classmethod
    def _find_paws(cls, family_path, ffound, group, group_created):
        """Go through a directory containing a family of paws and collect individual pseudopotentials"""
        paw_list = []
        for pawf in os.listdir(family_path):
            try:
                subfolder_path = os.path.join(family_path, pawf)
                potcar_path = os.path.join(subfolder_path, 'POTCAR')
                if os.path.isdir(subfolder_path) and os.path.exists(potcar_path):
                    ffound.append(pawf)
                    paw, paw_created = cls.get_or_create(subfolder_path)
                    # ~ paw._set_attr('family', famname)
                    # ~ upload = paw_created
                    # enforce group-wise uniqueness of symbols
                    in_group = False
                    if not group_created:
                        in_group = bool(cls.load_paw(group=group, symbol=paw.symbol, silent=True))
                    if not in_group:
                        paw_list.append((paw, paw_created, pawf))
            except Exception:  # pylint: disable=broad-except
                import sys
                err = sys.exc_info()[1]
                print 'WARNING: skipping ' + os.path.abspath(os.path.join(family_path, pawf))
                print '  ' + err.__class__.__name__ + ': ' + err.message
        return paw_list

    @property
    def xc_type(self):
        return self.get_attr('xc_type')

    @classmethod
    def get_or_create(cls, pawpath, store=False, psctr=None):
        """Fetch a paw data node from the DB if it exists, else create it"""
        if os.path.isdir(pawpath):
            potcar_path = os.path.join(pawpath, 'POTCAR')
            isdir = True
        else:
            potcar_path = pawpath
            isdir = False
        md5new = md5_file(potcar_path)
        paw = cls.query(dbattributes__key='md5', dbattributes__tval=md5new).first()
        created = False
        if not paw:
            if isdir:
                paw = cls.from_folder(pawpath)
            else:
                paw = cls.from_potcar(pawpath, ctrpath=psctr)
            if store:
                paw.store_all()
            created = True
        return paw, created

    @classmethod
    def from_potcar(cls, potpath, ctrpath=None):
        """Create a paw data node from a POTCAR file and optionally a PSCTR file"""
        res = cls()
        res.potcar = potpath
        if ctrpath:
            res.psctr = ctrpath
        return res

    @classmethod
    def from_folder(cls, pawpath):
        """Create a paw data node from a directory"""
        res = cls()
        abs_path = os.path.abspath(pawpath)
        psctr_path = os.path.join(abs_path, 'PSCTR')
        res.potcar = os.path.join(abs_path, 'POTCAR')
        if os.path.isfile(psctr_path):
            res.psctr = psctr_path
        return res

    @classmethod
    def _node_filter(cls, **kwargs):
        """
        Create and return a node filtering function

        :return: True if all kwargs match node attributes
        """

        def node_filter(node):
            for key, value in kwargs.iteritems():
                if not node.get_attr(key) == value:
                    return False
            return True

        return node_filter

    @classmethod
    def load_paw(cls, **kwargs):
        """
        Loads PawData nodes from the databank, use kwargs to filter.

        :return: a list of PawData instances
        :rtype: list
        :key str family: Filter by family
        :key str element: Filter by chemical symbol
        :key str symbol: Filter by PAW symbol (example: As vs. As_d)
        :raises ValueError: if no PAWs are found
        """
        usage_msg = 'use import_family or from_folder to import PAWs'
        error_msg = 'no PAWs found for the given kwargs!\n' + usage_msg
        group = kwargs.pop('group', None)
        family = kwargs.pop('family', None)
        silent = kwargs.pop('silent', None)
        if not (group or family):
            query_builder = QueryBuilder()
            query_builder.append(cls, tag='paw')
            filters = {}
            for key, value in kwargs.iteritems():
                filters['attributes.{}'.format(key)] = {'==': value}
            query_builder.add_filter('paw', filters)
            res = [i[0] for i in query_builder.all()]
        else:
            if family:
                group, created = cls.get_or_create_famgroup(family)
            elif group:
                created = not group.is_stored  # pylint: disable=protected-access
            try:
                paw_filter = cls._node_filter(**kwargs)
                res = filter(paw_filter, group.nodes)
            except ValueError as err:
                if silent:
                    res = []
                elif created:
                    raise NotExistent('No family with that name exists')
                else:
                    raise err

        if not res and not silent:
            raise ValueError(error_msg)
        return res

    def __repr__(self):
        try:
            sym = self.get_attr('symbol')
        except AttributeError:
            sym = '<symbol: (unset)>'
        return '<PawData: {s} uuid: {u} (pk: {p})>'.format(s=sym, u=self.uuid, p=self.pk)
