# pylint: disable=abstract-method
# explanation: pylint wrongly complains about (aiida) Node not implementing query
"""PAW Pseudopotential data node"""
import os

from aiida.orm import Data
from aiida.orm.querybuilder import QueryBuilder
from aiida.common.exceptions import NotExistent, UniquenessError
from aiida.common.utils import md5_file

from aiida_vasp.utils.io.potcar import PawParser as pcparser


class PawData(Data):
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
    def potcar(self, value):
        name = 'POTCAR'
        self.folder.insert_path(value, 'path/' + name)
        attr_dict = pcparser.parse_potcar(value)
        self._set_attr('md5', md5_file(value))
        for key, val in attr_dict.iteritems():
            self._set_attr(key, val)

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
        """Returns a PAW family group if it exists, otherwise
        raises an exception"""
        from aiida.orm import Group
        return Group.get(name=famname, type_string=cls.group_type)

    @classmethod
    def check_family(cls, name):
        """:py:method: checks wether a PAW family exists.
            :returns: True if exists, False otherwise"""
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
        from aiida.backends.utils import get_automatic_user

        group, group_created = Group.get_or_create(
            name=famname, type_string=cls.group_type)

        if group.user != get_automatic_user():
            raise UniquenessError("There is already a UpfFamily group "
                                  "with name {}, but it belongs to user {},"
                                  " therefore you cannot modify it".format(
                                      famname, group.user.email))
        return group, group_created

    @classmethod
    def get_paw_groups(cls, elements=None, symbols=None, user=None):
        """Find all paw groups containing potentials with the given attributes"""
        from aiida.orm import Group
        from aiida.backends.utils import get_automatic_user
        params = {
            'type_string': cls.group_type,
            'node_attributes': {
                'element': elements,
                'symbol': symbols
            }
        }
        if user:
            params['user'] = user
        else:
            params['user'] = get_automatic_user()

        res = Group.query(**params)
        groups = [(g.name, g) for g in res]
        # Sort by name
        groups.sort()
        # Return the groups, without name
        return [i[1] for i in groups]

    @classmethod
    def import_family(cls,
                      folder,
                      familyname=None,
                      family_desc=None,
                      store=True,
                      stop_if_existing=False):
        """Import a family from a folder like the ones distributed with VASP,
        usually named potpaw_XXX"""
        from aiida.common import aiidalogger

        ffound = []
        fupl = []
        family_path = os.path.abspath(folder)
        # ~ ffname = os.path.basename(
        # ~ os.path.dirname(folder)).replace('potpaw_', '')
        # ~ famname = familyname or ffname

        group, group_created = cls.get_or_create_famgroup(familyname)

        # Always update description, even if the group already existed
        group.description = family_desc

        paw_list = cls._find_paws(family_path, ffound, group, group_created)

        if stop_if_existing:
            for pawinfo in paw_list:
                if not pawinfo[1]:
                    raise ValueError("A PAW with identical MD5 to "
                                     '' + pawinfo[2] + " cannot be added with "
                                     "stop_if_existing")

        for pawinfo in paw_list:
            paw = pawinfo[0]
            created = pawinfo[1]
            path = pawinfo[2]
            if store:
                if created:
                    paw.store_all()
                    aiidalogger.debug("New node %s created for file %s",
                                      paw.uuid, path)
                    fupl.append(path)
                else:
                    aiidalogger.debug("Reusing node %s for file %s", paw.uuid,
                                      path)

        if store:
            if group_created:
                group.store()
                aiidalogger.debug("New PAW family goup %s created", group.uuid)
            group.add_nodes(i[0] for i in paw_list)
        else:
            print map(repr, [i[0] for i in paw_list])

        return ffound, fupl

    @classmethod
    def _find_paws(cls, family_path, ffound, group, group_created):
        """Go through a directory containing a family of paws and collect individual pseudopotentials"""
        paw_list = []
        for pawf in os.listdir(family_path):
            try:
                subfolder_path = os.path.join(family_path, pawf)
                potcar_path = os.path.join(subfolder_path, 'POTCAR')
                if os.path.isdir(subfolder_path) and os.path.exists(
                        potcar_path):
                    ffound.append(pawf)
                    paw, paw_created = cls.get_or_create(subfolder_path)
                    # ~ paw._set_attr('family', famname)
                    # ~ upload = paw_created
                    # enforce group-wise uniqueness of symbols
                    in_group = False
                    if not group_created:
                        in_group = bool(
                            cls.load_paw(
                                group=group, symbol=paw.symbol, silent=True))
                    if not in_group:
                        paw_list.append((paw, paw_created, pawf))
            except Exception:  # pylint: disable=broad-except
                import sys
                err = sys.exc_info()[1]
                print 'WARNING: skipping ' + os.path.abspath(
                    os.path.join(family_path, pawf))
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
        paw = cls.query(
            dbattributes__key='md5', dbattributes__tval=md5new).first()
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
        py:method:: load_paw([family=None][, element=None][, symbol=None])
        Load PawData nodes from the databank. Use kwargs to filter.

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
            query_builder.append(cls=cls)
            filters = {}
            for key, value in kwargs.iteritems():
                filters[key] = {'=': value}
            query_builder.append(filters=filters)
            res = list(query_builder.iterall())
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
        return '<PawData: {s} uuid: {u} (pk: {p})>'.format(
            s=sym, u=self.uuid, p=self.pk)
