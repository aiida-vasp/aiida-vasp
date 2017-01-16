import os
from aiida.orm import Data
from aiida.orm.querytool import QueryTool
from aiida.tools.codespecific.vasp.io.potcar import PawParser as pcparser
from aiida.common.exceptions import NotExistent, UniquenessError
from aiida.common.utils import md5_file


class PawData(Data):
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
        self.folder.insert_path(value, 'path/'+name)
        attr_dict = pcparser.parse_potcar(value)
        self._set_attr('md5', md5_file(value))
        for k, v in attr_dict.iteritems():
            self._set_attr(k, v)

    @property
    def psctr(self):
        return self.get_abs_path('PSCTR')

    @psctr.setter
    def psctr(self, value):
        name = 'PSCTR'
        self.folder.insert_path(value, 'path/'+name)

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
        '''Returns a PAW family group if it exists, otherwise
        raises an exception'''
        from aiida.orm import Group
        return Group.get(name=famname, type_string=cls.group_type)

    @classmethod
    def check_family(cls, name):
        ''':py:method: checks wether a PAW family exists.
            :returns: True if exists, False otherwise'''
        exists = False
        try:
            group = cls.get_famgroup(name)
            exists = bool(group)
        except NotExistent:
            exists = False
        return exists

    @classmethod
    def get_or_create_famgroup(cls, famname):
        '''Returns a PAW family group, creates it if it didn't exists'''
        from aiida.orm import Group
        from aiida.djsite.utils import get_automatic_user

        # TODO: maybe replace with Group.get_or_create?
        try:
            group = Group.get(name=famname, type_string=cls.group_type)
            group_created = False
        except NotExistent:
            group = Group(name=famname, type_string=cls.group_type,
                          user=get_automatic_user())
            group_created = True

        if group.user != get_automatic_user():
            raise UniquenessError("There is already a UpfFamily group "
                                  "with name {}, but it belongs to user {},"
                                  " therefore you cannot modify it".format(
                                      famname, group.user.email))
        return group, group_created

    @classmethod
    def get_paw_groups(cls, elements=set(), symbols=set(), user=None):
        from aiida.orm import Group
        from aiida.djsite.utils import get_automatic_user
        params = {'type_string': cls.group_type,
                  'node_attributes': {
                      'element': elements,
                      'symbol': symbols}
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
    def import_family(cls, folder, familyname=None,
                      family_desc=None, store=True, stop_if_existing=False):
        '''Import a family from a folder like the ones distributed with VASP,
        usually named potpaw_XXX'''
        from aiida.common import aiidalogger

        ffound = []
        fupl = []
        paw_list = []
        fp = os.path.abspath(folder)
        # ~ ffname = os.path.basename(
            # ~ os.path.dirname(folder)).replace('potpaw_', '')
        # ~ famname = familyname or ffname

        group, group_created = cls.get_or_create_famgroup(familyname)

        # Always update description, even if the group already existed
        group.description = family_desc

        for pawf in os.listdir(fp):
            try:
                ap = os.path.join(fp, pawf)
                pp = os.path.join(ap, 'POTCAR')
                if os.path.isdir(ap) and os.path.exists(pp):
                    ffound.append(pawf)
                    paw, paw_created = cls.get_or_create(ap)
                    # ~ paw._set_attr('family', famname)
                    # ~ upload = paw_created
                    # enforce group-wise uniqueness of symbols
                    in_group = False
                    if not group_created:
                        in_group = bool(cls.load_paw(group=group,
                                        symbol=paw.symbol, silent=True))
                    if not in_group:
                        paw_list.append((paw, paw_created, pawf))
            except:
                import sys
                e = sys.exc_info()[1]
                print 'WARNING: skipping ' + os.path.abspath(
                    os.path.join(fp, pawf))
                print '  ' + e.__class__.__name__ + ': ' + e.message

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
                    aiidalogger.debug("New node {} created for file {}".format(
                        paw.uuid, path))
                    fupl.append(path)
                else:
                    aiidalogger.debug("Reusing node {} for file {}".format(
                        paw.uuid, path))

        if store:
            if group_created:
                group.store()
                aiidalogger.debug("New PAW family goup {} created".format(
                    group.uuid))
            group.add_nodes(i[0] for i in paw_list)
        else:
            print map(repr, [i[0] for i in paw_list])

        return ffound, fupl

    @property
    def xc_type(self):
        return self.get_attr('xc_type')

    @classmethod
    def get_or_create(cls, pawpath, store=False, psctr=None):
        if os.path.isdir(pawpath):
            pp = os.path.join(pawpath, 'POTCAR')
            isdir = True
        else:
            pp = pawpath
            isdir = False
        md5new = md5_file(pp)
        paw = cls.query(dbattributes__key='md5',
                        dbattributes__tval=md5new).first()
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
        res = cls()
        res.potcar = potpath
        if ctrpath:
            res.psctr = ctrpath
        return res

    @classmethod
    def from_folder(cls, pawpath):
        res = cls()
        ap = os.path.abspath(pawpath)
        cp = os.path.join(ap, 'PSCTR')
        res.potcar = os.path.join(ap, 'POTCAR')
        if os.path.isfile(cp):
            res.psctr = cp
        return res

    @classmethod
    def _node_filter(cls, **kwargs):
        def node_filter(node):
            for k, v in kwargs.iteritems():
                if not (node.get_attr(k) == v):
                    return False
            return True
        return node_filter

    @classmethod
    def load_paw(cls, **kwargs):
        '''
        py:method:: load_paw([family=None][, element=None][, symbol=None])
        Load PawData nodes from the databank. Use kwargs to filter.

        :return: a list of PawData instances
        :rtype: list
        :key str family: Filter by family
        :key str element: Filter by chemical symbol
        :key str symbol: Filter by PAW symbol (example: As vs. As_d)
        :raises ValueError: if no PAWs are found
        '''
        usage_msg = 'use import_family or from_folder to import PAWs'
        error_msg = 'no PAWs found for the given kwargs!\n'+usage_msg
        group = kwargs.pop('group', None)
        family = kwargs.pop('family', None)
        silent = kwargs.pop('silent', None)
        if not (group or family):
            q = QueryTool()
            q.set_class(cls)
            for k, v in kwargs.iteritems():
                q.add_attr_filter(k, '=', v)
            res = list(q.run_query())
        else:
            if family:
                group, created = cls.get_or_create_famgroup(family)
            elif group:
                created = not group._is_stored
            try:
                paw_filter = cls._node_filter(**kwargs)
                res = filter(paw_filter, group.nodes)
            except ValueError as e:
                if silent:
                    res = []
                elif created:
                    raise NotExistent('No family with that name exists')
                else:
                    raise e

        if not res and not silent:
            raise ValueError(error_msg)
        return res

    def __repr__(self):
        try:
            sym = self.get_attr('symbol')
        except AttributeError:
            sym = '<symbol: (unset)>'
        return '<PawData: {s} uuid: {u} (pk: {p})>'.format(
            s=sym,
            u=self.uuid,
            p=self.pk)
