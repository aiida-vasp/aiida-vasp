import os
from aiida.orm import Data
from aiida.orm.querytool import QueryTool
from aiida.tools.codespecific.vasp.io.potcar import PawParser as pcparser



class PawData(Data):
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
        for k, v in attr_dict.iteritems():
            if k=='paw_date':
                self._set_attr(k, v.strftime('%Y-%m-%d'))
            else:
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
    def import_family(cls, folder, family_override=None):
        fp = os.path.abspath(folder)
        ffname = os.path.basename(os.path.dirname(folder)).replace('potpaw_', '')
        famname = family_override or ffname
        for pawf in os.listdir(fp):
            ap = os.path.join(fp, pawf)
            if os.path.isdir(ap):
                paw = cls.from_folder(ap)
                paw._set_attr('family', famname)
                if not cls.load_paw(family=famname, symbol=paw.symbol):
                    paw.store_all()
        @property
        def xc_type(self):
            return self.get_attr('xc_type')

    @classmethod
    def from_folder(cls, pawpath):
        res = cls()
        ap = os.path.abspath(pawpath)
        res.potcar = os.path.join(ap, 'POTCAR')
        res.psctr = os.path.join(ap, 'PSCTR')
        return res

    @classmethod
    def load_paw(cls, **kwargs):
        q = QueryTool()
        q.set_class(cls)
        for k, v in kwargs.iteritems():
            q.add_attr_filter(k, '=', v)
        return list(q.run_query())

    def __repr__(self):
        try:
            fam = self.get_attr('family')
        except AttributeError:
            fam = '<family: (unset)>'
        try:
            sym = self.get_attr('symbol')
        except AttributeError:
            sym = '<symbol: (unset)>'
        return '<PawData: {f}/{s} uuid: {u} (pk: {p})>'.format(
            f=fam,
            s=sym,
            u=self.uuid,
            p=self.pk)
