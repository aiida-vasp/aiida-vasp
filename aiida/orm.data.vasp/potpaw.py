import os
from aiida.orm import Data
from aiida.orm.querytool import QueryTool
from aiida.tools.codespecific.vasp.io.potcar import parser as pcparser


def import_family(folder):
    fp = os.path.abspath(folder)
    for pawf in os.listdir(fp):
        ap = os.path.join(fp, pawf)
        if os.path.isdir(ap):
            paw = PotpawData.from_folder(ap)
            if not PotpawData.load_paw(family=paw.family, symbol=paw.symbol):
                paw.store_all()


class PotpawData(Data):
    @property
    def symbol(self):
        return self.get_attr('symbol')

    @property
    def element(self):
        return self.get_attr('element')

    @property
    def potcar(self):
        return self.get_attr('potcar')

    @potcar.setter
    def potcar(self, value):
        name = 'POTCAR'
        self.folder.insert_path(value, 'path/'+name)
        self._set_attr('potcar', self.get_abs_path(name))
        attr_dict = pcparser.parse_potcar(value)
        for k, v in attr_dict.iteritems():
            self._set_attr(k, v)

    @property
    def psctr(self):
        return self.get_attr('psctr')

    @psctr.setter
    def psctr(self, value):
        name = 'PSCTR'
        self.folder.insert_path(value, 'path/'+name)
        self._set_attr('psctr', self.get_abs_path(name))

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
        return '<PotpawData: {f}/{s} uuid: {u} (pk: {p})>'.format(
            f=fam,
            s=sym,
            u=self.uuid,
            p=self.pk)
