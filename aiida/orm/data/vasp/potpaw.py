import os
from aiida.orm.data.singlefile import SinglefileData
from aiida.orm import Data
from aiida.orm.querytool import QueryTool


def import_family(folder):
    fp = os.path.abspath(folder)
    for pawf in os.listdir(fp):
        ap = os.path.join(fp, pawf)
        if os.path.isdir(ap):
            if not PotpawData.load_paw(family=paw.family, symbol=paw.symbol):
                paw = PotpawData.from_folder(ap)
                paw.store_all()


class PotpawData(Data):
    @property
    def symbol(self):
        return self.get_attr('symbol')
    @symbol.setter
    def symbol(self, value):
        self._set_attr('symbol', value)

    @property
    def kind(self):
        return self.get_attr('kind')
    @kind.setter
    def kind(self, value):
        self._set_attr('kind', value)

    @property
    def potcar(self):
        return self.get_attr('potcar')
    @potcar.setter
    def potcar(self, value):
        #~ name = 'POTCAR'.format(self.get_attr('symbol'))
        name = 'POTCAR'
        self.folder.insert_path(value, 'path/'+name)
        self._set_attr('potcar', self.get_abs_path(name))

    @property
    def psctr(self):
        return self.get_attr('psctr')
    @psctr.setter
    def psctr(self, value):
        #~ name = 'PSCTR_{}'.format(self.get_attr('symbol'))
        name = 'PSCTR'
        self.folder.insert_path(value, 'path/'+name)
        self._set_attr('psctr', self.get_abs_path(name))

    @property
    def family(self):
        return self.get_attr('family')
    @family.setter
    def family(self, value):
        self._set_attr('family', value)

    @classmethod
    def from_folder(cls, pawpath):
        res = cls()
        ap = os.path.abspath(pawpath)
        p = pawpath.split(os.path.sep)
        symbol = p[-1]
        family = p[-2]
        kind = symbol.split('_')[0]
        res.kind = kind
        res.symbol = symbol
        res.family = family
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
        return '<PotpawData: {f}/{s} uuid: {u} (pk: {p})>'.format(
                f=self.get_attr('family'),
                s=self.get_attr('symbol'),
                u=self.uuid,
                p=self.pk
                )
