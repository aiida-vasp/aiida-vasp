from aida.orm.data.singlefile import SinglefileData
from aiida.orm import Data


class PAWFileData(Data):
    @property
    def symbol(self):
        return self.get_attr('symbol')
    @value.setter
    def symbol(self, value):
        self._set_attr('symbol')
