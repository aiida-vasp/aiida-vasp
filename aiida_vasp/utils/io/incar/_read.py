import re

from ..parser import KeyValueParser

__all__ = ['IncarParser']

class IncarParser(KeyValueParser):
    '''
    parses INCAR files
    '''

    def __init__(self, filename):
        self.result = {}
        with open(filename) as incar:
            self.result = IncarParser.parse_incar(incar)

    @classmethod
    def parse_incar(cls, fobj_or_str):
        if isinstance(fobj_or_str, str):
            content = fobj_or_str
        else:
            content = fobj_or_str.read()
        kvd = dict(re.findall(cls.assignment, content))
        kvd = {k.lower(): v for k, v in kvd.iteritems()}
        return kvd
