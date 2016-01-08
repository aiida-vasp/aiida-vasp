__doc__ = '''
This module contains tools to read PotpawData attributes from POTCAR files.
'''
import re
import datetime as dt
from parser import BaseParser


class PotParser(KeyValueParser):
    '''
    contains regex and functions to find grammar elements
    for POTCAR files
    '''
    assignment = re.compile(r'(\w*)\s*=\s*([^;]*);?')
    comments = True

    @classmethod
    def single(cls, kv_list):
        k_list = [ssl[0] for subl in kv_list for ssl in subl]
        if len(filter(lambda x: x == 'TITEL', k_list)) == 1:
            return True
        else:
            return False

    @classmethod
    def title(cls, title):
        tl = title.split()
        fam = tl[0].replace('PAW_', '')
        sym = tl[1]
        date = dt.datetime.strptime(tl[2], '%d%b%Y').date()
        return fam, sym, date

    @classmethod
    def vrhfin(cls, vrhfin):
        vl = vrhfin.split()
        element = vl[0].strip(':')
        spconf = vl[1]
        return element, spconf

    @classmethod
    def parse_potcar(cls, filename):
        kv_list = cls.kv_list(filename)
        if not cls.single(kv_list):
            raise ValueError('not parsing concatenated POTCAR files')
        kv_dict = cls.kv_dict(kv_list)
        is_paw, cmt = cls.bool(kv_dict['LPAW'])
        is_ultrasoft, cmt = cls.bool(kv_dict['LULTRA'])
        if not is_paw:
            raise ValueError('POTCAR contains non-PAW potential')
        if is_ultrasoft:
            raise ValueError('POTCAR contains ultrasoft PP')

        attr_dict = {}
        fam, sym, date = cls.title(kv_dict['TITEL'])
        attr_dict['family'] = fam
        attr_dict['symbol'] = sym
        attr_dict['paw_date'] = date
        enmin, enmin_u, cmt = cls.float_unit(kv_dict['ENMIN'])
        attr_dict['enmin'] = enmin
        attr_dict['enmin_unit'] = enmin_u
        enmax, enmax_u, cmt = cls.float_unit(kv_dict['ENMAX'])
        attr_dict['enmax'] = enmax
        attr_dict['enmax_unit'] = enmax_u
        elem, spconf = cls.vrhfin(kv_dict['VRHFIN'])
        attr_dict['element'] = elem
        attr_dict['atomic_conf'] = spconf
        mass, cmt = cls.float(kv_dict['POMASS'])
        attr_dict['mass'] = mass
        val, cmt = cls.float(kv_dict['ZVAL'])
        attr_dict['valence'] = val
        xc_type, cmt = cls.string(kv_dict['LEXCH'])
        attr_dict['xc_type'] = xc_type

        return attr_dict
