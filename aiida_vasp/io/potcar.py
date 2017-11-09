"""
This module contains tools to read PawData attributes from POTCAR files.
"""
import re
import datetime as dt
import os

from .parser import KeyValueParser


class PotParser(KeyValueParser):
    """
    contains regex and functions to find grammar elements
    for POTCAR files
    """
    assignment = re.compile(r'(\w*)\s*=\s*([^;]*);?')
    comments = True

    @classmethod
    def single(cls, kv_list):
        k_list = [ssl[0] for subl in kv_list for ssl in subl]
        return bool([k for k in k_list if k == 'TITEL'])

    @classmethod
    def title(cls, title):
        """Parse the title line"""
        title_list = title.split()
        fam = title_list[0].replace('PAW_', '')
        sym = title_list[1]
        date = dt.datetime.strptime(title_list[2], '%d%b%Y').date()
        return fam, sym, date

    @classmethod
    def vrhfin(cls, vrhfin):
        vrhfin_list = vrhfin.split(':')
        element = vrhfin_list.pop(0).strip()
        spconf = vrhfin_list and vrhfin_list.pop(0).strip()
        return element, spconf

    @classmethod
    # pylint: disable=too-many-locals
    def parse_potcar(cls, filename):
        """Parse a VASP POTCAR file"""
        kv_list = cls.kv_list(filename)
        if not cls.single(kv_list):
            raise ValueError('not parsing concatenated POTCAR files')
        kv_dict = cls.kv_dict(kv_list)
        is_paw, _ = cls.bool(kv_dict.get('LPAW', 1))
        is_ultrasoft, _ = cls.bool(kv_dict['LULTRA'])
        if not is_paw:
            raise ValueError('POTCAR contains non-PAW potential')
        if is_ultrasoft:
            raise ValueError('POTCAR contains ultrasoft PP')

        attr_dict = {}
        fam, sym, date = cls.title(kv_dict['TITEL'])
        attr_dict['family'] = fam
        attr_dict['symbol'] = sym
        attr_dict['paw_date'] = date
        enmin, enmin_u, _ = cls.float_unit(kv_dict['ENMIN'])
        attr_dict['enmin'] = enmin
        attr_dict['enmin_unit'] = enmin_u
        enmax, enmax_u, _ = cls.float_unit(kv_dict['ENMAX'])
        attr_dict['enmax'] = enmax
        attr_dict['enmax_unit'] = enmax_u
        elem, spconf = cls.vrhfin(kv_dict['VRHFIN'])
        attr_dict['element'] = elem
        attr_dict['atomic_conf'] = spconf
        mass, _ = cls.float(kv_dict['POMASS'])
        attr_dict['mass'] = mass
        val, _ = cls.float(kv_dict['ZVAL'])
        attr_dict['valence'] = val
        xc_type, _ = cls.string(kv_dict['LEXCH'])
        attr_dict['xc_type'] = xc_type

        return attr_dict


class PawParser(KeyValueParser):
    """
    contains regex and functions to find grammar elements
    in POTCAR files in PAW libraries
    """
    assignment = re.compile(r'(\w*)\s*=\s*([^;]*);?')
    comments = True

    @classmethod
    def single(cls, kv_list):
        k_list = [ssl[0] for subl in kv_list for ssl in subl]
        return bool(len([k for k in k_list if k == 'TITEL']) == 1)

    @classmethod
    def title(cls, title):
        """Parse the title line"""
        title_list = title.split()
        fam = title_list.pop(0).replace('PAW_',
                                        '') if len(title_list) > 1 else 'none'
        fam = fam.replace('PAW', 'none')
        sym = title_list.pop(0) if title_list else ''
        date = title_list.pop(0) if title_list else 'none'
        return fam, sym, date

    @classmethod
    def vrhfin(cls, vrhfin):
        vrhfin_list = vrhfin.split(':')
        element = vrhfin_list.pop(0).strip() if vrhfin_list else ''
        spconf = vrhfin_list.pop(0).strip() if vrhfin_list else ''
        return element, spconf

    @classmethod
    # pylint: disable=too-many-locals
    def parse_potcar(cls, filename):
        """Parse a Vasp POTCAR file"""
        kv_list = cls.kv_list(filename)
        if not cls.single(kv_list):
            raise ValueError('not parsing concatenated POTCAR files')
        kv_dict = cls.kv_dict(kv_list)
        is_paw, _ = cls.bool(kv_dict.get('LPAW', 'T'))
        print kv_dict.get('LULTRA', 'F')
        is_ultrasoft, _ = cls.bool(kv_dict.get('LULTRA', 'F'))
        if not is_paw:
            raise ValueError('POTCAR contains non-PAW potential')
        if is_ultrasoft:
            raise ValueError('POTCAR contains ultrasoft PP')

        attr_dict = {}
        try:
            fam, sym, date = cls.title(kv_dict['TITEL'])
            attr_dict['family'] = fam
            attr_dict['symbol'] = sym
            attr_dict['paw_date'] = date
            enmin, enmin_u, _ = cls.float_unit(kv_dict['ENMIN'])
            attr_dict['enmin'] = enmin
            attr_dict['enmin_unit'] = enmin_u
            enmax, enmax_u, _ = cls.float_unit(kv_dict['ENMAX'])
            attr_dict['enmax'] = enmax
            attr_dict['enmax_unit'] = enmax_u
            elem, spconf = cls.vrhfin(kv_dict['VRHFIN'])
            attr_dict['element'] = elem
            attr_dict['atomic_conf'] = spconf
            mass, _ = cls.float(kv_dict['POMASS'])
            attr_dict['mass'] = mass
            val, _ = cls.float(kv_dict['ZVAL'])
            attr_dict['valence'] = val
            xc_type, _ = cls.string(kv_dict['LEXCH'])
            attr_dict['xc_type'] = xc_type
        except KeyError as err:
            msg = 'missing or misspelled keyword "' + err.message
            msg += '" in ' + os.path.abspath(filename)
            raise KeyError(msg)
        except Exception:
            import sys
            err = sys.exc_info()[1]
            msg = err.message
            msg += ' in file: ' + os.path.abspath(filename)
            raise err.__class__(msg)

        return attr_dict
