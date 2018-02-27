"""
Tools for parsing vasprun.xml files
"""

try:
    from lxml.objectify import parse
except ImportError:
    from xml.etree.ElementTree import parse
import datetime as dt
import numpy as np

from aiida_vasp.io.parser import BaseParser

DEFAULT_OPTIONS = {'quantities_to_parse': ['occupations', 'vrp_pdos', 'vrp_tdos']}


class VasprunParser(BaseParser):
    """
    parse xml into objecttree, provide convenience methods
    for parsing
    """

    PARSABLE_ITEMS = {
        'occupations': {
            'inputs': [],
            'parsers': ['vasprun.xml'],
            'nodeName': 'intermediate_data',
            'prerequisites': []
        },
        'vrp_pdos': {
            'inputs': [],
            'parsers': ['vasprun.xml'],
            'nodeName': 'intermediate_data',
            'prerequisites': []
        },
        'vrp_tdos': {
            'inputs': [],
            'parsers': ['vasprun.xml'],
            'nodeName': 'intermediate_data',
            'prerequisites': []
        },
    }

    def __init__(self, path, filename):
        super(VasprunParser, self).__init__()
        self._filepath = path
        self._filename = filename
        self._parsed_data = None
        self._parsable_items = VasprunParser.PARSABLE_ITEMS

        self.tree = parse(filename)

    def _parse_file(self, inputs):

        settings = inputs.gets('settings', DEFAULT_OPTIONS)

        result = {}
        for quantity in settings['quantities_to_parse']:
            if quantity in self._parsable_items:
                result[quantity] = getattr(self, quantity)()

        return result

    @property
    def program(self):
        return self._i('program')

    @property
    def version(self):
        return self._i('version').strip()

    @property
    def datetime(self):
        """Parse Date and time information into a Python datetime object."""
        date = self._i('date')
        time = self._i('time')
        dtstr = date + ' ' + time
        return dt.datetime.strptime(dtstr, '%Y %m %d %H:%M:%S')

    @property
    def cell(self):
        return self._varray('basis', path=self._fppath())

    @property
    def volume(self):
        return self._i('basis', path=self._fppath())

    @property
    def pos(self):
        return self._varray('positions', path=self._fppath())

    @property
    def root(self):
        return self.tree.getroot()

    @property
    def efermi(self):
        return self._i('efermi')

    @property
    def is_static(self):
        ibrion = self.param('IBRION', default=-1)
        nsw = self.param('NSW', default=0)
        return (ibrion == -1) or (nsw == 0)

    @property
    def is_md(self):
        ibrion = self.param('IBRION', default=-1)
        return ibrion == 0

    @property
    def is_relaxation(self):
        ibrion = self.param('IBRION', default=-1)
        nsw = self.param('NSW', default=0)
        return (ibrion in [1, 2, 3]) and (nsw > 0)

    @property
    def is_sc(self):
        icharg = self._i('ICHARG')
        return icharg < 10

    @property
    def occupations(self):
        eig = self._array(parent='calculation/eigenvalues')
        return eig['occ']

    @property
    def projected_occupations(self):
        eig = self._array(parent='calculation/projected/eigenvalues')
        return eig['occ']

    @property
    def bands(self):
        eig = self._array(parent='calculation/eigenvalues')
        return eig['eigene']

    @property
    def projected_bands(self):
        eig = self._array(parent='calculation/projected/eigenvalues')
        return eig['eigene']

    @property
    def vrp_tdos(self):
        return self._array(parent='dos/total')

    @property
    def vrp_pdos(self):
        """The partial DOS array"""
        try:
            dos = self._array(parent='dos/partial')
        except Exception:  # pylint: disable=broad-except
            dos = np.array([])
        return dos

    def param(self, key, default=None):
        path = '/parameters//'
        return self._i(key, path=path) or self._v(key, path=path) or default

    def _varray(self, key, path='//'):
        """Extract a <varray> tag"""
        tag = self.tag('varray', key, path)
        if tag is None:
            return None

        def split(string_):
            return string_.text.split()

        return np.array(map(split, tag.v), dtype=float)

    def _array(self, parent, key=None, path='//'):
        """extract an <array> tag"""
        pred = '[@name="%s"]' % key if key else ''
        tag = self.tree.find(path + parent + '/array%s' % pred)
        dims = [i.text for i in tag.findall('dimension')]

        def getdtf(field):
            d_type = field.attrib.get('type', float)
            if d_type == 'string':
                d_type = 'S128'
            return (field.text.strip(), d_type)

        dtyp = np.dtype([getdtf(f) for f in tag.findall('field')])
        ndim = len(dims)
        shape = []
        subset = tag.find('set')
        for _ in range(ndim - 1):
            if subset.find('set') is not None:
                shape.append(len(subset.findall('set')))
                subset = subset.find('set')
        ldim = subset.findall('r')
        mode = 'r'
        if not ldim:
            ldim = subset.findall('rc')
            mode = 'rc'
        shape.append(len(ldim))

        def split(string_):
            """Splits a string based on mode in ['r', 'rc']"""
            if mode == 'r':
                return tuple(string_.text.split())
            elif mode == 'rc':
                return tuple([x.text.strip() for x in string_.c])
            return None

        data = np.array(map(split, tag.iterfind('*//%s' % mode)), dtype=dtyp)
        return data.reshape(shape)

    def _i(self, key, path='//'):
        """Extract an <i> tag"""
        tag = self.tag('i', key, path)
        res = None
        if tag is not None:
            if tag.attrib.get('type') == 'logical':
                res = 'T' in tag.text
            elif tag.attrib.get('type') == 'int':
                res = int(tag.text)
            elif tag.attrib.get('type') == 'string':
                res = tag.text.strip()
            else:
                try:
                    res = int(tag.text)
                except ValueError:
                    try:
                        res = float(tag.text)
                    except ValueError:
                        res = tag.text.strip()
        return res

    def _v(self, key, path='//'):
        """Extract an <v> tag"""
        tag = self.tag('v', key, path)
        if tag is not None:
            dtype = tag.attrib.get('type', float)
            return np.array(tag.text.split(), dtype=dtype)
        return None

    @staticmethod
    def _fppath():
        return '//structure[@name="finalpos"]//'

    def tag(self, tag, key, path='//'):
        path = '{p}{t}[@name="{n}"]'.format(p=path, t=tag, n=key)
        return self.tree.find(path)
