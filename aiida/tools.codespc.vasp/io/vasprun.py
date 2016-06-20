__doc__ = '''
Tools for parsing vasprun.xml files
'''

try:
    from lxml.objectify import parse
except:
    from xml.etree.ElementTree import parse

import datetime as dt
import numpy as np


class VasprunParser(object):
    '''
    parse xml into objecttree, provide convenience methods
    for parsing
    '''
    def __init__(self, fname):
        super(VasprunParser, self).__init__()
        self.tree = parse(fname)

    @property
    def program(self):
        return self._i('program')

    @property
    def version(self):
        return self._i('version').strip()

    @property
    def datetime(self):
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
        return ibrion == -1

    @property
    def is_md(self):
        ibrion = self.param('IBRION', default=-1)
        return ibrion not in [-1, 1, 2]

    @property
    def is_relaxation(self):
        ibrion = self.param('IBRION', default=-1)
        return ibrion in [1, 2]

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
    def tdos(self):
        return self._array(parent='dos/total')

    @property
    def pdos(self):
        try:
            dos = self._array(parent='dos/partial')
        except:
            dos = np.array([])
        return dos

    def param(self, key, default=None):
        path = '/parameters//'
        return self._i(key, path=path) or self._v(key, path=path) or default

    def _varray(self, key, path='//'):
        tag = self.tag('varray', key, path)
        if not tag:
            return None

        def split(s):
            return s.text.split()
        return np.array(map(split, tag.v), dtype=float)

    def _array(self, parent, key=None, path='//'):
        pred = key and '[@name="%s"]' % key or ''
        tag = self.tree.find(path+parent+'/array%s' % pred)
        dims = [i.text for i in tag.findall('dimension')]

        def getdtf(field):
            dt = field.attrib.get('type', float)
            if dt == 'string':
                dt = 'S128'
            return (field.text.strip(), dt)

        dtyp = np.dtype([getdtf(f) for f in tag.findall('field')])
        ndim = len(dims)
        shape = []
        subset = tag.find('set')
        for d in range(ndim-1):
            if subset.find('set'):
                shape.append(len(subset.findall('set')))
                subset = subset.find('set')
        ldim = subset.findall('r')
        mode = 'r'
        if not ldim:
            ldim = subset.findall('rc')
            mode = 'rc'
        shape.append(len(ldim))

        def split(s):
            if mode == 'r':
                return tuple(s.text.split())
            elif mode == 'rc':
                return tuple(map(lambda x: x.text.strip(), s.c))

        data = np.array(map(split, tag.iterfind('*//%s' % mode)), dtype=dtyp)
        return data.reshape(shape)

    def _i(self, key, path='//'):
        tag = self.tag('i', key, path)
        res = None
        if tag:
            if tag.attrib.get('type') == 'logical':
                res = 'T' in tag.text
            elif isinstance(tag.pyval, str):
                res = tag.text.strip()
            else:
                res = tag.pyval
        return res

    def _v(self, key, path='//'):
        tag = self.tag('v', key, path)
        if tag:
            dtype = tag.attrib.get('type', float)
            return np.array(tag.text.split(), dtype=dtype)
        else:
            return None

    def _fppath(self):
        return '//structure[@name="finalpos"]//'

    def tag(self, tag, key, path='//'):
        path = '{p}{t}[@name="{n}"]'.format(p=path, t=tag, n=key)
        return self.tree.find(path)
