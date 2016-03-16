from aiida.orm import DataFactory
import os
from os.path import dirname, realpath, join
import numpy as np


def subpath(*args):
    return realpath(join(dirname(__file__), *args))


class Common(object):
    @classmethod
    def structure(cls):
        larray = np.array([[0, .5, .5],
                           [.5, 0, .5],
                           [.5, .5, 0]])
        alat = 6.058
        structure = DataFactory('structure')(cell=larray*alat)
        structure.append_atom(position=[0, 0, 0], symbols='In')
        structure.append_atom(position=[.25, .25, .25], symbols='As')

        return structure

    @classmethod
    def cif(cls):
        cifpath = realpath(join(dirname(__file__),
                                'data', 'EntryWithCollCode43360.cif'))
        cif = DataFactory('cif').get_or_create(cifpath)[0]
        return cif

    @classmethod
    def import_paw(cls):
        DataFactory('vasp.paw').import_family(
            subpath('LDA'), familyname='TEST')

    @classmethod
    def paw_in(cls):
        return DataFactory('vasp.paw').load_paw(element='In')[0]

    @classmethod
    def paw_as(cls):
        return DataFactory('vasp.paw').load_paw(element='As')[0]

    @classmethod
    def kpoints_mesh(cls):
        kp = DataFactory('array.kpoints')()
        kp.set_kpoints_mesh([2, 2, 2])
        return kp

    @classmethod
    def kpoints_mesh_res(cls):
        res = np.array([[2., 2., 2.], [0., 0., 0.]])
        return res

    @classmethod
    def kpoints_list(cls):
        kp = DataFactory('array.kpoints')()
        kp.set_kpoints([[0., 0., 0.], [0., 0., .5]], weights=[1., 1.])
        return kp

    @classmethod
    def kpoints_list_res(cls):
        kres = np.array([[0., 0., 0.], [0., 0., .5]])
        wres = np.array([1., 1.])
        return (kres, wres)

    @classmethod
    def settings(cls):
        settings = {
            'gga': 'PE',
            'gga_compat': False,
            'lorbit': 11,
            'sigma': .05
        }
        return DataFactory('parameter')(dict=settings)

    @classmethod
    def charge_density(cls):
        return DataFactory('vasp.chargedensity')(
            file=subpath('data', 'CHGCAR')
        )

    @classmethod
    def charge_density_res(cls):
        with open(subpath('data', 'CHGCAR'), 'r') as chg:
            res = chg.read()
        return res

    @classmethod
    def wavefunctions(cls):
        return DataFactory('vasp.wavefun')(
            file=subpath('data', 'WAVECAR')
        )

    @classmethod
    def wavefunctions_res(cls):
        with open(subpath('data', 'WAVECAR'), 'r') as wav:
            res = wav.read()
        return res

    @classmethod
    def retrieved_nscf(cls):
        ret = DataFactory('folder')()
        for fname in os.listdir(subpath('data', 'retrieved_nscf', 'path')):
            ret.add_path(subpath('data', 'retrieved_nscf', 'path', fname), '')
        return ret
