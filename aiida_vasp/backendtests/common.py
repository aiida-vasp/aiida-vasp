"""Common tasks used in setting up tests"""
import os
from os.path import dirname, realpath, join

import numpy as np
from aiida.orm import DataFactory


def subpath(*args):
    return realpath(join(dirname(__file__), *args))


class Common(object):
    """Class encapsulating some common functionality"""

    @staticmethod
    def structure():
        """Create a StructureData object which can be used in tests"""
        larray = np.array([[0, .5, .5], [.5, 0, .5], [.5, .5, 0]])
        alat = 6.058
        structure = DataFactory('structure')(cell=larray * alat)
        structure.append_atom(position=[0, 0, 0], symbols='In')
        structure.append_atom(position=[.25, .25, .25], symbols='As')

        return structure

    @staticmethod
    def cif():
        cifpath = realpath(join(dirname(__file__), 'data', 'EntryWithCollCode43360.cif'))
        cif = DataFactory('cif').get_or_create(cifpath)[0]
        return cif

    @staticmethod
    def import_paw():
        DataFactory('vasp.paw').import_family(subpath('LDA'), familyname='TEST', family_desc='stub ' 'Potpaw Family for testing purposes')

    @staticmethod
    def paw_in():
        return DataFactory('vasp.paw').load_paw(element='In')[0]

    @staticmethod
    def paw_as():
        return DataFactory('vasp.paw').load_paw(element='As')[0]

    @staticmethod
    def kpoints_mesh():
        kpoints = DataFactory('array.kpoints')()
        kpoints.set_kpoints_mesh([2, 2, 2])
        return kpoints

    @staticmethod
    def kpoints_mesh_res():
        res = np.array([[2., 2., 2.], [0., 0., 0.]])
        return res

    @staticmethod
    def kpoints_list():
        kpoints = DataFactory('array.kpoints')()
        kpoints.set_kpoints([[0., 0., 0.], [0., 0., .5]], weights=[1., 1.])
        return kpoints

    @staticmethod
    def kpoints_list_res():
        kres = np.array([[0., 0., 0.], [0., 0., .5]])
        wres = np.array([1., 1.])
        return (kres, wres)

    @staticmethod
    def parameters():
        """Create ParameterData object with common parameters"""
        parameters = {'gga': 'PE', 'gga_compat': False, 'lorbit': 11, 'sigma': .05}
        return DataFactory('parameter')(dict=parameters)

    @staticmethod
    def charge_density():
        return DataFactory('vasp.chargedensity')(file=subpath('data', 'CHGCAR'))

    @staticmethod
    def charge_density_res():
        with open(subpath('data', 'CHGCAR'), 'r') as chg:
            res = chg.read()
        return res

    @staticmethod
    def wavefunctions():
        return DataFactory('vasp.wavefun')(file=subpath('data', 'WAVECAR'))

    @staticmethod
    def wavefunctions_res():
        with open(subpath('data', 'WAVECAR'), 'r') as wav:
            res = wav.read()
        return res

    @staticmethod
    def win():
        return DataFactory('parameter')(dict={'Test': True})

    @staticmethod
    def win_res():
        return 'Test = T'

    @staticmethod
    def wdat():
        return DataFactory('vasp.archive')()

    @staticmethod
    def retrieved_nscf():
        """Create a test retrieve folder."""
        ret = DataFactory('folder')()
        for fname in os.listdir(subpath('data', 'retrieved_nscf', 'path')):
            ret.add_path(subpath('data', 'retrieved_nscf', 'path', fname), '')
        return ret
