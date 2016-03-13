from aiida.orm import DataFactory
from os.path import dirname, realpath, join


class Common(object):
    @classmethod
    def structure(cls):
        larray = np.array([[0, .5, .5],
                           [.5, 0, .5],
                           [.5, .5, 0]])
        alat = 6.058
        self.structure = DataFactory('structure')(cell=larray*alat)
        self.structure.append_atom(position=[0, 0, 0], symbols='In')
        self.structure.append_atom(position=[.25, .25, .25], symbols='As')

        return structure

    @classmethod
    def cif(cls):
        cifpath = realpath(join(dirname(__file__),
                                'data', 'EntryWithCollCode43360.cif'))
        self.cif = DataFactory('cif').get_or_create(cifpath)[0]
        reurn cif

    @classmethod
    def import_paw(cls):
        DataFactory('vasp.paw').import_family(
            realpath(join(dirname(__file__), 'LDA')), familyname='TEST')

    @classmethod
    def kpoints_mesh(cls):
        kp = DataFactory('array.kpoints')()
        kp.set_kpoints_mesh([4,4,4])

    @classmethod
    def kpoints_list(cls):
        kp = DataFactory('array.kpoints')()
        kp.set_kpoints([[0., 0., 0.], [0., 0., .5]], weights=[1., 1.])

    @classmethod
    def settings(cls):
        settings = {
            'gga': 'PE',
            'gga_compat': False,
            'lorbit': 11,
            'sigma': .05
        }
        return DataFactory('parameter')(dict=settings)
