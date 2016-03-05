from base import VaspCalcBase, Input
from aiida.common.utils import classproperty
from aiida.orm import DataFactory


class ScfCalculation(VaspCalcBase):
    '''
    Limited VASP calculation, can only take kpoints grid,
    returns only CHGCAR and WAVECAR
    '''
    settings = Input(types='parameter')
    structure = Input(types=['structure', 'cif'])
    paw = Input(types='vasp.paw', param='kind')
    kpoints = Input(types='array.kpoints')
    default_parser = 'vasp.scf'

    def _prepare_for_submission(self, tempfolder, inputdict):
        '''changes the retrieve_list to retrieve only files needed to continue with nonscf runs'''
        calcinfo = super(ScfCalculation, self)._prepare_for_submission(tempfolder, inputdict)
        calcinfo.retrieve_list = ['CHGCAR', 'WAVECAR', 'IBZKPT', 'OUTCAR']
        return calcinfo

    def write_incar(self, inputdict, dst):
        '''
        converts from settings node (ParameterData) to INCAR format
        and writes to dst
        :param
            inputdict: required by baseclass
            dst: absolute path of the file to write to
        '''
        from incar import dict_to_incar
        with open(dst, 'w') as incar:
            incar.write(dict_to_incar(self.inp.settings.get_dict()))

    def write_poscar(self, inputdict, dst):
        '''
        converts from structures node (StructureData) to POSCAR format
        and writes to dst
        :param
            inputdict: required by baseclass
            dst: absolute path of the file to write to
        '''
        from ase.io.vasp import write_vasp
        with open(dst, 'w') as poscar:
            write_vasp(poscar, self.inp.structure.get_ase(), vasp5=True)

    def write_potcar(self, inputdict, dst):
        '''
        concatenatest multiple paw files into a POTCAR
        :param
            inputdict: required by baseclass
            dst: absolute path of the file to write to
        '''
        import subprocess32 as sp
        catcom = ['cat']
        # ~ structure = inputdict['structure']
        # ~ structure = self.inp.structure
        # order the symbols according to order given in structure
        if 'elements' not in self.attrs():
            self._prestore()
        for kind in self.elements:
            paw = inputdict[self._get_paw_linkname(kind)]
            catcom.append(paw.get_abs_path('POTCAR'))
        # cat the pawdata nodes into the file
        with open(dst, 'w') as pc:
            sp.check_call(catcom, stdout=pc)

    def write_kpoints(self, inputdict, dst):
        '''
        converts from kpoints node (KpointsData) to KPOINTS format
        and writes to dst
        :param
            inputdict: required by baseclass
            dst: absolute path of the file to write to
        '''
        from base import kpmtemp, kplitemp, kpltemp
        kp = self.inp.kpoints
        mesh, offset = kp.get_kpoints_mesh()
        with open(dst, 'w') as kpoints:
            kps = kpmtemp.format(N=mesh, s=offset)
            kpoints.write(kps)

    def write_additional(self, tempfolder, inputdict):
        super(ScfCalculation, self).write_additional(
            tempfolder, inputdict)

    def verify_inputs(self, inputdict, *args, **kwargs):
        # ~ notset_msg = 'input not set: %s'
        super(ScfCalculation, self).verify_inputs(inputdict, *args, **kwargs)
        self.check_input(inputdict, 'settings')
        self.check_input(inputdict, 'structure')
        if 'elements' not in self.attrs():
            self._prestore()
        for kind in self.elements:
            self.check_input(inputdict, self._get_paw_linkname(kind))
        self.check_input(inputdict, 'kpoints')
        try:
            inputdict['kpoints'].get_kpoints_mesh()
        except AttributeError as e:
            print('ScfCalculation can only be used with kpoints mesh')
            raise e

    @classmethod
    def _get_paw_linkname(cls, kind):
        return 'paw_%s' % kind

    def _prestore(self):
        '''
        set attributes prior to storing
        '''
        super(ScfCalculation, self)._prestore()
        self._set_attr('elements', list(set(
            self.inp.structure.get_ase().get_chemical_symbols())))

    @property
    def _settings(self):
        return {k.lower(): v for
                k, v in self.inp.settings.get_dict().iteritems()}

    @classmethod
    def new_settings(self, **kwargs):
        return DataFactory('parameter')(**kwargs)

    @classmethod
    def new_structure(self, **kwargs):
        return DataFactory('structure')(**kwargs)

    @classmethod
    def new_kpoints(self, **kwargs):
        return DataFactory('array.kpoints')(**kwargs)

    @classmethod
    def load_paw(self, *args, **kwargs):
        return self.Paw.load_paw(*args, **kwargs)[0]

    @classproperty
    def Paw(self):
        return DataFactory('vasp.paw')

    @property
    def elements(self):
        return self.get_attr('elements')

    def _init_internal_params(self):
        '''
        let the metaclass py:class:`~aiida.orm.calculation.job.vasp.base.CalcMeta` ref CalcMeta pick up internal parameters from the class body
        and insert them
        '''
        super(ScfCalculation, self)._init_internal_params()
        self._update_internal_params()
