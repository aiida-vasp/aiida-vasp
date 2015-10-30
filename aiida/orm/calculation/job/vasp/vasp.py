#encoding: utf-8
from aiida.orm import JobCalculation, DataFactory
#~ from aiida.orm.data.parameter import ParameterData, SinglefileData
from aiida.common.datastructures import CalcInfo
from aiida.common.utils import classproperty
#~ from aiida.orm.data.vasp import VaspStructureData, VaspPotentialData, VaspKPointData

ParameterData = DataFactory('parameter')
SinglefileData = DataFactory('singlefile')

__copyright__ = u'Copyright (c), 2015, Rico Häuselmann'

class VaspCalculation(JobCalculation):
    '''
    Perform a calculation with VASP
    '''
    def _init_internal_params(self):
        super(VaspCalculation, self)._init_internal_params()
        self._default_parser = 'vasp.vasp'
        self._INPUT_FILE_NAME = 'INCAR'
        self._OUTPUT_FILE_NAME = 'OUTCAR'
        self._DEFAULT_OUTPUT_FILE = self._OUTPUT_FILE_NAME
        self._DEFAULT_INPUT_FILE = self._INPUT_FILE_NAME

    @classproperty
    def _use_methods(cls):
        retdict = JobCalculation._use_methods
        retdict.update({
            'settings': {
                'valid_types': ParameterData,
                'additional_parameter': None,
                'linkname': 'settings',
                'docstring': 'Settings that go into the INCAR file for vasp',
            },
            'structure': {
                'valid_types': ParameterData,
                'additional_parameter': None,
                'linkname': 'structure_in',
                'docstring': 'Structure Data as in POSCAR file',
            },
            'potentials': {
                'valid_types': SinglefileData,
                'additional_parameter': None,
                'linkname': 'potentials_in',
                'docstring': 'Potentials as in POTCAR file',
            },
            'kpoints': {
                'valid_types': ParameterData,
                'additional_parameter': None,
                'linkname': 'kpoints_in',
                'docstring': 'k points as in KPOINTS file',
            },
        })
        return retdict

    def _prepare_for_submission(self, tempfolder, inputdict):
        try:
            import pymatgen as pmg
        except:
            print 'import pymatgen as pmg failed'

        incar_file = tempfolder.get_abs_path('INCAR')
        incar_data = inputdict['settings']
        incar_object = pmg.io.vasp.Incar.from_dict(incar_data.get_dict())
        incar_object.write_file(incar_file)

        poscar_file = tempfolder.get_abs_path('POSCAR')
        poscar_data = inputdict['structure_in']
        poscar_object = pmg.io.vasp.Poscar.from_dict(poscar_data.get_dict())
        poscar_object.write_file(poscar_file)

        potcar_file = tempfolder.get_abs_path('POTCAR')
        potcar_data = inputdict['potentials_in'].get_file_abs_path()
        #~ potcar_object = pmg.io.vasp.Potcar.from_dict(potcar_data.get_dict())
        #~ potcar_object.write_file(potcar_file)
        with open(potcar_data) as potcar_in:
            with open(potcar_file, 'w') as potcar_out:
                potcar_out.write(potcar_in.read())

        kpoints_file = tempfolder.get_abs_path('KPOINTS')
        kpoints_data = inputdict['kpoints_in']
        kpoints_object = pmg.io.vasp.Kpoints.from_dict(kpoints_data.get_dict())
        kpoints_object.write_file(kpoints_file)

        calcinfo = CalcInfo()
        calcinfo.uuid = self.uuid
        calcinfo.local_copy_list = []
        calcinfo.remote_copy_list = []
        #~ calcinfo.stdin_name = self._INPUT_FILE_NAME
        #~ calcinfo.stdout_name = self._OUTPUT_FILE_NAME
        calcinfo.retrieve_list = ['CHG', 'CHGCAR', 'CONTCAR', 'DOSCAR',
                                  'EIGENVAL', 'OSZICAR', 'OUTCAR', 'PCDAT',
                                  'PROCAR', 'WAVECAR', 'XDATCAR', 'vasprun.xml']
        return calcinfo
