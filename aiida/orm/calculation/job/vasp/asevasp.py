#encoding: utf-8
from aiida.orm import JobCalculation, DataFactory
#~ from aiida.orm.data.parameter import ParameterData, SinglefileData
from aiida.common.datastructures import CalcInfo
from aiida.common.utils import classproperty
#~ from aiida.orm.data.vasp import VaspStructureData, VaspPotentialData, VaspKPointData
import numpy as np
import os
import pymatgen as pmg
from aiida.orm.calculation.job.vasp.incar import _incarify

ParameterData = DataFactory('parameter')
StructureData = DataFactory('structure')
SinglefileData = DataFactory('singlefile')

__copyright__ = u'Copyright (c), 2015, Rico Häuselmann'

class AsevaspCalculation(JobCalculation):
    '''
    Perform a calculation with VASP
    '''
    def _init_internal_params(self):
        super(AsevaspCalculation, self)._init_internal_params()
        self._default_parser = 'vasp.asevasp'
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
                'valid_types': StructureData,
                'additional_parameter': None,
                'linkname': 'structure_in',
                'docstring': 'Structure Data as in POSCAR file',
            },
            'potentials': {
                'valid_types': (SinglefileData, ParameterData),
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
            'chgcar': {
                'valid_types': (None, SinglefileData),
                'additional_parameter': None,
                'linkname': 'chgcar_in,
                'docstring': 'CHGCAR file from a selfconsistent run (only used in non-SC runs)'
            }
        })
        return retdict

    def _prepare_for_submission(self, tempfolder, inputdict):
        try:
            import ase
        except:
            print 'ase needs to be installed'

        # ---------- get input data ----------
        incar_data = inputdict['settings']
        poscar_data = inputdict['structure_in']
        potcar_data = inputdict['potentials_in']
        kpoints_data = inputdict['kpoints_in']
        chgcar_data = inputdict['chgcar_in']

        # ---------- set up an ase Vasp calculation ----------
        from ase.calculators.vasp import Vasp as AseVasp
        asevasp = AseVasp(xc='PBE')
        # push incar info
        for k,v in incar_data.get_dict().iteritems():
            try:
                asevasp.set(**{k.lower() : v})
            except TypeError as te:
                print 'ASE: '+te.message
                asevasp.string_params[k.lower()] = _incarify(v)

        # push kpoints info
        kpoints_dict = kpoints_data.get_dict()
        if kpoints_dict.get('generation_style') == 'Gamma':
            asevasp.set(gamma=True)
        kp = np.array(kpoints_dict['kpoints'])
        if len(kp.flatten()) == 3:
            kp = kp.flatten()
        asevasp.set(kpts=kp)

        # check wether potcar is parameter
        if isinstance(potcar_data, ParameterData):
            potcar_dict = potcar_data.get_dict()
            if potcar_dict.get('potpaw_path'):
                os.environ['VASP_PP_PATH'] = potcar_dict['potpaw_path']
            else:
                if not os.environ['VASP_PP_PATH']:
                    raise InputValidationError(
                        'No VASP_PP_PATH env variable is set. '+\
                        'Please provide a potpaw_path parameter with the path containing '+\
                        'your potpaw_XXX directory in the ParameterData passed to '+\
                        'use_potentials or export VASP_PP_PATH=<path containing potpaw> '+\
                        'or consider passing a POTCAR in a SinglefileData instance.')
            asevasp.set(setups=potcar_dict['special_symbols'])

        # push structure info
        asevasp.set_atoms(poscar_data.get_ase())

        asevasp.initialize(asevasp.atoms)

        # ---------- write input files ----------
        os.chdir(tempfolder.get_abs_path(''))
        from ase.io import vasp as vio
        vio.write_vasp('POSCAR', asevasp.atoms, direct=True, vasp5=True)
        asevasp.write_incar(asevasp.atoms)
        asevasp.write_kpoints()

        if isinstance(potcar_data, ParameterData):
            asevasp.write_potcar()
        else:
            import shutil
            potcar_src = potcar_data.get_file_abs_path()
            potcar_dst = tempfolder.get_abs_path('POTCAR')
            shutil.copy(potcar_src, potcar_dst)
            #~ with open(potcar_data.get_file_abs_path()) as potcar_in:
                #~ with open(potcar_file, 'w') as potcar_out:
                    #~ potcar_out.write(potcar_in.read())

        icharg = incar_data.get_dict().get('icharg')
        if not icharg:
            icharg = incar_data.get_dict().get('ICHARG')
        if icharg < 10:
            if chgcar_data:
                self.logger.warn('ICHARG tag in INCAR file is {}, input CHGCAR not used!'.format(icharg))
        elif icharg >= 10:
            if not chgcar_data:
                raise InputValidationError('ICHARG tag in INCAR file is {}, you must give a CHGCAR file!'.format(icharg))
            else:
                import shutil
                chgcar_src = chgcar_data.get_file_abs_path()
                chgcar_dst = tempfolder.get_abs_path('CHGCAR')
                shutil.copy(chgcar_src, chgcar_dst)
                #~ with open(chgcar_data.get_file_abs_path()) as chgcar_in:
                    #~ with open(chgcar_file, 'w') as chgcar_out:
                        #~ chgcar_out.write(chgcar_in.read())

        # ---------- set up and return calcinfo ----------
        calcinfo = CalcInfo()
        calcinfo.uuid = self.uuid
        calcinfo.local_copy_list = []
        calcinfo.remote_copy_list = []
        #~ calcinfo.stdin_name = self._INPUT_FILE_NAME
        #~ calcinfo.stdout_name = self._OUTPUT_FILE_NAME
        calcinfo.retrieve_list = ['INCAR', 'KPOINTS', 'POSCAR', 'POTCAR', 'CHG', 'CHGCAR', 'CONTCAR', 'DOSCAR',
                                  'EIGENVAL', 'OSZICAR', 'OUTCAR', 'PCDAT',
                                  'PROCAR', 'WAVECAR', 'XDATCAR', 'vasprun.xml']
        return calcinfo

