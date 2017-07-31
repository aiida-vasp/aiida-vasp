try:
    from collections import ChainMap
except ImportError:
    from chainmap import ChainMap

from base import BasicCalculation
from aiida.orm import DataFactory


class VaspCalculation(BasicCalculation):
    '''
    General-purpose VASP calculation. By default retrieves only the 'OUTCAR', 'vasprun.xml', 'EIGENVAL', 'DOSCAR' and Wannier90 input / output files, but additional retrieve files can be specified via the 'settings['ADDITIONAL_RETRIEVE_LIST']' input.
    '''
    default_parser = 'vasp.vasp'
    charge_density = Input(types='vasp.chargedensity',
                           doc='chargedensity node: should be obtained \n' +
                           'from the output of a selfconsistent ' +
                           'NscfCalculation (written to CHGCAR)')
    wavefunctions = Input(types='vasp.wavefun',
                          doc='wavefunction node: to speed up convergence ' +
                          'for continuation jobs')

    def _prepare_for_submission(self, tempfolder, inputdict):
        '''add EIGENVAL, DOSCAR, and all files starting with wannier90 to
        the list of files to be retrieved.'''
        calcinfo = super(VaspCalculation, self)._prepare_for_submission(
            tempfolder, inputdict)
        calcinfo.retrieve_list.extend(['EIGENVAL', 'DOSCAR'])
        calcinfo.retrieve_list.append(('wannier90*', '.', 0))
        calcinfo.retrieve_list = list(set(calcinfo.retrieve_list))
        return calcinfo

    def verify_inputs(self, inputdict, *args, **kwargs):
        super(VaspCalculation, self).verify_inputs(inputdict, *args, **kwargs)
        self.check_input(inputdict, 'kpoints', self._need_kp)
        self.check_input(inputdict, 'charge_density', self._need_chgd)
        self.check_input(inputdict, 'wavefunctions', self._need_wfn)

    def _need_kp(self):
        '''
        return wether an input kpoints node is needed or not.
        :return output:
            True if input kpoints node is needed
            (py:method::VaspCalculation.use_kpoints),
            False otherwise
        needs 'parameters' input to be set
        (py:method::VaspCalculation.use_parameters)
        '''
        if 'kspacing' in self._parameters and 'kgamma' in self._parameters:
            return False
        else:
            return True

    def _need_chgd(self):
        '''
        Test wether an charge_densities input is needed or not.
        :return output:
            True if a chgcar file must be used
            (py:method::NscfCalculation.use_charge_densities),
            False otherwise
        needs 'parameters' input to be set
        (py:method::NscfCalculation.use_parameters)
        '''
        ichrg_d = self._need_wfn() and 0 or 2
        icharg = self._parameters.get('icharg', ichrg_d)
        if icharg in [1, 11]:
            return True
        else:
            return False

    def _need_wfn(self):
        '''
        Test wether a wavefunctions input is needed or not.
        :return output:
            True if a wavecar file must be
            used (py:method::NscfCalculation.use_wavefunctions),
            False otherwise
        needs 'parameters' input to be set
        (py:method::NscfCalculation.use_parameters)
        '''
        istrt_d = self.get_inputs_dict().get('wavefunctions') and 1 or 0
        istart = self._parameters.get('istart', istrt_d)
        if istart in [1, 2, 3]:
            return True
        else:
            return False

    def write_additional(self, tempfolder, inputdict):
        '''write CHGAR and WAVECAR files if needed'''
        super(VaspCalculation, self).write_additional(
            tempfolder, inputdict)
        if self._need_chgd():
            chgcar = tempfolder.get_abs_path('CHGCAR')
            self.write_chgcar(inputdict, chgcar)
        if self._need_wfn():
            wavecar = tempfolder.get_abs_path('WAVECAR')
            self.write_wavecar(inputdict, wavecar)

    def write_chgcar(self, inputdict, dst):
        import shutil
        shutil.copyfile(self.inp.charge_density.get_file_abs_path(), dst)

    def write_wavecar(self, inputdict, dst):
        import shutil
        shutil.copyfile(self.inp.wavefunctions.get_file_abs_path(), dst)
