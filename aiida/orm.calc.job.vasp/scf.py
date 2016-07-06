from base import BasicCalculation


class ScfCalculation(BasicCalculation):
    '''
    Limited VASP calculation, can only take kpoints grid,
    returns only CHGCAR and WAVECAR
    '''
    default_parser = 'vasp.scf'

    def _prepare_for_submission(self, tempfolder, inputdict):
        '''
        changes the retrieve_list to retrieve only
        files needed to continue with nonscf runs,
        CHGCAR, WAVECAR, IBZKPT
        '''
        calcinfo = super(
            ScfCalculation, self)._prepare_for_submission(
                tempfolder, inputdict)
        calcinfo.retrieve_list.extend(['CHGCAR', 'WAVECAR', 'IBZKPT'])
        return calcinfo

    def write_additional(self, tempfolder, inputdict):
        super(ScfCalculation, self).write_additional(
            tempfolder, inputdict)

    def verify_inputs(self, inputdict, *args, **kwargs):
        '''
        make sure k-points are given as a mesh
        '''
        super(ScfCalculation, self).verify_inputs(inputdict, *args, **kwargs)
        kpa = inputdict['kpoints'].get_attrs()
        if 'mesh' not in kpa:
            raise AttributeError(
                'ScfCalculation can only be used with kpoints mesh')
        return True

    def _init_internal_params(self):
        '''
        let the metaclass py:class:`aiida.orm.calculation.job.vasp.base.CalcMeta`
        ref CalcMeta pick up internal parameters from the class body
        and insert them
        '''
        super(ScfCalculation, self)._init_internal_params()
        self._update_internal_params()
