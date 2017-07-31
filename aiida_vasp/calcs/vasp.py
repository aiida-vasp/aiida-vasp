try:
    from collections import ChainMap
except ImportError:
    from chainmap import ChainMap

from base import VaspCalcBase
from nscf import NscfCalculation
from aiida.orm import DataFactory
from aiida.common.utils import classproperty


class VaspCalculation(NscfCalculation):
    '''
    General-purpose VASP calculation. By default retrieves only the 'OUTCAR', 'vasprun.xml', 'EIGENVAL', 'DOSCAR' and Wannier90 input / output files, but additional retrieve files can be specified via the 'settings['ADDITIONAL_RETRIEVE_LIST']' input.
    '''
    default_parser = 'vasp.vasp'

    def verify_inputs(self, inputdict, *args, **kwargs):
        # ~ notset_msg = 'input not set: %s'
        super(VaspCalculation, self).verify_inputs(inputdict, *args, **kwargs)
        self.check_input(inputdict, 'parameters')
        self.check_input(inputdict, 'structure')
        if 'elements' not in self.attrs():
            self._prestore()
        for kind in self.elements:
            self.check_input(inputdict, self._get_paw_linkname(kind))
        self.check_input(inputdict, 'kpoints', self._need_kp)

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
