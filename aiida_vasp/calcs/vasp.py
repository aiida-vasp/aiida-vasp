#encoding: utf-8
# pylint: disable=abstract-method
# explanation: pylint wrongly complains about (aiida) Node not implementing query
"""VASP - Calculation: Generic run using pymatgen_aiida for file preparation"""
try:
    from collections import ChainMap
except ImportError:
    from chainmap import ChainMap

from aiida.plugins import DataFactory

from aiida_vasp.io.incar import IncarIo
from aiida_vasp.io.potcar import MultiPotcarIo
from aiida_vasp.io.poscar import PoscarParser
from aiida_vasp.io.kpoints import KpParser
from aiida_vasp.utils.aiida_utils import get_data_node
from aiida_vasp.calcs.base import VaspCalcBase
from aiida_vasp.utils.inheritance import update_docstring

PARAMETER_CLS = DataFactory('dict')
SINGLEFILE_CLS = DataFactory('singlefile')

_IMMIGRANT_EXTRA_KWARGS = """
vasp.vasp specific kwargs:

:param use_chgcar: bool, if True, read the CHGCAR file (has to exist) and convert it to an input node.
:param use_wavecar: bool, if True, read the WAVECAR file (has to exist) and convert it to an input node.
"""


@update_docstring('immigrant', _IMMIGRANT_EXTRA_KWARGS, append=True)
class VaspCalculation(VaspCalcBase):
    """
    General-purpose VASP calculation.

    By default retrieves only the 'OUTCAR', 'vasprun.xml', 'EIGENVAL', 'DOSCAR' and Wannier90 input / output files,
    but additional retrieve files can be specified via the ``settings['ADDITIONAL_RETRIEVE_LIST']`` input.

    Floating point precision for writing POSCAR files can be adjusted using ``settings['poscar_precision']``, default: 10

    The following assumes you are familiar with the AiiDA data structures and how to set up and run an AiiDA calculation in general.

    Example usage::

        from aiida.orm import CalculationFactory, DataFactory
        from aiida.work import submit

        proc = CalculationFactory('vasp.vasp').process()
        inputs = proc.get_inputs_template()

        inputs.parameter = <Dict with INCAR params>
        inputs.structure = <StructureData> or <CifData>
        inputs.kpoints = <KpointsData>
        inputs.settings = <Dict with parser settings etc>
        inputs.potential = DataFactory('vasp.potcar').get_potcars_from_structure(structure, ...)
        inputs.code = <Code representing vasp on your cluster>
        inputs._options = <Computer, resources, etc, AiiDA specific stuff>

        submit(proc, **inputs)

    Example Low-Level usage::

        ## assuming already set up incar (Dict), structure, kpoints, settings, etc
        calc = CalculationFactory('vasp.vasp')()
        calc.use_parameter(incar)
        calc.use_structure(structure)
        calc.use_kpoints(kpoints)
        calc.use_settings(settings)
        calc.use_potential(potential_kind_1, kind=<kind 1>)
        calc.use_potential(potential_kind_2, kind=<kind 2>)
        ## unspecific to this calculation: set computer, resources, etc
        calc.submit()

    """

    _DEFAULT_PARAMETERS = {}
    _ALWAYS_RETRIEVE_LIST = ['CONTCAR', 'OUTCAR', 'vasprun.xml', 'EIGENVAL', 'DOSCAR', ('wannier90*', '.', 0)]
    _query_type_string = 'vasp.vasp'
    _plugin_type_string = 'vasp.vasp'

    @classmethod
    def define(cls, spec):
        super(VaspCalculation, cls).define(spec)
        spec.input('parameters', valid_type=get_data_class('dict'), help='The VASP input parameters (INCAR).')
        spec.input('structure', valid_type=(get_data_class('dict'), get_data_class('cif')), help='The input structure (POSCAR).')
        # Need dynamic on this as it should also accept a parameter `kind`
        spec.input_namespace('potential', valid_type=get_data_class('vasp.potcar'), help='The potentials (POTCAR).', dynamic=True)
        spec.input('kpoints', valid_type=get_data_class('array.kpoints'), help='The kpoints to use (KPOINTS).')
        spec.input('charge_density', valid_type=get_data_class('vasp.chargedensity'), required=False, help='The charge density. (CHGCAR)')
        spec.input(
            'wavefunctions', valid_type=get_data_class('vasp.wavefun'), required=False, help='The wave function coefficients. (WAVECAR)')
        spec.input(
            'settings', valid_type=get_data_class('dict'), required=False, help='Additional parameters not related to VASP itself.')

    def prepare_for_submission(self, tempfolder):
        """Add EIGENVAL, DOSCAR, and all files starting with wannier90 to the list of files to be retrieved."""
        calcinfo = super(VaspCalculation, self).prepare_for_submission(tempfolder)
        try:
            additional_retrieve_list = self.inputs.settings.get_attr('ADDITIONAL_RETRIEVE_LIST')
        except (KeyError, AttributeError):
            additional_retrieve_list = []
        calcinfo.retrieve_list = list(set(self._ALWAYS_RETRIEVE_LIST + additional_retrieve_list))
        return calcinfo

    def verify_inputs(self):
        super(VaspCalculation, self).verify_inputs()
        if 'elements' not in self.attrs():
            self._prestore()

    def _prestore(self):
        """Set attributes prior to storing."""
        super(VaspCalculation, self)._prestore()
        self._set_attr('elements', ordered_unique_list(self.inputs.structure.get_ase().get_chemical_symbols()))

    @classmethod
    def _get_potential_linkname(cls, kind):
        """Required for storing multiple input potential nodes."""
        return 'potential_%s' % kind

    @property
    def _parameters(self):
        all_parameters = ChainMap(self.inputs.parameters.get_dict(), self._DEFAULT_PARAMETERS)
        return {k.lower(): v for k, v in all_parameters.items()}

    @property
    def elements(self):
        return self.get_attr('elements')

    def _need_kp(self):
        """
        Return wether an input kpoints node is needed or not.

        :return output:
            True if input kpoints node is needed
            (py:method::VaspCalculation.use_kpoints),
            False otherwise
        needs 'parameters' input to be set
        (py:method::VaspCalculation.use_parameters)
        """
        return not bool('kspacing' in self._parameters and 'kgamma' in self._parameters)

    def _need_chgd(self):
        """
        Test wether an charge_densities input is needed or not.

        :return output:
            True if a chgcar file must be used
            (py:method::NscfCalculation.use_charge_densities),
            False otherwise
        needs 'parameters' input to be set
        (py:method::NscfCalculation.use_parameters)
        """
        ichrg_d = 0 if self._need_wfn() else 2
        icharg = self._parameters.get('icharg', ichrg_d)
        return bool(icharg in [1, 11])

    def _need_wfn(self):
        """
        Test wether a wavefunctions input is needed or not.

        :return output:
            True if a wavecar file must be
            used (py:method::NscfCalculation.use_wavefunctions),
            False otherwise
        needs 'parameters' input to be set
        (py:method::NscfCalculation.use_parameters)
        """
        istrt_d = 1 if self.get_inputs_dict().get('wavefunctions') else 0
        istart = self._parameters.get('istart', istrt_d)
        return bool(istart in [1, 2, 3])

    def _structure(self):
        """Get the input structure as sorted pymatgen structure object."""
        structure = self.inputs.structure
        if not hasattr(structure, 'get_pymatgen'):
            structure = get_data_node('structure', ase=structure.get_ase())
        return structure

    def write_additional(self, tempfolder):
        """Write CHGAR and WAVECAR files if needed."""
        super(VaspCalculation, self).write_additional(tempfolder)
        if self._need_chgd():
            chgcar = tempfolder.get_abs_path('CHGCAR')
            self.write_chgcar(chgcar)
        if self._need_wfn():
            wavecar = tempfolder.get_abs_path('WAVECAR')
            self.write_wavecar(wavecar)

    def write_incar(self, dst):  # pylint: disable=unused-argument
        """
        Converts from parameters node (Dict) to INCAR format and writes to dst.

        Unless otherwise specified, the values specified in _DEFAULT_PARAMETERS are also written to the INCAR file.

        :param dst: absolute path of the file to write to
        """
        incar_dict = ChainMap(self.inputs.parameters.get_dict(), self._DEFAULT_PARAMETERS)
        incar_io = IncarIo(incar_dict=incar_dict)
        incar_io.write(dst)

    def write_poscar(self, dst):  # pylint: disable=unused-argument
        """
        Converts from structures node (StructureData) to POSCAR format and writes to dst.

        :param dst: absolute path of the file to write to
        """
        settings = self.inputs.get('settings')
        settings = settings.get_dict() if settings else {}
        poscar_precision = settings.get('poscar_precision', 10)
        writer = PoscarParser(data=self._structure(), precision=poscar_precision)
        writer.get_quantity('poscar-structure', {})
        writer.write(dst)

    def write_potcar(self, dst):
        """
        Concatenates multiple POTCAR files into one in the same order as the elements appear in POSCAR.

        :param dst: absolute path of the file to write to
        """
        structure = self._structure()
        pot_key = self._get_potential_linkname
        potentials = {symbol: self.inputs.get(pot_key(symbol)) for symbol in structure.get_kind_names()}
        multi_potcar = MultiPotcarIo.from_structure(structure, potentials)
        multi_potcar.write(dst)

    def write_kpoints(self, dst):  # pylint: disable=unused-argument
        """
        Converts from kpoints node (KpointsData) to KPOINTS format and writes to dst.

        :param dst: absolute path of the file to write to
        """
        kpoints = self.inputs.kpoints
        kpoint_parser = KpParser(data=kpoints)
        kpoint_parser.write(dst)

    def write_chgcar(self, dst):  # pylint: disable=unused-argument
        with open(dst, 'w') as out_handle, self.inputs.charge_density.open() as in_handle:
            out_handle.writelines(line for line in in_handle)

    def write_wavecar(self, dst):  # pylint: disable=unused-argument
        with open(dst, 'w') as out_handle, self.inputs.wavefunctions.open() as in_handle:
            out_handle.writelines(line for line in in_handle)

    @classmethod
    def _immigrant_add_inputs(cls, transport, remote_path, sandbox_path, builder, **kwargs):
        from aiida_vasp.calcs.immigrant import get_chgcar_input, get_wavecar_input
        add_wavecar = kwargs.get('use_wavecar') or bool(builder.parameters.get_dict().get('istart', 0))
        add_chgcar = kwargs.get('use_chgcar') or builder.parameters.get_dict().get('icharg', -1) in [1, 11]
        if add_chgcar:
            transport.get(remote_path.join('CHGCAR').strpath, sandbox_path.strpath)
            builder.charge_density = get_chgcar_input(sandbox_path)
        if add_wavecar:
            transport.get(remote_path.join('WAVECAR').strpath, sandbox_path.strpath)
            builder.wavefunctions = get_wavecar_input(sandbox_path)


def ordered_unique_list(in_list):
    """List unique elements in input list, in order of first occurrence"""
    out_list = []
    for i in in_list:
        if i not in out_list:
            out_list.append(i)
    return out_list
