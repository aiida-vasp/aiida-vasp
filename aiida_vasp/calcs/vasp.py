#encoding: utf-8
# pylint: disable=abstract-method
# explanation: pylint wrongly complains about (aiida) Node not implementing query
"""VASP - Calculation: Generic run using pymatgen_aiida for file preparation"""
try:
    from collections import ChainMap
except ImportError:
    from chainmap import ChainMap

from aiida.orm import DataFactory

from aiida_vasp.io.incar import IncarIo
from aiida_vasp.io.potcar import MultiPotcarIo
from aiida_vasp.io.poscar import PoscarParser
from aiida_vasp.io.kpoints import KpParser
from aiida_vasp.utils.aiida_utils import get_data_node
from aiida_vasp.calcs.base import VaspCalcBase, Input

PARAMETER_CLS = DataFactory('parameter')
SINGLEFILE_CLS = DataFactory('singlefile')


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

        inputs.parameter = <ParameterData with INCAR params>
        inputs.structure = <StructureData> or <CifData>
        inputs.kpoints = <KpointsData>
        inputs.settings = <ParameterData with parser settings etc>
        inputs.potential = DataFactory('vasp.potcar').get_potcars_from_structure(structure, ...)
        inputs.code = <Code representing vasp on your cluster>
        inputs._options = <Computer, resources, etc, AiiDA specific stuff>

        submit(proc, **inputs)

    Example Low-Level usage::

        ## assuming already set up incar (ParameterData), structure, kpoints, settings, etc
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

    default_parser = 'vasp.vasp'
    parameters = Input(types='parameter', doc='VASP INCAR parameters.')
    structure = Input(types=['structure', 'cif'])
    potential = Input(types='vasp.potcar', param='kind')
    kpoints = Input(types='array.kpoints')
    settings = Input(types='parameter', doc='Additional settings for the calculation.')
    charge_density = Input(
        types='vasp.chargedensity',
        doc='chargedensity node: should be obtained from the output of a selfconsistent SCF calculation (written to CHGCAR)')
    wavefunctions = Input(types='vasp.wavefun', doc='wavefunction node: to speed up convergence for continuation jobs')

    _DEFAULT_PARAMETERS = {}
    _ALWAYS_RETRIEVE_LIST = ['CONTCAR', 'OUTCAR', 'vasprun.xml', 'EIGENVAL', 'DOSCAR', ('wannier90*', '.', 0)]
    _query_type_string = 'vasp.vasp'
    _plugin_type_string = 'vasp.vasp'

    def _prepare_for_submission(self, tempfolder, inputdict):
        """Add EIGENVAL, DOSCAR, and all files starting with wannier90 to the list of files to be retrieved."""
        calcinfo = super(VaspCalculation, self)._prepare_for_submission(tempfolder, inputdict)
        try:
            additional_retrieve_list = inputdict['settings'].get_attr('ADDITIONAL_RETRIEVE_LIST')
        except (KeyError, AttributeError):
            additional_retrieve_list = []
        calcinfo.retrieve_list = list(set(self._ALWAYS_RETRIEVE_LIST + additional_retrieve_list))
        return calcinfo

    def verify_inputs(self, inputdict, *args, **kwargs):
        super(VaspCalculation, self).verify_inputs(inputdict, *args, **kwargs)
        self.check_input(inputdict, 'parameters')
        self.check_input(inputdict, 'structure')
        if 'elements' not in self.attrs():
            self._prestore()
        for kind in self._structure().get_kind_names():
            self.check_input(inputdict, self._get_potential_linkname(kind))
        self.check_input(inputdict, 'kpoints', self._need_kp)
        self.check_input(inputdict, 'charge_density', self._need_chgd)
        self.check_input(inputdict, 'wavefunctions', self._need_wfn)

    def _prestore(self):
        """Set attributes prior to storing."""
        super(VaspCalculation, self)._prestore()
        self._set_attr('elements', ordered_unique_list(self.inp.structure.get_ase().get_chemical_symbols()))

    @classmethod
    def _get_potential_linkname(cls, kind):
        """Required for storing multiple input potential nodes."""
        return 'potential_%s' % kind

    @property
    def _parameters(self):
        all_parameters = ChainMap(self.inp.parameters.get_dict(), self._DEFAULT_PARAMETERS)
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
        structure = self.inp.structure
        if not hasattr(structure, 'get_pymatgen'):
            structure = get_data_node('structure', ase=structure.get_ase())
        return structure

    def write_additional(self, tempfolder, inputdict):
        """Write CHGAR and WAVECAR files if needed."""
        super(VaspCalculation, self).write_additional(tempfolder, inputdict)
        if self._need_chgd():
            chgcar = tempfolder.get_abs_path('CHGCAR')
            self.write_chgcar(inputdict, chgcar)
        if self._need_wfn():
            wavecar = tempfolder.get_abs_path('WAVECAR')
            self.write_wavecar(inputdict, wavecar)

    def write_incar(self, inputdict, dst):  # pylint: disable=unused-argument
        """
        Converts from parameters node (ParameterData) to INCAR format and writes to dst.

        Unless otherwise specified, the values specified in _DEFAULT_PARAMETERS are also written to the INCAR file.

        :param inputdict: required by baseclass
        :param dst: absolute path of the file to write to
        """
        incar_dict = ChainMap(self.inp.parameters.get_dict(), self._DEFAULT_PARAMETERS)
        incar_io = IncarIo(incar_dict=incar_dict)
        incar_io.write(dst)

    def write_poscar(self, inputdict, dst):  # pylint: disable=unused-argument
        """
        Converts from structures node (StructureData) to POSCAR format and writes to dst.

        :param inputdict: required by baseclass
        :param dst: absolute path of the file to write to
        """
        settings = inputdict.get('settings')
        settings = settings.get_dict() if settings else {}
        poscar_precision = settings.get('poscar_precision', 10)
        writer = PoscarParser(data=self._structure(), precision=poscar_precision)
        writer.get_quantity('poscar-structure', {})
        writer.write(dst)

    def write_potcar(self, inputdict, dst):
        """
        Concatenates multiple POTCAR files into one in the same order as the elements appear in POSCAR.

        :param inputdict: required by baseclass
        :param dst: absolute path of the file to write to
        """
        structure = self._structure()
        pot_key = self._get_potential_linkname
        potentials = {symbol: inputdict[pot_key(symbol)] for symbol in structure.get_kind_names()}
        multi_potcar = MultiPotcarIo.from_structure(structure, potentials)
        multi_potcar.write(dst)

    def write_kpoints(self, inputdict, dst):  # pylint: disable=unused-argument
        """
        Converts from kpoints node (KpointsData) to KPOINTS format and writes to dst.

        :param inputdict: required by baseclass
        :param dst: absolute path of the file to write to
        """
        kpoints = self.inp.kpoints

        kpoint_parser = KpParser(data=kpoints)
        kpoint_parser.write(dst)

    def write_chgcar(self, inputdict, dst):  # pylint: disable=unused-argument
        import shutil
        shutil.copyfile(self.inp.charge_density.get_file_abs_path(), dst)

    def write_wavecar(self, inputdict, dst):  # pylint: disable=unused-argument
        import shutil
        shutil.copyfile(self.inp.wavefunctions.get_file_abs_path(), dst)

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
