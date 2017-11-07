#encoding: utf-8
# pylint: disable=abstract-method
# explanation: pylint wrongly complains about (aiida) Node not implementing query
"""VASP - Calculation: Generic run using pymatgen for file preparation"""
try:
    from collections import ChainMap
except ImportError:
    from chainmap import ChainMap

from aiida.orm import DataFactory

from aiida_vasp.calcs.base import VaspCalcBase, Input
from aiida_vasp.io.pymatgen.vasprun import get_data_node

PARAMETER_CLS = DataFactory('parameter')
SINGLEFILE_CLS = DataFactory('singlefile')


class VaspCalculation(VaspCalcBase):
    """
    General-purpose VASP calculation.

    By default retrieves only the 'OUTCAR', 'vasprun.xml', 'EIGENVAL', 'DOSCAR' and Wannier90 input / output files,
    but additional retrieve files can be specified via the 'settings['ADDITIONAL_RETRIEVE_LIST']' input.
    """

    default_parser = 'vasp.vasp'
    parameters = Input(types='parameter', doc='VASP INCAR parameters.')
    structure = Input(types=['structure', 'cif'])
    paw = Input(types='vasp.paw', param='kind')
    kpoints = Input(types='array.kpoints')
    settings = Input(
        types='parameter', doc='Additional settings for the calculation.')
    charge_density = Input(
        types='vasp.chargedensity',
        doc=
        'chargedensity node: should be obtained from the output of a selfconsistent SCF calculation (written to CHGCAR)'
    )
    wavefunctions = Input(
        types='vasp.wavefun',
        doc='wavefunction node: to speed up convergence for continuation jobs')

    _DEFAULT_PARAMETERS = {}
    _ALWAYS_RETRIEVE_LIST = [
        'OUTCAR', 'vasprun.xml', 'EIGENVAL', 'DOSCAR', ('wannier90*', '.', 0)
    ]

    def _prepare_for_submission(self, tempfolder, inputdict):
        """add EIGENVAL, DOSCAR, and all files starting with wannier90 to
        the list of files to be retrieved."""
        calcinfo = super(VaspCalculation, self)._prepare_for_submission(
            tempfolder, inputdict)
        try:
            additional_retrieve_list = inputdict['settings'].get_attr(
                'ADDITIONAL_RETRIEVE_LIST')
        except (KeyError, AttributeError):
            additional_retrieve_list = []
        calcinfo.retrieve_list = list(
            set(self._ALWAYS_RETRIEVE_LIST + additional_retrieve_list))
        return calcinfo

    def verify_inputs(self, inputdict, *args, **kwargs):
        super(VaspCalculation, self).verify_inputs(inputdict, *args, **kwargs)
        self.check_input(inputdict, 'parameters')
        self.check_input(inputdict, 'structure')
        if 'elements' not in self.attrs():
            self._prestore()
        for kind in self.elements:
            self.check_input(inputdict, self._get_paw_linkname(kind))
        self.check_input(inputdict, 'kpoints', self._need_kp)
        self.check_input(inputdict, 'charge_density', self._need_chgd)
        self.check_input(inputdict, 'wavefunctions', self._need_wfn)

    def _prestore(self):
        """
        set attributes prior to storing
        """
        super(VaspCalculation, self)._prestore()
        self._set_attr(
            'elements',
            ordered_unique_list(
                self.inp.structure.get_ase().get_chemical_symbols()))

    @classmethod
    def _get_paw_linkname(cls, kind):
        """required for storing multiple input paw nodes"""
        return 'paw_%s' % kind

    @property
    def _parameters(self):
        all_parameters = ChainMap(self.inp.parameters.get_dict(),
                                  self._DEFAULT_PARAMETERS)
        return {k.lower(): v for k, v in all_parameters.items()}

    @property
    def elements(self):
        return self.get_attr('elements')

    def _need_kp(self):
        """
        return wether an input kpoints node is needed or not.
        :return output:
            True if input kpoints node is needed
            (py:method::VaspCalculation.use_kpoints),
            False otherwise
        needs 'parameters' input to be set
        (py:method::VaspCalculation.use_parameters)
        """
        return not bool('kspacing' in self._parameters
                        and 'kgamma' in self._parameters)

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

    def write_additional(self, tempfolder, inputdict):
        """write CHGAR and WAVECAR files if needed"""
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
        from aiida_vasp.io.incar import dict_to_incar
        with open(dst, 'w') as incar:
            incar.write(
                dict_to_incar(
                    ChainMap(self.inp.parameters.get_dict(),
                             self._DEFAULT_PARAMETERS)))

    def write_poscar(self, inputdict, dst):  # pylint: disable=unused-argument
        """
        converts from structures node (StructureData) to POSCAR format
        and writes to dst

        :param inputdict: required by baseclass
        :param dst: absolute path of the file to write to
        """
        from pymatgen.io.vasp.inputs import Poscar
        structure = self.inp.structure
        if not hasattr(structure, 'get_pymatgen'):
            structure = get_data_node('structure', ase=structure.get_ase())
        pmg_structure = structure.get_pymatgen()
        pmg_structure.sort()
        writer = Poscar(pmg_structure)
        writer.write_file(dst)

    def write_potcar(self, inputdict, dst):
        """
        Concatenates multiple paw files into a POTCAR

        :param inputdict: required by baseclass
        :param dst: absolute path of the file to write to
        """
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
        with open(dst, 'w') as potcar_f:
            sp.check_call(catcom, stdout=potcar_f)

    def write_kpoints(self, inputdict, dst):  # pylint: disable=unused-argument
        """
        converts from kpoints node (KpointsData) to KPOINTS format
        and writes to dst

        :param inputdict: required by baseclass
        :param dst: absolute path of the file to write to
        """
        kpoints = self.inp.kpoints
        if kpoints.get_attrs().get('mesh'):
            self._write_kpoints_mesh(dst)
        elif kpoints.get_attrs().get('array|kpoints'):
            self._write_kpoints_list(dst)
        else:
            raise AttributeError('you supplied an empty kpoints node')

    def _write_kpoints_mesh(self, dst):
        """Write kpoints in mesh format to the destination file `dst`"""
        kpoints = self.inp.kpoints
        mesh, offset = kpoints.get_kpoints_mesh()
        kpmtemp = ("Automatic mesh\n"
                   "0\n"
                   "Gamma\n"
                   "{N[0]} {N[1]} {N[2]}\n"
                   "{s[0]} {s[1]} {s[2]}\n")
        with open(dst, 'w') as kpoints:
            kps = kpmtemp.format(N=mesh, s=offset)
            kpoints.write(kps)

    def _write_kpoints_list(self, dst):
        """Write a list of kpoints to the destination file `dst`"""
        kpoints = self.inp.kpoints
        if 'array|weights' in kpoints.get_attrs():
            kpl, weights = kpoints.get_kpoints(also_weights=True)
        else:
            kpl = kpoints.get_kpoints()
            weights = [1.] * kpl.shape[0]
        kpoint_weights = list(zip(kpl, weights))

        kpls = '\n'.join([
            '{k[0]} {k[1]} {k[2]} {w}'.format(k=k, w=w)
            for k, w in kpoint_weights
        ])
        kps = ("Explicit list\n"
               "{N}\n"
               "Direct\n"
               "{klist}\n").format(
                   N=len(kpoint_weights), klist=kpls)
        with open(dst, 'w') as kpoints:
            kpoints.write(kps)

    def write_chgcar(self, inputdict, dst):  # pylint: disable=unused-argument
        import shutil
        shutil.copyfile(self.inp.charge_density.get_file_abs_path(), dst)

    def write_wavecar(self, inputdict, dst):  # pylint: disable=unused-argument
        import shutil
        shutil.copyfile(self.inp.wavefunctions.get_file_abs_path(), dst)


def ordered_unique_list(in_list):
    """List unique elements in input list, in order of first occurrence"""
    out_list = []
    for i in in_list:
        if i not in out_list:
            out_list.append(i)
    return out_list
