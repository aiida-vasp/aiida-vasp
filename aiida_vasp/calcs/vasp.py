"""
VASP calculation.

-----------------
The calculation class that prepares a specific VASP calculation.
"""
#encoding: utf-8
# pylint: disable=abstract-method
# explanation: pylint wrongly complains about (aiida) Node not implementing query
import os
from aiida.common.exceptions import InputValidationError, ValidationError

from aiida_vasp.parsers.content_parsers.incar import IncarParser
from aiida_vasp.parsers.content_parsers.potcar import MultiPotcarIo
from aiida_vasp.parsers.content_parsers.poscar import PoscarParser
from aiida_vasp.parsers.content_parsers.kpoints import KpointsParser
from aiida_vasp.utils.aiida_utils import get_data_node, get_data_class
from aiida_vasp.calcs.base import VaspCalcBase


class VaspCalculation(VaspCalcBase):
    """
    General-purpose VASP calculation.

    ---------------------------------
    By default retrieves only the 'OUTCAR', 'vasprun.xml', 'EIGENVAL', 'DOSCAR'
    and Wannier90 input / output objects. These objects are deleted after parsing.
    Additional retrieve objects can be specified via the
    ``settings['ADDITIONAL_RETRIEVE_TEMPORARY_LIST']`` input. In addition, if you want to keep
    any objects after parsing, put them in ``settings['ADDITIONAL_RETRIEVE_LIST']`` which is empty
    by default.

    Floating point precision for writing POSCAR objects can be adjusted using
    ``settings['poscar_precision']``, default: 10

    The following assumes you are familiar with the AiiDA data structures and
    how to set up and run an AiiDA calculation in general.

    Example usage::

        from aiida.orm import CalculationFactory, DataFactory
        from aiida.work import submit

        proc = CalculationFactory('vasp.vasp').process()
        inputs = proc.get_inputs_template()
        inputs.parameter = <Dict with INCAR params>
        inputs.structure = <StructureData>
        inputs.kpoints = <KpointsData>
        inputs.settings = <Dict with parser settings etc.>
        inputs.potential = DataFactory('vasp.potcar').get_potcars_from_structure(structure, ...)
        inputs.code = <Code representing vasp on your cluster>

        submit(proc, **inputs)

    Which is very similar to the workchain example.

    Since we do not want the content parsers to know about the AiiDA infrastructure,
    i.e. processes etc. we have no access to the exit codes defined on the CalcJob.
    We thus have to deal with failures in parsing directly in the write calls here.

    """

    _VASP_OUTPUT = 'vasp_output'
    _ALWAYS_RETRIEVE_LIST = ['CONTCAR', 'OUTCAR', 'vasprun.xml', 'EIGENVAL', 'DOSCAR', 'wannier90*', _VASP_OUTPUT]
    _query_type_string = 'vasp.vasp'
    _plugin_type_string = 'vasp.vasp'

    @classmethod
    def define(cls, spec):
        super(VaspCalculation, cls).define(spec)
        # Define the inputs.
        # options is passed automatically.
        spec.input(
            'parameters',
            valid_type=get_data_class('core.dict'),
            help='The VASP input parameters (INCAR).',
        )
        spec.input(
            'dynamics',
            valid_type=get_data_class('core.dict'),
            help='The VASP parameters related to ionic dynamics, e.g. flags to set the selective dynamics',
            required=False,
        )
        spec.input(
            'structure',
            valid_type=(get_data_class('core.structure'), get_data_class('core.cif')),
            help='The input structure (POSCAR).',
        )

        # Need namespace on this as it should also accept keys that are of `kind`. These are unknown
        # until execution.
        spec.input_namespace(
            'potential',
            valid_type=get_data_class('vasp.potcar'),
            help='The potentials (POTCAR).',
            dynamic=True,
        )
        spec.input(
            'kpoints',
            valid_type=get_data_class('core.array.kpoints'),
            help='The kpoints to use (KPOINTS).',
        )
        spec.input(
            'charge_density',
            valid_type=get_data_class('vasp.chargedensity'),
            required=False,
            help='The charge density. (CHGCAR)',
        )
        spec.input(
            'wavefunctions',
            valid_type=get_data_class('vasp.wavefun'),
            required=False,
            help='The wave function coefficients. (WAVECAR)',
        )
        spec.input(
            'settings',
            valid_type=get_data_class('core.dict'),
            required=False,
            help='Additional parameters not related to VASP itself.',
        )
        spec.input(
            'metadata.options.parser_name',
            default='vasp.vasp',
        )

        # Define outputs.
        # remote_folder and retrieved are passed automatically
        spec.output(
            'misc',
            valid_type=get_data_class('core.dict'),
            help='The output parameters containing smaller quantities that do not depend on system size.',
        )
        spec.output(
            'structure',
            valid_type=get_data_class('core.structure'),
            required=False,
            help='The output structure.',
        )
        spec.output(
            'kpoints',
            valid_type=get_data_class('core.array.kpoints'),
            required=False,
            help='The output k-points.',
        )
        spec.output(
            'trajectory',
            valid_type=get_data_class('core.array.trajectory'),
            required=False,
            help='The output trajectory data.',
        )
        spec.output(
            'chgcar',
            valid_type=get_data_class('vasp.chargedensity'),
            required=False,
            help='The output charge density CHGCAR file.',
        )
        spec.output(
            'wavecar',
            valid_type=get_data_class('vasp.wavefun'),
            required=False,
            help='The output plane wave coefficients file.',
        )
        spec.output(
            'bands',
            valid_type=get_data_class('core.array.bands'),
            required=False,
            help='The output band structure.',
        )
        spec.output(
            'forces',
            valid_type=get_data_class('core.array'),
            required=False,
            help='The output forces.',
        )
        spec.output(
            'stress',
            valid_type=get_data_class('core.array'),
            required=False,
            help='The output stress.',
        )
        spec.output(
            'dos',
            valid_type=get_data_class('core.array'),
            required=False,
            help='The output dos.',
        )
        spec.output(
            'energies',
            valid_type=get_data_class('core.array'),
            required=False,
            help='The output total energies.',
        )
        spec.output(
            'projectors',
            valid_type=get_data_class('core.array'),
            required=False,
            help='The output projectors of decomposition.',
        )
        spec.output(
            'dielectrics',
            valid_type=get_data_class('core.array'),
            required=False,
            help='The output dielectric functions.',
        )
        spec.output(
            'born_charges',
            valid_type=get_data_class('core.array'),
            required=False,
            help='The output Born effective charges.',
        )
        spec.output(
            'hessian',
            valid_type=get_data_class('core.array'),
            required=False,
            help='The output Hessian matrix.',
        )
        spec.output(
            'dynmat',
            valid_type=get_data_class('core.array'),
            required=False,
            help='The output dynamical matrix.',
        )
        spec.output(
            'charge_density',
            valid_type=get_data_class('core.array'),
            required=False,
            help='The output charge density.',
        )
        spec.output(
            'magnetization_density',
            valid_type=get_data_class('core.array'),
            required=False,
            help='The output magnetization density.',
        )
        spec.output(
            'site_magnetization',
            valid_type=get_data_class('core.dict'),
            required=False,
            help='The output of the site magnetization',
        )
        spec.output_namespace('custom_outputs', required=False, dynamic=True)
        spec.exit_code(0, 'NO_ERROR', message='the sun is shining')
        spec.exit_code(
            350,
            'ERROR_NO_RETRIEVED_FOLDER',
            message='the retrieved folder data node could not be accessed.',
            invalidates_cache=True,
        )
        spec.exit_code(
            351,
            'ERROR_NO_RETRIEVED_TEMPORARY_FOLDER',
            message='the retrieved_temporary folder data node could not be accessed.',
            invalidates_cache=True,
        )
        spec.exit_code(
            352,
            'ERROR_CRITICAL_MISSING_OBJECT',
            message='an object that is marked by the parser as critical is missing.',
            invalidates_cache=True,
        )
        spec.exit_code(
            333,
            'ERROR_VASP_DID_NOT_EXECUTE',
            message='VASP did not produce any output and did likely not execute properly.',
            invalidates_cache=True,
        )

        # 700 series of the errors catches VASP execution related problems
        spec.exit_code(
            700,
            'ERROR_DID_NOT_FINISH',
            message='Calculation did not reach the end of execution.',
            invalidates_cache=True,
        )
        spec.exit_code(
            701,
            'ERROR_ELECTRONIC_NOT_CONVERGED',
            message='The electronic structure is not converged.',
        )
        spec.exit_code(
            702,
            'ERROR_IONIC_NOT_CONVERGED',
            message='The ionic relaxation is not converged.',
        )
        spec.exit_code(
            703,
            'ERROR_VASP_CRITICAL_ERROR',
            message='VASP calculation encountered a critical error: {error_message}.',
        )
        spec.exit_code(
            704,
            'ERROR_DIAGNOSIS_OUTPUTS_MISSING',
            message=
            'Outputs for diagnosis are missing, please make sure `run_status` and `notifications` quantities are requested for parsing.',
        )

        spec.exit_code(
            1001,
            'ERROR_PARSING_OBJECT_FAILED',
            message='parsing an object has failed.',
        )
        spec.exit_code(
            1002,
            'ERROR_NOT_ABLE_TO_PARSE_QUANTITY',
            message='the parser is not able to parse the {quantity} quantity',
        )
        spec.exit_code(
            1003,
            'ERROR_RECOVERY_PARSING_OF_XML_FAILED',
            message=
            'the vasprun.xml was truncated and recovery parsing failed to parse at least one of the requested quantities: {quantities}, '
            'very likely the VASP calculation did not run properly',
        )
        spec.exit_code(
            1004,
            'ERROR_NOT_ABLE_TO_CREATE_NODE',
            message='the parser is not able to compose one or more output nodes: {nodes}',
        )

    def prepare_for_submission(self, folder):
        """
        Add all objects to the list of objects to be retrieved.

        Notice that we here utilize both the retrieve batch of objects, which are always stored after retrieval and
        the temporary retrieve list which is automatically cleared after parsing.
        """
        calcinfo = super().prepare_for_submission(folder)

        # Combine stdout and stderr into vasp_output so that the stream parser can parse it later.
        calcinfo.codes_info[0].stdout_name = self._VASP_OUTPUT
        calcinfo.codes_info[0].join_files = True

        # Still need the exceptions in case settings is not defined on inputs
        # Check if we want to store all always retrieve objects
        try:
            store = self.inputs.settings.base.attributes.get('ALWAYS_STORE', default=True)
        except AttributeError:
            store = True
        try:
            additional_retrieve_list = self.inputs.settings.base.attributes.get('ADDITIONAL_RETRIEVE_LIST', default=[])
        except AttributeError:
            additional_retrieve_list = []
        try:
            additional_retrieve_temp_list =\
                self.inputs.settings.base.attributes.get('ADDITIONAL_RETRIEVE_TEMPORARY_LIST', default=[])  # pylint: disable=invalid-name
        except AttributeError:
            additional_retrieve_temp_list = []
        if store:
            calcinfo.retrieve_list = list(set(self._ALWAYS_RETRIEVE_LIST + additional_retrieve_list))
            calcinfo.retrieve_temporary_list = additional_retrieve_temp_list  # pylint: disable=invalid-name
        else:
            calcinfo.retrieve_temporary_list = list(set(self._ALWAYS_RETRIEVE_LIST + additional_retrieve_temp_list))  # pylint: disable=invalid-name
            calcinfo.retrieve_list = additional_retrieve_list
        try:
            provenance_exclude_list = self.inputs.settings.base.attributes.get('PROVENANCE_EXCLUDE_LIST', default=[])
        except AttributeError:
            provenance_exclude_list = []
        # Always include POTCAR in the exclude list (not added to the repository, regardless of store)
        calcinfo.provenance_exclude_list = list(set(provenance_exclude_list + ['POTCAR']))

        return calcinfo

    def verify_inputs(self):
        super().verify_inputs()
        _parameters = self.inputs.parameters.get_dict()
        _lorbit = _parameters.get('lorbit', 0)
        if 'settings' in self.inputs:
            _settings = self.inputs.settings.get_dict()
        else:
            _settings = {}
        _site_magnetization = _settings.get('parser_settings', {}).get('add_site_magnetization', False)
        if _site_magnetization and _lorbit < 10:
            raise InputValidationError(f'Site magnetization requires "LORBIT>=10", value given {_lorbit}')

    @property
    def _parameters(self):
        """Make sure all parameters are lowercase."""
        all_parameters = self.inputs.parameters.get_dict()
        try:
            return {k.lower(): v for k, v in all_parameters.items()}
        except KeyError:
            return {}

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
        return not bool('kspacing' in self._parameters or 'kgamma' in self._parameters)

    def _need_chgcar(self):
        """
        Test wether an charge_densities input is needed or not.

        :return output:
            True if CHGCAR must be present
            (py:method::NscfCalculation.use_charge_densities),
            False otherwise
        needs 'parameters' input to be set
        (py:method::NscfCalculation.use_parameters)
        """
        ichrg_d = 0 if self._need_wavecar() else 2
        icharg = self._parameters.get('icharg', ichrg_d)
        return bool(icharg in [1, 11])

    def _need_wavecar(self):
        """
        Test wether a wavefunctions input is needed or not.

        :return output:
            True if WAVECAR must be present
            used (py:method::NscfCalculation.use_wavefunctions),
            False otherwise
        needs 'parameters' input to be set
        (py:method::NscfCalculation.use_parameters)
        """
        istrt_d = 1 if self.inputs.get('wavefunctions') else 0
        istart = self._parameters.get('istart', istrt_d)
        return bool(istart in [1, 2, 3])

    def _structure(self):
        """
        Get the input structure as AiiDa StructureData.

        This is required in order to support CifData as input as well.
        """
        structure = self.inputs.structure
        if not hasattr(structure, 'get_pymatgen'):
            structure = get_data_node('core.structure', ase=structure.get_ase())
        return structure

    def write_additional(self, folder, calcinfo):
        """Write CHGAR and WAVECAR if needed."""
        super().write_additional(folder, calcinfo)
        # a list of object names to be copied
        remote_copy_fnames = [os.path.split(entry[1])[1] for entry in calcinfo.remote_copy_list]
        if self._need_chgcar():
            # If we restart, we do not require inputs, but we should have a basic check
            # that CHGCAR is present
            if not self._is_restart():
                self.write_chgcar('CHGCAR', calcinfo)
            else:
                remote_folder = self.inputs.restart_folder
                if 'CHGCAR' not in remote_copy_fnames:
                    raise FileNotFoundError('Could not find CHGCAR in {}'.format(remote_folder.get_remote_path()))
        if self._need_wavecar():
            # If we restart, we do not require inputs, but we should have a basic check
            # that WAVECAR is present
            if not self._is_restart():
                self.write_wavecar('WAVECAR', calcinfo)
            else:
                remote_folder = self.inputs.restart_folder
                if 'WAVECAR' not in remote_copy_fnames:
                    raise FileNotFoundError('Could not find WAVECAR in {}'.format(remote_folder.get_remote_path()))

    def write_incar(self, dst, validate_tags=True):  # pylint: disable=unused-argument
        """
        Write the INCAR.

        Passes the parameters node (Dict) from to the INCAR parser for
        preparation and writes to dst.

        :param dst: absolute path of the object to write to
        """
        try:
            incar_parser = IncarParser(data=self.inputs.parameters, validate_tags=validate_tags)
            incar_parser.write(dst)
        except SystemExit as parser_error:
            raise ValidationError('The INCAR content did not pass validation.') from parser_error

    def write_poscar(self, dst):  # pylint: disable=unused-argument
        """
        Write the POSCAR.

        Passes the structures node (StructureData) to the POSCAR parser for
        preparation and writes to dst.

        :param dst: absolute path of the object to write to
        """
        settings = self.inputs.get('settings')
        settings = settings.get_dict() if settings else {}
        poscar_precision = settings.get('poscar_precision', 10)
        dynamics = self.inputs.get('dynamics')
        options = None
        if dynamics is not None:
            dynamics = dynamics.get_dict()
            positions_dof = dynamics.get('positions_dof')
            if positions_dof is not None:
                options = {'positions_dof': positions_dof}

        try:
            poscar_parser = PoscarParser(data=self._structure(), precision=poscar_precision, options=options)
            poscar_parser.write(dst)
        except SystemExit as parser_error:
            raise ValidationError('The POSCAR content did not pass validation.') from parser_error

    def write_potcar(self, dst):
        """
        Concatenates multiple POTCARs into one in the same order as the elements appear in POSCAR.

        :param dst: absolute path of the object to write to
        """
        structure = self._structure()
        multi_potcar = MultiPotcarIo.from_structure(structure, self.inputs.potential)
        multi_potcar.write(dst)

    def write_kpoints(self, dst):  # pylint: disable=unused-argument
        """
        Write the KPOINTS.

        Passes the kpoints node (KpointsData) to the KPOINTS parser for
        preparation and writes to dst.

        :param dst: absolute path of the object to write to
        """
        try:
            kpoint_parser = KpointsParser(data=self.inputs.kpoints)
            kpoint_parser.write(dst)
        except SystemExit as parser_error:
            raise ValidationError('The KPOINTS content did not pass validation.') from parser_error

    def write_chgcar(self, dst, calcinfo):  # pylint: disable=unused-argument
        charge_density = self.inputs.charge_density
        calcinfo.local_copy_list.append((charge_density.uuid, charge_density.filename, dst))

    def write_wavecar(self, dst, calcinfo):  # pylint: disable=unused-argument
        wave_functions = self.inputs.wavefunctions
        calcinfo.local_copy_list.append((wave_functions.uuid, wave_functions.filename, dst))

    @classmethod
    def immigrant(cls, code, remote_path, **kwargs):
        """
        Returns VaspImmigrant class and associated inputs. This method will be obsolete at v3.0

        :param code: a Code instance for the code originally used.
        :param remote_path: The directory on the code's computer in which the simulation was run.
        :param metadata: dict.
        :param settings: dict. Used for non-default parsing instructions, etc.
        :param potential_family: str.
        :param potential_mapping: dict.
        :param use_wavecar: bool. Try to read WAVECAR.
        :param use_chgcar bool. Try to read CHGCAR.
        """

        from aiida_vasp.calcs.immigrant import VaspImmigrant  # pylint: disable=import-outside-toplevel

        proc_cls = VaspImmigrant
        builder = proc_cls.get_builder_from_folder(code, str(remote_path), **kwargs)
        options = {'max_wallclock_seconds': 1, 'resources': {'num_machines': 1, 'num_mpiprocs_per_machine': 1}}
        builder.metadata = kwargs.get('metadata', {'options': options})
        options = builder.metadata.get('options', options)
        max_wallclock_seconds = options.get('max_wallclock_seconds', 1)
        resources = options.get('resources', {'num_machines': 1, 'num_mpiprocs_per_machine': 1})
        builder.metadata['options']['max_wallclock_seconds'] = max_wallclock_seconds  # pylint: disable=no-member
        builder.metadata['options']['resources'] = resources  # pylint: disable=no-member
        settings = kwargs.get('settings', {})
        settings.update({'import_from_path': str(remote_path)})
        builder.settings = get_data_node('core.dict', dict=settings)

        return proc_cls, builder


def ordered_unique_list(in_list):
    """List unique elements in input list, in order of first occurrence."""
    out_list = []
    for i in in_list:
        if i not in out_list:
            out_list.append(i)
    return out_list


def ordered_unique_symbols(structure):
    """
    Return a list of ordered unique symbols in the structure
    """
    symbols = [structure.get_kind(kindname).symbol for kindname in structure.get_site_kindnames()]
    return ordered_unique_list(symbols)
