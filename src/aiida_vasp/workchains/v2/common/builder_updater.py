"""

This module provides a set of updater classes and utility functions for constructing and managing
AiiDA process builders for VASP-based workflows. The updaters encapsulate logic for applying
presets, setting calculation options, managing input sets, and updating workflow-specific
settings for various VASP workchains, including standard calculations, relaxations, NEB,
convergence tests, and band structure calculations.

Key Classes:

- VaspPresetConfig: Handles loading and managing preset configurations from YAML files.
- BaseBuilderUpdater: Base class for builder updaters, providing common methods for builder manipulation.
- VaspBuilderUpdater: Updater for standard VASP calculations.
- VaspNEBUpdater: Updater for NEB (nudged elastic band) calculations.
- VaspRelaxUpdater: Updater for relaxation workflows.
- VaspMultiStageRelaxUpdater: Updater for multi-stage relaxation workflows.
- VaspConvUpdater: Updater for convergence testing workflows.
- VaspBandUpdater: Updater for band structure workflows.
- VaspHybridBandUpdater: Updater for hybrid functional band structure workflows.

Key Utilities:

- update_dict_node: Utility to safely update AiiDA Dict nodes.
- builder_to_dict: Converts a builder to a Python dictionary for inspection.
- incar_dict_to_relax_settings: Extracts relaxation settings from INCAR parameters.
- is_specified: Checks if any values are set in a ProcessBuilderNamespace.

The module is designed to facilitate programmatic and reproducible setup of VASP workflows
in AiiDA, supporting both interactive and automated use cases.

"""

from __future__ import annotations

import logging
from copy import deepcopy
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Union
from warnings import warn

from aiida import orm
from aiida.common.exceptions import InputValidationError
from aiida.common.extendeddicts import AttributeDict
from aiida.engine import run_get_node, submit
from aiida.engine.processes.builder import ProcessBuilder, ProcessBuilderNamespace
from aiida.plugins import WorkflowFactory
from ase.visualize import view
from yaml import safe_load

from aiida_vasp.utils.opthold import ConvOptions
from aiida_vasp.workchains.v2.bands import BandOptions
from aiida_vasp.workchains.v2.relax import VaspMultiStageRelaxWorkChain

from ..inputset.base import convert_lowercase
from ..inputset.pmgset import PymatgenInputSet
from ..inputset.vaspsets import VASPInputSet
from ..relax import RelaxOptions
from .transform import neb_interpolate

DEFAULT_PRESET = 'VaspPreset'
DEFAULT_INPUTSET = 'UCLRelaxSet'


__all__ = (
    'VaspBandUpdater',
    'VaspBuilderUpdater',
    'VaspConvUpdater',
    'VaspHybridBandUpdater',
    'VaspNEBUpdater',
    'VaspPresetConfig',
    'VaspRelaxUpdater',
)


def get_library_path() -> Path:
    """
    Get the path where the YAML files are stored within this package.

    :returns: Path to the library directory containing YAML configuration files
    :rtype: pathlib.Path
    """
    return Path(__file__).parent


# Template for setting options
OPTIONS_TEMPLATES = {
    'SGE': {
        'resources': {'tot_num_mpiprocs': 1, 'parallel_env': 'mpi'},
        'max_wallclock_seconds': 3600,
        'import_sys_environment': False,
    },
    'FW': {
        'resources': {
            'tot_num_mpiprocs': 1,
        },
        'max_wallclock_seconds': 3600,
    },
    'SLURM': {
        'resources': {
            'num_machines': 1,
        },
        'max_wallclock_seconds': 3600,
        'import_sys_environment': False,
    },
    'ARCHER2': {
        'resources': {
            'tot_num_mpiprocs': 128,
            'num_machines': 1,
        },
        'max_wallclock_seconds': 3600,
        'import_sys_environment': False,
        'mpirun_extra_params': ['--distribution=block:block', '--hint=nomultithread'],
        'account': 'e05-power-dos',
        'queue_name': 'standard',
        'qos': 'standard',
    },
}


@dataclass
class VaspPresetConfig:
    """Class to store the preset for VaspBuilderUpdater"""

    name: str
    inputset: str
    default_code: str
    code_specific: dict = field(default_factory=dict)
    default_options: dict = field(default_factory=dict)
    default_settings: dict = field(default_factory=dict)
    default_inputset_overrides: dict = field(default_factory=dict)
    default_relax_settings: dict = field(default_factory=dict)
    default_band_settings: dict = field(default_factory=dict)

    @classmethod
    def from_file(cls, fname: str) -> VaspPresetConfig:
        """
        Load preset configuration from a YAML file.

        Searches for the configuration file in the package library path and user's
        home directory (~/.aiida-vasp).

        :param fname: Name of the configuration file (without .yaml extension)
        :type fname: str

        :returns: VaspPresetConfig instance loaded from file
        :rtype: VaspPresetConfig

        :raises RuntimeError: If the preset definition file cannot be found
        """
        _load_paths = (get_library_path(), Path('~/.aiida-vasp').expanduser())
        for parent in _load_paths:
            target_path = parent / (fname + '.yaml')
            if target_path.is_file():
                break
        if target_path is None:
            raise RuntimeError(f'Cannot find preset definition for {fname}')

        with open(target_path, encoding='utf-8', mode='r') as fhandle:
            data = safe_load(fhandle)
        return cls(**data)

    def get_code_specific_options(self, code: str, namespace: str) -> dict[str, Any]:
        """
        Return code-specific options for a given namespace.

        If code-specific options exist, they are merged with the default options
        for the namespace, with code-specific options taking precedence.

        :param code: Name/identifier of the computational code
        :type code: str
        :param namespace: Configuration namespace (e.g., 'options', 'settings')
        :type namespace: str

        :returns: Dictionary containing the merged options
        :rtype: dict
        """
        if code in self.code_specific:
            if namespace in self.code_specific[code]:
                code_specific = self.code_specific[code][namespace]
                default = getattr(self, f'default_{namespace}', {})
                if default is None:
                    default = {}
                default = deepcopy(default)
                default.update(code_specific)
                return default
        return deepcopy(getattr(self, f'default_{namespace}'))


class BaseBuilderUpdater:
    """Base class for builder updater"""

    def __init__(
        self,
        preset_name: str | None = None,
        builder: ProcessBuilder | None = None,
        verbose: bool = False,
        inputset_name: str | None = None,
        set_name: str | None = None,
    ) -> None:
        """Instantiate a pipeline"""
        # Configure the builder

        assert hasattr(self, 'WF_ENTRYPOINT'), 'WF_ENTRYPOINT must be specified by the class'
        self.verbose = verbose
        if builder is None:
            builder = WorkflowFactory(self.WF_ENTRYPOINT).get_builder()
        self._builder = builder
        if preset_name is None:
            preset_name = DEFAULT_PRESET
        self.preset_name = preset_name
        self.preset = VaspPresetConfig.from_file(preset_name)
        if set_name is not None:
            inputset_name = set_name
            warn("The 'set_name' parameter is deprecated, use 'inputset_name' instead.")
        self.inputset_name = inputset_name if inputset_name is not None else self.preset.inputset

    @property
    def builder(self) -> ProcessBuilder:
        """
        The builder to be used for launching the calculation.

        :returns: Process builder instance
        :rtype: ProcessBuilder
        """
        return self._builder

    def submit(self) -> orm.WorkChainNode:
        """
        Submit the workflow to the daemon and return the workchain node.

        :returns: The submitted workchain node
        :rtype: orm.WorkChainNode
        """
        return submit(self.builder)

    def run_get_node(self, verbose: bool = True) -> orm.WorkChainNode:
        """
        Run the workflow with the current python process.

        :param verbose: If True, print debugging information for failed calculations
        :type verbose: bool

        :returns: Tuple containing the workflow outputs and the workchain node
        :rtype: orm.WorkChainNode
        """
        output = run_get_node(self.builder)
        # Verbose output (for debugging)
        if not output.node.is_finished_ok and verbose:
            for node in output.node.called_descendants:
                if isinstance(node, orm.CalcJobNode):
                    stdout = node.called[0].outputs.retrieved.get_object_content('vasp_output')
                    print(node, 'STDOUT:', stdout)
                    print(node, 'Retrieved files:', node.retrieved.list_object_names())
                    script = node.base.repository.get_object_content('_aiidasubmit.sh')
                    print(node, 'Submission script:', script)
                    print(node, 'Exit_message', node.exit_message)
        return output

    def _get_help(self, namespace: str, print_to_stdout: bool = True, inout: str = 'inputs') -> str | None:
        """
        Return the help message for a given namespace.

        The `.` syntax for the namespace is supported for nested namespaces.

        :param namespace: Namespace path (e.g., 'vasp.parameters')
        :type namespace: str
        :param print_to_stdout: Whether to print help to stdout or return it
        :type print_to_stdout: bool
        :param inout: Whether to get help for 'inputs' or 'outputs'
        :type inout: str

        :returns: Help message if print_to_stdout is False, otherwise None
        :rtype: str or None
        """
        levels = namespace.split('.')
        data_dict = self.builder._process_spec.get_description()[inout]
        for key in levels:
            data_dict = data_dict[key]

        if print_to_stdout is True:
            print(data_dict.get('help', 'No help message information found'))
        else:
            return data_dict.get('help', 'No help message information found')

    def get_output_help(self, namespace: str, print_to_stdout: bool = True) -> str | None:
        """
        Return the help message for a given output namespace.

        :param namespace: Output namespace path
        :type namespace: str
        :param print_to_stdout: Whether to print help to stdout or return it
        :type print_to_stdout: bool

        :returns: Help message if print_to_stdout is False, otherwise None
        :rtype: str or None
        """
        self._get_help(namespace, print_to_stdout=print_to_stdout, inout='outputs')

    def get_input_help(self, namespace: str, print_to_stdout: bool = True) -> str | None:
        """
        Return the help message for a given input namespace.

        :param namespace: Input namespace path
        :type namespace: str
        :param print_to_stdout: Whether to print help to stdout or return it
        :type print_to_stdout: bool

        :returns: Help message if print_to_stdout is False, otherwise None
        :rtype: str or None
        """
        self._get_help(namespace, print_to_stdout=print_to_stdout, inout='inputs')


class VaspBuilderUpdater(BaseBuilderUpdater):
    WF_ENTRYPOINT = 'vasp.v2.vasp'
    DEFAULT_INPUTSET = DEFAULT_INPUTSET

    def __init__(
        self,
        preset_name: str | None = None,
        builder: ProcessBuilder | None = None,
        root_namespace: ProcessBuilderNamespace | None = None,
        code: str | None = None,
        verbose: bool = False,
        inputset_name: str | None = None,
    ) -> None:
        """
        Initialise the update object.

        :param builder: The ``ProcessBuilder`` or ``ProcessBuilderNamespace`` to be used for setting
          standard VaspWorkChain inputs.

        :param root_namespace: The namespace to be assumed to be the *root*, e.g. where the input structure
          should be specified. The v2 series of workchain in aiida-vasp usually has the StructureData input
          port at the top level interface, although there are a few exceptions.
        :param preset_name: The name of the Preset to be used for the updater.
        :param code: The code to be used for the calculation. If not specified, the default code from the
          preset will be used.
        :param verbose: If True, print additional information during the update.
        :param set_name: The name of the input set to be used. If not specified, the default input set from the preset
          will be used.

        returns: An instance of VaspBuilderUpdater with the specified preset and builder.
        """
        super().__init__(preset_name=preset_name, builder=builder, verbose=verbose, inputset_name=inputset_name)
        # Define the root namespace - e.g. the VaspWorkChain namespace where structure should be specified
        if root_namespace is None:
            self.root_namespace = self._builder
        else:
            self.root_namespace = root_namespace

        self.namespace_vasp = self._builder
        self.code = self.preset.default_code if code is None else code

    @property
    def reference_structure(self) -> orm.StructureData:
        """
        Reference structure used for setting kpoints and other calculations.

        :returns: The structure data node used as reference
        :rtype: orm.StructureData
        """
        return self.root_namespace.structure

    def clear(self) -> None:
        """
        Clear all nodes set in the VASP and root namespaces.

        Resets parameters, options, settings, kpoints, potential family/mapping,
        structure, and metadata label to None.
        """
        self.namespace_vasp.parameters = None
        self.namespace_vasp.options = None
        self.namespace_vasp.settings = None
        self.namespace_vasp.kpoints = None
        self.namespace_vasp.potential_family = None
        self.namespace_vasp.potential_mapping = None

        self.root_namespace.structure = None
        self.root_namespace.metadata.label = None

    def apply_preset(
        self,
        initial_structure: orm.StructureData,
        code: str | None = None,
        label: str | None = None,
        overrides: dict[str, Any] | None = None,
        inputset_name: str | None = None,
    ) -> VaspBuilderUpdater:
        """
        Apply the complete preset configuration to the builder.

        This method applies the input set, sets the computational code, options,
        settings, and label according to the preset configuration.

        :param initial_structure: Structure to be used for the calculation
        :type initial_structure: orm.StructureData
        :param code: Computational code to use (defaults to preset default)
        :type code: str or None
        :param label: Label for the calculation (defaults to structure label)
        :type label: str or None
        :param overrides: Dictionary of parameter overrides
        :type overrides: dict or None
        :param inputset_name: Name of input set to use (defaults to preset default)
        :type inputset_name: str or None

        :returns: Self for method chaining
        :rtype: VaspBuilderUpdater
        """
        if code is None:
            code = self.code
            logging.info(f'Using code {code}')
        self.use_inputset(
            initial_structure,
            set_name=self.inputset_name if inputset_name is None else inputset_name,
            overrides=overrides,
            apply_preset=True,
            code=code,
        )
        self.set_code(code=code)
        self.set_options(code=code, apply_preset=True)
        self.set_settings(code=code, apply_preset=True)
        self.set_label(label)
        return self

    def use_inputset(
        self,
        structure: orm.StructureData,
        set_name: str | None = None,
        overrides: dict[str, Any] | None = None,
        apply_preset: bool = False,
        code: str | None = None,
        structure_port_name: str = 'structure',
        pmg_kwargs: dict[str, Any] | None = None,
    ) -> VaspBuilderUpdater:
        """
        Update the inputs ports for the VASP calculation.

        :param structure: The structure to be used for the calculation.
        :param set_name: The name of the input set to be used.
        :param overrides: Any overrides to be applied to the input set.
        :param apply_preset: Whether to apply the preset options.
        :param code: The code to be used for the calculation.
        :param structure_node_name: The name of in put port where the structure should be set.
        :param pmg_kwargs: Additional kwargs to be passed to pymatgen's InputSet when using a pymatgen
            inputset.

        :returns : self, the VaspBuilderUpdater instance with the input set applied.
        """
        # Use the default inputset name if not defined
        if set_name is None:
            set_name = self.DEFAULT_INPUTSET
        if overrides is None:
            overrides = {}

        if apply_preset:
            if code is None:
                code = self.preset.default_code
            overrides_ = convert_lowercase(self.preset.get_code_specific_options(code, 'inputset_overrides'))
            overrides_.update(overrides)
        else:
            overrides_ = overrides

        if set_name in PymatgenInputSet.KNOWN_SETS:
            inset = PymatgenInputSet(set_name, overrides=overrides_, verbose=self.verbose, pmg_kwargs=pmg_kwargs)
            # PymatgenInputSet uses explicit kpoints
            kpt = inset.get_kpoints(structure)
            # Use kpoints mesh or kspacing provided by pymatgen inputset
            if kpt is not None:
                self.namespace_vasp.kpoints = kpt
            else:
                self.namespace_vasp.kpoints_spacing = orm.Float(inset.get_kpoints_spacing(structure))
        else:
            inset = VASPInputSet(set_name, overrides=overrides_, verbose=self.verbose)
            self.namespace_vasp.kpoints_spacing = orm.Float(inset.get_kpoints_spacing())

        self.namespace_vasp.parameters = orm.Dict(dict={'incar': inset.get_input_dict(structure, raw_python=True)})
        try:
            self.namespace_vasp.potential_family = orm.Str(inset.get_potcar_family())
        except InputValidationError:
            warn(
                f'Error validating potential family {inset.get_potcar_family()} for input set {set_name}. '
                'Potential family will not be set. '
            )
        self.namespace_vasp.potential_mapping = orm.Dict(dict=inset.get_pp_mapping(structure))
        setattr(self.root_namespace, structure_port_name, structure)
        return self

    def set_kspacing(self, kspacing: float) -> VaspBuilderUpdater:
        """
        Set the k-point spacing and remove any existing k-point mesh.

        :param kspacing: K-point spacing value in inverse Angstroms
        :type kspacing: float

        :returns: Self for method chaining
        :rtype: VaspBuilderUpdater
        """
        self.namespace_vasp.kpoints_spacing = orm.Float(kspacing)
        if self.namespace_vasp.kpoints:
            del self.namespace_vasp.kpoints
        return self

    def set_potential_family(self, family: str) -> VaspBuilderUpdater:
        """
        Set the potential family for the VASP calculation.

        :param family: Name of the potential family
        :type family: str

        :returns: Self for method chaining
        :rtype: VaspBuilderUpdater
        """
        self.namespace_vasp.potential_family = orm.Str(family)
        return self

    def set_potential_mapping(self, mapping: dict[str, str]) -> VaspBuilderUpdater:
        """
        Set the potential mapping for the VASP calculation.

        :param mapping: Dictionary mapping element symbols to potential names
        :type mapping: dict[str, str]

        :returns: Self for method chaining
        :rtype: VaspBuilderUpdater
        """
        self.namespace_vasp.potential_mapping = orm.Dict(dict=mapping)
        return self

    update_kspacing = set_kspacing

    @property
    def parameters(self) -> Union[orm.Dict, None]:
        """
        Return the parameters node containing INCAR settings.

        :returns: Parameters node or None if not set
        :rtype: orm.Dict or None
        """
        return self.namespace_vasp.parameters

    @property
    def settings(self) -> Union[orm.Dict, None]:
        """
        Return the settings node for VASP calculation options.

        :returns: Settings node or None if not set
        :rtype: orm.Dict or None
        """
        return self.namespace_vasp.settings

    @property
    def options(self) -> Union[orm.Dict, None]:
        """
        Return the computational options node.

        :returns: Options node or None if not set
        :rtype: orm.Dict or None
        """
        return self.namespace_vasp.options

    def set_code(self, code: str | orm.Code | None = None) -> VaspBuilderUpdater:
        """
        Set the Code for the VASP calculation.

        :param code: Code identifier string or Code node (defaults to preset default)
        :type code: str, orm.Code, or None

        :returns: Self for method chaining
        :rtype: VaspBuilderUpdater
        """
        if code is None:
            code = self.preset.default_code
        if isinstance(code, str):
            code = orm.load_code(code)

        self.namespace_vasp.code = code
        return self

    def update_code(self, code: str | orm.Code) -> VaspBuilderUpdater:
        warn('update_code is deprecated, use set_code instead', DeprecationWarning)
        return self.set_code(code)

    def set_incar(self, *args: Any, **kwargs: Any) -> VaspBuilderUpdater:
        """
        Update INCAR parameters for the VASP calculation.

        :param args: Positional arguments passed to dict constructor
        :param kwargs: INCAR parameter key-value pairs

        :returns: Self for method chaining
        :rtype: VaspBuilderUpdater
        """
        if self.namespace_vasp.parameters is None:
            self.namespace_vasp.parameters = orm.Dict(dict={'incar': {}})
        content = dict(*args, **kwargs)
        node = update_dict_node(self.namespace_vasp.parameters, content, 'incar')
        self.namespace_vasp.parameters = node
        return self

    def update_incar(self, *args: Any, **kwargs: Any) -> VaspBuilderUpdater:
        warn('update_incar is deprecated, use set_incar instead', DeprecationWarning)
        return self.set_incar(*args, **kwargs)

    def set_options(
        self, *args: Any, code: str | None = None, apply_preset: bool = False, **kwargs: Any
    ) -> VaspBuilderUpdater:
        """
        Set computational options for the VASP calculation.

        :param args: Positional arguments passed to dict constructor
        :param code: Code name for code-specific options
        :type code: str or None
        :param apply_preset: Whether to apply preset-defined options
        :type apply_preset: bool
        :param kwargs: Option key-value pairs

        :returns: Self for method chaining
        :rtype: VaspBuilderUpdater
        """
        if apply_preset:
            if code is None:
                code = self.preset.default_code
            odict = self.preset.get_code_specific_options(code, 'options')
            odict.update(dict(*args, **kwargs))
        else:
            odict = dict(*args, **kwargs)

        if self.namespace_vasp.options is None:
            self.namespace_vasp.options = orm.Dict(odict)
        else:
            self.namespace_vasp.options = update_dict_node(self.namespace_vasp.options, odict)
        return self

    def update_options(self, *args: Any, **kwargs: Any) -> VaspBuilderUpdater:
        warn('update_options is deprecated, use set_options instead', DeprecationWarning)
        return self.set_options(*args, **kwargs)

    def set_kpoints_mesh(self, mesh: list[int], offset: list[float] = (0.0, 0.0, 0.0)) -> VaspBuilderUpdater:
        """
        Set explicit k-points mesh for the calculation.
        The plugin generates the KPOINTS file with a Gamma-centered mesh.
        Monkhorst-Pack meshes can be applied by using the offset parameter, e.g. (0.5, 0.5, 0.5)

        :param mesh: K-point mesh dimensions [nx, ny, nz]
        :type mesh: List[int]
        :param offset: K-point mesh offset [ox, oy, oz]
        :type offset: List[float]

        :returns: Self for method chaining
        :rtype: VaspBuilderUpdater
        """
        kpoints = orm.KpointsData()
        kpoints.set_cell_from_structure(self.reference_structure)
        kpoints.set_kpoints_mesh(mesh, list(offset))
        self.namespace_vasp.kpoints = kpoints
        try:
            del self.namespace_vasp.kpoints_spacing
        except KeyError:
            pass
        return self

    def update_kpoints_mesh(self, mesh: list[int], offset: list[float]) -> VaspBuilderUpdater:
        warn('update_kpoints_mesh is deprecated, use set_kpoints_mesh instead', DeprecationWarning)
        return self.set_kpoints_mesh(mesh, offset)

    def set_settings(
        self, *args: Any, code: str | None = None, apply_preset: bool = False, **kwargs: Any
    ) -> VaspBuilderUpdater:
        """
        Set the 'settings' input port.

        :param args: Positional arguments passed to dict constructor
        :param code: Code name for code-specific settings
        :type code: str or None
        :param apply_preset: Whether to apply preset-defined settings
        :type apply_preset: bool
        :param kwargs: Setting key-value pairs

        :returns: Self for method chaining
        :rtype: VaspBuilderUpdater
        """
        if apply_preset:
            if code is None:
                code = self.preset.default_code
            sdict = self.preset.get_code_specific_options(code, 'settings')
            # Apply use supplied contents
            sdict.update(dict(*args, **kwargs))
        else:
            sdict = dict(*args, **kwargs)

        if self.namespace_vasp.settings is None:
            self.namespace_vasp.settings = orm.Dict(sdict)
        else:
            self.namespace_vasp.settings = update_dict_node(self.namespace_vasp.settings, sdict)
        return self

    def update_settings(self, *args: Any, **kwargs: Any) -> VaspBuilderUpdater:
        warn('update_settings is deprecated, use set_settings instead', DeprecationWarning)
        return self.set_settings(*args, **kwargs)

    def set_label(self, label: str | None = None) -> VaspBuilderUpdater:
        """
        Set the top-level label for the calculation.

        :param label: Label string (defaults to structure label if available)
        :type label: str or None

        :returns: Self for method chaining
        :rtype: VaspBuilderUpdater
        """
        if label is None:
            # Default to the label of the structure if available
            if 'structure' in self.root_namespace:
                label = self.root_namespace.structure.label
        self.root_namespace.metadata.label = label
        return self

    def update_label(self, label: str | None = None) -> VaspBuilderUpdater:
        warn('update_label is deprecated, use set_label instead', DeprecationWarning)
        return self.set_label(label)

    def set_resources(self, *args: Any, **kwargs: Any) -> VaspBuilderUpdater:
        """
        Update computational resources in the options.
        NOTE: The available options can be found in the documentation of the Calculation class. These are
        identical to those used in the metadata.options namespace.

        :param args: Positional arguments passed to dict constructor
        :param kwargs: Resource key-value pairs

        :returns: Self for method chaining
        :rtype: VaspBuilderUpdater

        :raises RuntimeError: If options are not set before calling this method
        """
        if self.namespace_vasp.options is None:
            raise RuntimeError('Please set the options before setting resources')
        resources = dict(self.namespace_vasp.options['resources'])
        resources.update(*args, **kwargs)
        self.namespace_vasp.options = update_dict_node(self.namespace_vasp.options, resources, 'resources')
        return self

    def update_resources(self, *args: Any, **kwargs: Any) -> VaspBuilderUpdater:
        warn('update_resources is deprecated, use set_resources instead', DeprecationWarning)
        return self.set_resources(*args, **kwargs)

    def _set_options(
        self,
        option_class: type,
        option_name: str,
        target_namespace: ProcessBuilder | ProcessBuilderNamespace,
        **kwargs: Any,
    ) -> VaspBuilderUpdater:
        """
        Set options using a specific option class.

        :param option_class: Class used to validate and structure options
        :param option_name: Name of the option attribute in the target namespace
        :type option_name: str
        :param target_namespace: Namespace where options should be set
        :type target_namespace: ProcessBuilder or ProcessBuilderNamespace
        :param kwargs: Option key-value pairs

        :returns: Self for method chaining
        """
        if getattr(target_namespace, option_name) is None:
            current_option = option_class()
        else:
            current_option = option_class(**getattr(target_namespace, option_name).get_dict())
        for key, value in kwargs.items():
            setattr(current_option, key, value)
        setattr(target_namespace, option_name, current_option.aiida_dict())
        return self


class VaspNEBUpdater(VaspBuilderUpdater):
    WF_ENTRYPOINT = 'vasp.neb'

    @property
    def reference_structure(self):
        """
        Return the reference structure for NEB calculations.

        :returns: Initial structure for NEB calculation
        :rtype: orm.StructureData
        """
        return self.namespace_vasp.initial_structure

    def apply_preset(
        self,
        structure_init: orm.StructureData,
        structure_final: orm.StructureData,
        code: str | None = None,
        label: str | None = None,
        interpolate: bool = True,
        nimages: int = 5,
        **kwargs: Any,
    ) -> VaspNEBUpdater:
        super().apply_preset(structure_init, code, label, **kwargs)
        self.set_final_structure(structure_final)
        if interpolate:
            self.set_interpolated_images(nimages)
        else:
            logging.info('Not interpolating images, please set with .set_neb_image(images)')
        self.update_incar(images=nimages)

        return self

    def use_inputset(
        self,
        initial_structure: orm.StructureData,
        set_name: str | None = None,
        overrides: dict[str, Any] | None = None,
        apply_preset: bool = False,
        code: str | None = None,
    ) -> VaspNEBUpdater:
        super().use_inputset(
            structure=initial_structure,
            set_name=set_name,
            overrides=overrides,
            apply_preset=apply_preset,
            code=code,
            structure_port_name='initial_structure',
        )
        return self

    def set_label(self, label: str | None = None) -> VaspNEBUpdater:
        """
        Set the toplevel label, default to the label of the structure"""
        if label is None:
            label = self.root_namespace.initial_structure.label
        self.root_namespace.metadata.label = label
        return self

    def set_final_structure(self, final_structure: orm.StructureData) -> VaspNEBUpdater:
        """
        Set the final structure for NEB calculation.

        :param final_structure: Final structure for the NEB path
        :type final_structure: orm.StructureData

        :returns: Self for method chaining
        :rtype: VaspNEBUpdater
        """
        self.namespace_vasp.final_structure = final_structure
        return self

    def set_neb_images(self, images: list | dict | AttributeDict) -> VaspNEBUpdater:
        """
        Set the intermediate NEB images.

        :param images: List of structures or dictionary mapping image names to structures
        :type images: list, dict, or AttributeDict

        :returns: Self for method chaining
        :rtype: VaspNEBUpdater
        """
        if isinstance(images, list):
            output = {f'image_{i:02d}': image for i, image in enumerate(images)}
        elif isinstance(images, (dict, AttributeDict)):
            output = images
        self.namespace_vasp.neb_images = output
        return self

    def set_interpolated_images(self, nimages: int) -> VaspNEBUpdater:
        """
        Generate and set interpolated images between initial and final structures.

        This requires the initial and final structures to be set already.
        Also updates the final image with PBC issues fixed.

        :param nimages: Number of intermediate images to generate
        :type nimages: int

        :returns: Self for method chaining
        :rtype: VaspNEBUpdater
        """

        initial = self.namespace_vasp.initial_structure
        final = self.namespace_vasp.final_structure
        assert initial
        assert final
        # Generate interpolated images and fix PBC issues if applicable
        interpolated = neb_interpolate(initial, final, orm.Int(nimages))
        images = {key: value for key, value in interpolated.items() if not ('init' in key or 'final' in key)}
        self.namespace_vasp.neb_images = images
        # Update the final image - make sure that is atoms are not wrapped around
        self.set_final_structure(interpolated['image_final'])
        return self

    def view_images(self, *args: Any, **kwargs: Any) -> None:
        """
        Visualize the NEB images using ASE viewer.

        Displays all images including initial, intermediate, and final structures.

        Hint: In a notebook environment, you can pass "viewer='weas'" to use weas-widget viewer.
        This requires the ase-weas-widget package to be installed.
        """
        view(
            map(
                lambda x: x.get_ase(),
                [self.builder.initial_structure, *self.builder.neb_images.values(), self.builder.final_structure],
            ),
            *args,
            **kwargs,
        )


class VaspRelaxUpdater(VaspBuilderUpdater):
    """
    An updater for VaspRelaxWorkChain
    """

    WF_ENTRYPOINT = 'vasp.v2.relax'

    def __init__(
        self,
        preset_name: str | None = None,
        builder: ProcessBuilder | None = None,
        override_vasp_namespace: ProcessBuilderNamespace | None = None,
        namespace_relax: ProcessBuilderNamespace | None = None,
        code: str | None = None,
    ) -> None:
        super().__init__(preset_name=preset_name, builder=builder, code=code, root_namespace=builder)
        # The primary VASP namespace is under builder.vasp
        if override_vasp_namespace is None:
            self.namespace_vasp = self._builder.vasp
        else:
            self.namespace_vasp = override_vasp_namespace

        if namespace_relax is None:
            self.namespace_relax = self._builder
        else:
            self.namespace_relax = namespace_relax

    def use_inputset(self, *args: Any, set_name: str | None = None, **kwargs: Any) -> VaspRelaxUpdater:
        super().use_inputset(*args, set_name=set_name, **kwargs)
        if set_name is not None:
            if set_name in PymatgenInputSet.KNOWN_SETS:
                incar_in = self.namespace_vasp.parameters['incar']
                incar_out, relax_update = incar_dict_to_relax_settings(incar_in)
                self.set_incar(**incar_out)
                self.set_relax_settings(**relax_update)

    def apply_preset(
        self,
        structure: orm.StructureData,
        code: str | None = None,
        label: str | None = None,
        **kwargs: Any,
    ) -> VaspRelaxUpdater:
        out = super().apply_preset(structure, code, label, **kwargs)
        self.set_relax_settings()
        return out

    def set_relax_settings(self, **kwargs: Any) -> VaspRelaxUpdater:
        """
        Set/update RelaxOptions controlling the operation of the workchain.

        :param kwargs: Relaxation option key-value pairs

        :returns: Self for method chaining
        :rtype: VaspRelaxUpdater
        """
        self._set_options(RelaxOptions, 'relax_settings', self.namespace_relax, **kwargs)
        return self

    update_relax_settings = set_relax_settings

    def clear_relax_settings(self) -> VaspRelaxUpdater:
        """
        Reset any existing relax options to defaults.

        :returns: Self for method chaining
        :rtype: VaspRelaxUpdater
        """
        self.namespace_relax.relax_settings = RelaxOptions().aiida_dict()
        return self

    def clear(self) -> VaspRelaxUpdater:
        """
        Clear all settings including relax-specific settings.

        :returns: Self for method chaining
        :rtype: VaspRelaxUpdater
        """
        super().clear()
        self.clear_relax_settings()
        return self


class VaspMultiStageRelaxUpdater(VaspRelaxUpdater):
    """
    An updater for VaspRelaxWorkChain
    """

    WF_ENTRYPOINT = 'vasp.v2.staged_relax'

    def __init__(
        self,
        preset_name: str | None = None,
        builder: ProcessBuilder | None = None,
        override_vasp_namespace: ProcessBuilderNamespace | None = None,
        namespace_relax: ProcessBuilderNamespace | None = None,
        code: str | None = None,
    ) -> None:
        if builder is None:
            builder = VaspMultiStageRelaxWorkChain.get_builder()
        if override_vasp_namespace is None:
            override_vasp_namespace = builder.relax.vasp
        if namespace_relax is None:
            namespace_relax = builder.relax
        super().__init__(
            preset_name=preset_name,
            builder=builder,
            code=code,
            override_vasp_namespace=override_vasp_namespace,
            namespace_relax=namespace_relax,
        )


class VaspConvUpdater(VaspBuilderUpdater):
    """Update for VaspConvergenceWorkChain"""

    WF_ENTRYPOINT = 'vasp.v2.converge'

    def apply_preset(
        self, initial_structure: orm.StructureData, code: str | None = None, label: str | None = None, **kwargs: Any
    ) -> VaspBuilderUpdater:
        super().apply_preset(initial_structure, code, label, **kwargs)
        self.set_conv_settings()
        return self

    def set_conv_settings(self, **kwargs: Any) -> VaspConvUpdater:
        """
        Set the convergence testing settings.

        :param kwargs: Convergence option key-value pairs

        :returns: Self for method chaining
        :rtype: VaspConvUpdater
        """
        self._set_options(ConvOptions, 'conv_settings', self.builder, **kwargs)
        return self


class VaspBandUpdater(VaspBuilderUpdater):
    """Updater for VaspBandsWorkChain"""

    WF_ENTRYPOINT = 'vasp.v2.bands'

    def __init__(
        self,
        preset_name: str | None = None,
        builder: ProcessBuilder | None = None,
        override_vasp_namespace: ProcessBuilderNamespace | None = None,
        code: str | None = None,
    ) -> None:
        super().__init__(preset_name=preset_name, builder=builder, code=code, root_namespace=builder)
        # The primary VASP namespace is under builder.vasp
        if override_vasp_namespace is None:
            self.namespace_vasp = self.builder.scf
        else:
            self.namespace_vasp = override_vasp_namespace

    def get_relax_updater(self) -> VaspRelaxUpdater:
        """
        Return the relax updater for this band structure calculation.

        The relax updater can be used to populate the `.relax` namespace which will
        trigger the relaxation of the structure before band structure calculation.

        :returns: VaspRelaxUpdater instance configured for this band calculation
        :rtype: VaspRelaxUpdater
        """
        # Apply relax settings if requested
        relax = VaspRelaxUpdater(
            preset_name=self.preset_name,
            builder=self.builder,
            namespace_relax=self.builder.relax,
            override_vasp_namespace=self.builder.relax.vasp,
            code=self.code,
        )
        return relax

    def apply_preset(
        self, structure: orm.StructureData, run_relax: bool = False, label: str | None = None, **kwargs: Any
    ) -> VaspBandUpdater:
        super().apply_preset(structure, label=label, **kwargs)

        # Specify the relaxation and NAC namespace
        if run_relax:
            relax_upd = self.get_relax_updater()
            relax_upd.apply_preset(structure, label=label, **kwargs)
        self.set_band_settings()
        return self

    def set_band_settings(self, **kwargs: Any) -> VaspBandUpdater:
        """
        Set band structure calculation specific settings.

        :param kwargs: Band calculation option key-value pairs

        :returns: Self for method chaining
        :rtype: VaspBandUpdater
        """
        self._set_options(BandOptions, 'band_settings', self.root_namespace, **kwargs)
        return self


class VaspHybridBandUpdater(VaspBandUpdater):
    """Updater for VaspHybridBandsWorkChain"""

    WF_ENTRYPOINT = 'vasp.v2.hybrid_bands'


# class VaspAutoPhononUpdater(VaspBuilderUpdater):
#     """Updater for VaspAutoPhononWorkChain"""

#     WF_ENTRYPOINT = 'vasp.v2.phonopy'

#     def __init__(self, builder: ProcessBuilder):
#         """Initialise with an existing ProcessBuilder for VaspAutoPhononWorkChain"""
#         super().__init__(builder.singlepoint, root_namespace=builder)

#     def set_phonon_settings(self, options):
#         """
#         Update the phonon-related options

#         example::

#           {
#             'primitive_matrix': 'auto',
#             'supercell_matrix': [2, 2, 2],    # Supercell matrix
#             'mesh': 30,                       # Mesh for DOS/PDOS/thermal properties
#           }


#         """
#         self.root_namespace.phonon_settings = orm.Dict(options)
#         return self

#     def update_from_config(self, structure: orm.StructureData, config: dict):
#         """
#         Update the builder from a configuration dictionary.

#         The dictionary must has a ``singlepoint`` key holding the configurations for singlepoint
#         calculations, and a ``phonon_options`` for Phonopy options to be used.
#         The ``relax`` and ``nac`` keys are optional.
#         """

#         super().update_from_config(structure, config['singlepoint'])

#         # Specify the relaxation and NAC namespace
#         if 'relax' in config:
#             relax_upd = VaspRelaxUpdater(
#                 self.root_namespace,
#                 namespace_relax=self.root_namespace.relax,
#                 override_vasp_namespace=self.root_namespace.relax.vasp,
#             )
#             relax_upd.update_from_config(structure, config['relax'])

#         if 'nac' in config:
#             nac_upd = VaspBuilderUpdater(self.root_namespace.nac, root_namespace=self.root_namespace)
#             nac_upd.update_from_config(structure, config['nac'])

#         # Update the phonon settings
#         self.set_phonon_settings(config['phonon_settings'])
#         return self

#     def set_kpoints_mesh(self, mesh, offset) -> None:
#         """Use mesh for kpoints"""
#         kpoints = orm.KpointsData()
#         # Use the reference supercell structure
#         kpoints.set_cell_from_structure(self.reference_structure)
#         kpoints.set_kpoints_mesh(mesh, offset)
#         self.namespace_vasp.kpoints = kpoints
#         if self.namespace_vasp.kpoints_spacing:
#             del self.namespace_vasp.kpoints_spacing
#         return self

#     def _get_singlepoint_supercell(self) -> orm.StructureData:
#         """Obtain the supercell for the singlepoint calculation"""
#         import numpy as np
#         from ase.build import make_supercell

#         ref = self.root_namespace.structure.get_ase()

#         # The sueprcell matrix should be a vector or a matrix
#         mat = np.array(self.root_namespace.phonon_settings['supercell_matrix'])
#         if mat.size == 3:
#             mat = np.diag(mat)

#         # Convention of phonopy - the supercell matrix is the transpose of that would be used
#         # for ase
#         return orm.StructureData(ase=make_supercell(ref, mat.T))

#     def show_builder(self):
#         """Print stuff defined in the builder"""
#         pprint(builder_to_dict(self.root_namespace, unpack=True))


def is_specified(port_namespace: ProcessBuilderNamespace) -> bool:
    """
    Check if there is anything specified under a PortNamespace.

    :param port_namespace: Namespace to check for specified values
    :type port_namespace: ProcessBuilderNamespace

    :returns: True if any values are specified in the namespace
    :rtype: bool
    """
    return any(map(bool, port_namespace.values()))


def update_dict_node(
    node: orm.Dict,
    content: dict[str, Any],
    namespace: str | None = None,
    reuse_if_possible: bool = True,
) -> orm.Dict:
    """
    Update a Dict node with new content.

    Optionally updates a specific namespace within the Dict node.
    If the node is stored and immutable, creates a new node with updated content.

    :param node: The Dict node to update
    :type node: orm.Dict
    :param content: Dictionary content to merge into the node
    :type content: dict
    :param namespace: Optional namespace key within the Dict to update
    :type namespace: str or None
    :param reuse_if_possible: Whether to reuse the existing node if content is unchanged
    :type reuse_if_possible: bool

    :returns: Updated Dict node (may be the same or a new node)
    :rtype: orm.Dict
    """
    # Get pure-python dictionary
    dtmp = node.get_dict()
    dtmp_backup = None
    if reuse_if_possible and node.is_stored:
        dtmp_backup = deepcopy(dtmp)
    if namespace:
        dtmp.get(namespace, {}).update(content)
    else:
        dtmp.update(content)
    if node.is_stored:
        # There is no need to update the node if the content is the same as before
        if reuse_if_possible and dtmp == dtmp_backup:
            return node
        # The content is different, but the node is immutable, so we create a new node
        return orm.Dict(dict=dtmp)
    node.update_dict(dtmp)
    return node


def builder_to_dict(builder: ProcessBuilder, unpack: bool = True) -> dict[str, Any]:
    """
    Convert a builder to a dictionary and optionally unpack certain nodes.

    When unpacked, the resulting dictionary cannot be used for `submit`/`run`.
    The primary usage of the resulting dictionary is for pretty printing.

    :param builder: Process builder to convert
    :type builder: ProcessBuilder
    :param unpack: Whether to unpack Dict/List nodes to Python objects
    :type unpack: bool

    :returns: Dictionary representation of the builder
    :rtype: dict
    """
    data = {}
    for key, value in builder._data.items():
        if hasattr(value, '_data'):
            value_ = builder_to_dict(builder[key])
        if unpack:
            if isinstance(value, orm.Dict):
                value_ = value.get_dict()
            elif isinstance(value, orm.List):
                value_ = value.get_list()
            else:
                value_ = value
        data[key] = value_
    return data


def incar_dict_to_relax_settings(incar_in: dict[str, Any]) -> tuple[dict[str, Any], dict[str, Any]]:
    """
    Convert INCAR tags to relax_settings and remove them from INCAR.

    Extracts relaxation-specific INCAR parameters (NSW, IBRION, EDIFFG) and
    converts them to equivalent relax_settings options.

    :param incar_in: Input dictionary containing INCAR parameters
    :type incar_in: dict

    :returns: Tuple of (updated_incar_dict, relax_settings_dict)
    :rtype: tuple
    """
    # Convert INCAR tags to relax_settings
    updated = {}
    incar_out = dict(incar_in)
    nsw = incar_out['incar'].pop('nsw', None)
    if nsw is not None:
        updated['steps'] = nsw
    # Convert ibrion
    ibrion = incar_out['incar'].pop('ibrion', None)
    if ibrion == 1:
        updated['algo'] = 'rd'
    if ibrion == 2:
        updated['algo'] = 'cg'
    # Convert ediffg
    ediffg = incar_out['incar'].pop('ediffg', None)
    if ediffg is not None:
        updated['force_cutoff'] = ediffg
    return incar_out, updated
