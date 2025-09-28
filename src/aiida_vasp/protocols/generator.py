"""
Input generators based on protocols

This module aimed at interactive post-generation update for the builder created
by `.get_builder_from_protocol` method of various workchain classes.
"""

from __future__ import annotations

from copy import deepcopy
from dataclasses import dataclass, field
from itertools import chain
from pathlib import Path
from typing import Any

from aiida import orm
from aiida.engine import run_get_node, submit
from aiida.engine.processes.builder import ProcessBuilderNamespace
from aiida.plugins import WorkflowFactory
from yaml import safe_load

from aiida_vasp.utils.dict_merge import recursive_merge

__all__ = [
    'VaspBandsInputGenerator',
    'VaspConvergenceInputGenerator',
    'VaspHybridBandsInputGenerator',
    'VaspInputGenerator',
    'VaspRelaxInputGenerator',
]


DEFAULT_PRESET = 'default_preset'
DEFAULT_PROTOCOL = 'balanced'


def get_library_path() -> Path:
    """
    Get the path where the YAML files are stored within this package.

    :returns: Path to the library directory containing YAML configuration files
    :rtype: pathlib.Path
    """
    return Path(__file__).parent / 'presets'


def list_protocol_presets() -> list[Path]:
    """
    List all available presets in the package.
    """
    _load_paths = (get_library_path(), Path('~/.aiida-vasp/protocol_presets').expanduser())
    presets = []
    for parent in _load_paths:
        files = chain(parent.glob('*.yaml'), parent.glob('*.yml'))
        for file in files:
            presets.append(file.absolute())
    return presets


@dataclass
class PresetConfig:
    """Class to store the preset for inputs"""

    name: str
    default_protocol: str
    default_code: str
    code_specific: dict = field(default_factory=dict)
    default_options: dict = field(default_factory=dict)
    default_settings: dict = field(default_factory=dict)
    protocol_overrides: dict = field(default_factory=dict)
    default_relax_settings: dict = field(default_factory=dict)
    default_band_settings: dict = field(default_factory=dict)

    @classmethod
    def from_file(cls, fname: str) -> PresetConfig:
        """
        Load preset configuration from a YAML file.

        Searches for the configuration file in the package library path and user's
        home directory (~/.aiida-vasp/protocol_presets/).

        :param fname: Name of the configuration file (without .yaml extension)
        :type fname: str

        :returns: ProtocolPresetConfig instance loaded from file
        :rtype: ProtocolPresetConfig

        :raises RuntimeError: If the preset definition file cannot be found
        """
        _load_paths = (get_library_path(), Path('~/.aiida-vasp/protocol_presets').expanduser())
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
        return deepcopy(
            getattr(
                self,
                f'default_{namespace}',
                {},
            )
        )


class BaseInputGenerator:
    """
    BaseClass for all protocol builder updaters

    The protocol updater serves two purposes:
    - Generating a builder based on a user-defined "preset", e.g. with options and overrides pre-loaded
    - Allow interactive modifications of common parameters such as incar tag's, resources and options.
    """

    WF_ENTRYPOINT = 'vasp.vasp'

    def __init__(
        self,
        preset_name: str = 'default',
        protocol: str | None = None,
        verbose: bool = False,
    ) -> None:
        """Instantiate a pipeline"""
        # Configure the builder

        assert hasattr(self, 'WF_ENTRYPOINT'), 'WF_ENTRYPOINT must be specified by the class'
        self.verbose = verbose
        # Initialise the preset
        self.preset_name = preset_name
        self.preset = PresetConfig.from_file(preset_name)
        self.protocol = protocol if protocol is not None else self.preset.default_protocol
        self.builder = None

    def get_builder(self, structure, code=None, protocol=None, overrides=None, **kwargs):
        """
        Generate builder base on a given structure and overrides (if supplied)
        """
        protocol = protocol or self.protocol
        overrides = overrides or {}
        code = code or self.preset.default_code
        options = kwargs.pop('options', {})
        options = recursive_merge(self.preset.get_code_specific_options(code, 'options'), options)

        builder = WorkflowFactory(self.WF_ENTRYPOINT).get_builder_from_protocol(
            code=orm.load_code(code),
            structure=structure,
            protocol=protocol,
            overrides=overrides,
            options=options,
            **kwargs,
        )
        self.builder = builder
        # Apply other settings
        self.set_settings(self.preset.get_code_specific_options(code, 'settings'))
        # Apply other settings
        self.set_incar(self.preset.get_code_specific_options(code, 'incar'))
        return builder

    @property
    def reference_structure(self):
        return self.builder.structure

    def set_incar(self, incar_updates=None, update_all=True, ports=None, namespace='incar', **kwargs):
        """
        Set incar dictionary
        """
        if incar_updates is None and not kwargs:
            return self

        if update_all:
            ports_nodes = recursive_search_dict_with_key(self.builder, 'incar')
        else:
            ports = ports or ['parameters']
            ports_nodes = [[port, self._get_port_node(port)] for port in ports]
        updates = deepcopy(incar_updates or {})
        updates.update(kwargs)
        for port, node in ports_nodes:
            self._update_dict_node(port, updates, dict_node=node, namespace=namespace)
        return self

    def set_options(self, option_updates=None, ports=None, update_all=True, **kwargs):
        """Set the options input port"""
        if option_updates is None and not kwargs:
            return self
        if update_all:
            calc_namespaces = []
            for port, namespace in recursive_search_port_basename(self.builder, 'calc'):
                if 'metadata' in namespace and 'options' in namespace['metadata']:
                    calc_namespaces.append([port, namespace])
        else:
            ports = ports or ['calc']
            calc_namespaces = [[port, self._get_port_node(port)] for port in ports]
        updates = option_updates or {}
        # Use recursive merge so existing nested keys will not be replaced
        updates = recursive_merge(updates, kwargs)
        # Update the options
        for port, namespace in calc_namespaces:
            # Here the port is only updated if the parent namespace is not empty or it is marked as 'required'
            # This is because `options`` is a special none-db port which may exist even if 'populate_defaults' is
            # set to False for namespaces that is optional. Otherwise, these optional namespace becomes 'defined'
            # , triggering its validation and then fails (as other 'required' fields are not defined inside the
            # namespace)
            if has_content(namespace) or namespace._port_namespace._required:
                namespace['metadata']['options'] = recursive_merge(dict(namespace['metadata']['options']), updates)
        return self

    def set_resources(self, resources_updates=None, ports=None, update_all=True, **kwargs):
        """Set the options input port"""
        if resources_updates is None and not kwargs:
            return self
        if update_all:
            calc_namespaces = []
            for port, namespace in recursive_search_port_basename(self.builder, 'calc'):
                if 'metadata' in namespace and 'options' in namespace['metadata']:
                    calc_namespaces.append([port, namespace])
        else:
            ports = ports or ['calc']
            calc_namespaces = [[port, self._get_port_node(port)] for port in ports]
        # Update the resources
        updates = deepcopy(resources_updates or {})
        updates.update(kwargs)
        for port, namespace in calc_namespaces:
            # Here the port is only updated if the parent namespace is not empty or it is marked as 'required'
            # This is because `options`` is a special none-db port which may exist even if 'populate_defaults' is
            # set to False for namespaces that is optional. Otherwise, these optional namespace becomes 'defined'
            # , triggering its validation and then fails (as other 'required' fields are not defined inside the
            # namespace)
            if has_content(namespace) or namespace._port_namespace._required:
                namespace['metadata']['options']['resources'].update(updates)
        return self

    def _update_ports_by_base_name(
        self, value, port_basename, ports=None, update_all=True, merge=False, skip_empty=True
    ):
        """Update a port by basename"""
        if update_all:
            port_and_nodes = recursive_search_port_basename(self.builder, port_basename)
        else:
            ports = ports or [port_basename]
            port_and_nodes = [[port, self._get_port_node(port)] for port in ports]
        # Update the options
        for port, node in port_and_nodes:
            if merge and isinstance(node, orm.Dict):
                self._set_node_to_port(port, update_dict_node(node, value))
            else:
                self._set_node_to_port(port, value)
        return self

    def _update_dict_node(self, port, update: dict, dict_node=None, namespace=None, reuse_if_possible=True):
        """ """
        if not update:
            return
        dict_node = dict_node or self._get_port_node(port)
        updated = update_dict_node(dict_node, update, namespace=namespace, reuse_if_possible=reuse_if_possible)
        self._set_node_to_port(port, updated)

    def _set_generic_port_by_dict(self, _port_name, value=None, ports=None, update_all=True, skip_empty=True, **kwargs):
        """Set a generic port by a value or kwargs"""
        if value is None and not kwargs:
            return self
        value = value or {}
        value = deepcopy(value)
        value.update(kwargs)
        self._update_ports_by_base_name(
            value, _port_name, ports=ports, update_all=update_all, skip_empty=skip_empty, merge=True
        )

    def _get_port_node(self, port):
        """Return the node corresponds to specific port"""
        parts = port.split('.')
        item = self.builder
        for part in parts:
            item = item.get(part)
        return item

    def _set_node_to_port(self, port, node: orm.Data):
        """Set a node to a specific port of hte builder"""
        if node is None:
            return
        parts = port.split('.')
        item = self.builder
        for part in parts[:-1]:
            item = item[part]
        setattr(item, parts[-1], node)

    def __repr__(self):
        string = f'{self.__class__.__name__}(protocol={self.protocol}, preset_name={self.preset_name})'
        if self.builder is not None:
            string += f'\nBuilder: {self.builder}'
        return string

    def _repr_pretty_(self, p, _=None) -> str:
        """Pretty representation for in the IPython console and notebooks."""

        string = f'{self.__class__.__name__}(protocol={self.protocol}, preset_name={self.preset_name})'
        p.text(string)
        if self.builder is not None:
            p.text('\nWith Builder:\n')
            self.builder._repr_pretty_(p, _)

    def set_kspacing(self, value, ports=None, update_all=True):
        """Update the kpoints spacing"""
        self._update_ports_by_base_name(value, 'kpoints_spacing', ports=ports, update_all=update_all)
        return self

    def set_kpoints_mesh(self, mesh: list[int], offset=(0.0, 0.0, 0.0), ports=None, update_all=True):
        """Set kpoints mesh"""
        kpoints = orm.KpointsData()
        kpoints.set_cell_from_structure(self.reference_structure)
        kpoints.set_kpoints_mesh(mesh, list(offset))
        self._update_ports_by_base_name(kpoints, 'kpoints', ports=ports, update_all=update_all)
        return self

    def set_label(self, label=None):
        """Alias to set the self.builder.metadata.label"""
        label = label or self.reference_structure.label
        self.builder.metadata.label = label
        return self

    def set_potential_family(self, value, ports=None, update_all=True):
        """Update the potential family"""
        self._update_ports_by_base_name(value, 'potential_family', ports=ports, update_all=update_all)
        return self

    def set_potential_mapping(self, value=None, ports=None, update_all=True, **kwargs):
        """Set the potential mapping"""
        self._set_generic_port_by_dict('potential_mapping', value=value, ports=ports, update_all=update_all, **kwargs)
        return self

    def set_code(self, value, ports=None, update_all=True):
        """Update the code node"""
        if isinstance(value, str):
            value = orm.load_code(value)
        self._update_ports_by_base_name(value, 'code', ports=ports, update_all=update_all)

    def set_settings(self, value, ports=None, update_all=True, **kwargs):
        """Update the `settings` port."""
        self._set_generic_port_by_dict('settings', value=value, ports=ports, update_all=update_all, **kwargs)

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
                    stdout = node.outputs.retrieved.get_object_content('vasp_output')
                    print(node, 'STDOUT:', stdout)
                    print(node, 'Retrieved files:', node.outputs.retrieved.list_object_names())
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


class VaspInputGenerator(BaseInputGenerator):
    """
    Updater for VaspWorkChain's builder
    """

    pass


class VaspRelaxInputGenerator(BaseInputGenerator):
    """
    Updater for VaspRelaxWorkChain's builder
    """

    WF_ENTRYPOINT = 'vasp.relax'

    def set_relax_settings(self, value=None, **kwargs):
        """Set the `relax_settings` port"""
        self._set_generic_port_by_dict('relax_settings', ports=['relax_settings'], value=value, **kwargs)
        return self

    def get_builder(self, structure, code=None, protocol=None, overrides=None, **kwargs):
        builder = super().get_builder(structure=structure, code=code, protocol=protocol, overrides=overrides, **kwargs)
        pdict = builder.vasp.parameters.get_dict()
        pdict['incar'].pop('nsw', None)
        pdict['incar'].pop('ibrion', None)
        pdict['incar'].pop('isif', None)
        # Case if the the parameters is stored
        if builder.vasp.parameters.is_stored:
            builder.vasp.parameters = pdict
        return builder


class VaspBandsInputGenerator(BaseInputGenerator):
    """
    Updater for VaspBandsWorkChain's builder
    """

    WF_ENTRYPOINT = 'vasp.bands'

    def set_band_settings(self, value=None, **kwargs):
        """Set the `band_settings` port"""
        self._set_generic_port_by_dict('band_settings', ports=['band_settings'], value=value, **kwargs)
        return self

    def set_settings(self, *args, **kwargs):
        """Set the settings port"""
        return super().set_settings(*args, ports=['scf.settings'], update_all=False, **kwargs)

    def get_builder(self, structure, code=None, protocol=None, overrides=None, run_relax=True, **kwargs):
        """
        Generate builder base on a given structure and overrides (if supplied)
        """
        protocol = protocol or self.protocol
        overrides = overrides or {}
        code = code or self.preset.default_code
        options = kwargs.pop('options', {})
        options = recursive_merge(self.preset.get_code_specific_options(code, 'options'), options)

        builder = WorkflowFactory(self.WF_ENTRYPOINT).get_builder_from_protocol(
            code=orm.load_code(code),
            structure=structure,
            protocol=protocol,
            overrides=overrides,
            options=options,
            run_relax=run_relax,
            **kwargs,
        )
        self.builder = builder
        # Apply other settings
        self.set_settings(self.preset.get_code_specific_options(code, 'settings'))
        # Apply other settings
        self.set_incar(self.preset.get_code_specific_options(code, 'incar'))
        return builder


class VaspConvergenceInputGenerator(BaseInputGenerator):
    """Updater for VaspConvergenceWorkChain"""

    WF_ENTRYPOINT = 'vasp.converge'

    def set_conv_settings(self, value=None, **kwargs):
        """Set the `conv_settings` port"""
        self._set_generic_port_by_dict('conv_settings', ports=['conv_settings'], value=value, **kwargs)
        return self


class VaspHybridBandsInputGenerator(VaspBandsInputGenerator):
    """Update for VaspHybridBandsWorkChain"""

    WF_ENTRYPOINT = 'vasp.hybrid_bands'


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
        left = dtmp.get(namespace, {})
    else:
        left = dtmp
    left = recursive_merge(left, content)
    # If namepsace is supplied, only update the target namespace inside the dict
    if namespace:
        dtmp[namespace] = left
    else:
        dtmp = left
    if node.is_stored:
        # There is no need to update the node if the content is the same as before
        if reuse_if_possible and dtmp == dtmp_backup:
            return node
        # The content is different, but the node is immutable, so we create a new node
        return orm.Dict(dict=dtmp)
    node.set_dict(dtmp)
    return node


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


def recursive_search_dict_with_key(namespace, search_key):
    """Recursively search for Dict node with certain key"""
    ports = []
    for port_key in namespace._valid_fields:
        value = namespace.get(port_key)
        if isinstance(value, orm.Dict):
            if search_key in value.get_dict():
                ports.append([port_key, value])
        if isinstance(value, ProcessBuilderNamespace):
            ports.extend(
                [
                    [port_key + '.' + sub_key, sub_value]
                    for sub_key, sub_value in recursive_search_dict_with_key(value, search_key)
                ]
            )
    return ports


def recursive_search_port_basename(namespace, basename):
    """Recursively search for Dict node with certain key"""
    ports = []
    for port_key in namespace._valid_fields:
        value = namespace.get(port_key)
        if port_key == basename:
            ports.append([port_key, value])
        if isinstance(value, ProcessBuilderNamespace):
            ports.extend(
                [
                    [port_key + '.' + sub_key, sub_value]
                    for sub_key, sub_value in recursive_search_port_basename(value, basename)
                ]
            )
    return ports


def has_content(mapping):
    """Check if a dictionary is all empty"""
    for key, value in mapping.items():
        if hasattr(value, 'items'):
            if has_content(value) is True:
                return True
        else:
            return True
    return False
