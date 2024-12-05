"""
One liner input generator for aiida-vasp
"""

import logging
from copy import deepcopy
from dataclasses import dataclass, field
from pathlib import Path
from typing import List, Optional, Union
from warnings import warn

from aiida import orm
from aiida.common.extendeddicts import AttributeDict
from aiida.engine.processes.builder import ProcessBuilder, ProcessBuilderNamespace
from yaml import safe_load

from ..inputset.base import convert_lowercase
from ..inputset.vaspsets import VASPInputSet
from ..relax import RelaxOptions

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


def get_library_path():
    """Get the path where the YAML files are stored within this package"""
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
    def from_file(cls, fname):
        """Load from file"""

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

    def get_code_specific_options(self, code: str, namespace: str) -> dict:
        """Return code specific options, if exists"""
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
        preset_name: Union[None, str] = None,
        builder: Union[ProcessBuilder, None] = None,
        verbose=False,
        inputset=None,
    ):
        """Instantiate a pipeline"""
        # Configure the builder
        from aiida.plugins import WorkflowFactory

        assert hasattr(self, 'WF_ENTRYPOINT'), 'WF_ENTRYPOINT must be specified by the class'
        self.verbose = verbose
        if builder is None:
            builder = WorkflowFactory(self.WF_ENTRYPOINT).get_builder()
        self._builder = builder
        if preset_name is None:
            preset_name = DEFAULT_PRESET
        self.preset_name = preset_name
        self.preset = VaspPresetConfig.from_file(preset_name)
        self.inputset = inputset if inputset is not None else self.preset.inputset

    @property
    def builder(self) -> ProcessBuilder:
        """The builder to be used for launching the calculation"""
        return self._builder

    def submit(self) -> orm.WorkChainNode:
        """Submit the workflow to the daemon and return the workchain node"""
        from aiida.engine import submit

        return submit(self.builder)

    def run_get_node(self, verbose=True) -> orm.WorkChainNode:
        """Run the workflow with the current python process"""
        from aiida.engine import run_get_node

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

    def _get_help(self, namespace: str, print_to_stdout=True, inout='inputs'):
        """
        Return the help message for a given namespace
        The `.` syntax for the namespace is supported.
        """
        levels = namespace.split('.')
        data_dict = self.builder._process_spec.get_description()[inout]
        for key in levels:
            data_dict = data_dict[key]

        if print_to_stdout is True:
            print(data_dict.get('help', 'No help message information found'))
        else:
            return data_dict.get('help', 'No help message information found')

    def get_output_help(self, namespace: str, print_to_stdout=True):
        """Return the help message for a given namespace"""
        self._get_help(namespace, print_to_stdout=print_to_stdout, inout='outputs')

    def get_input_help(self, namespace: str, print_to_stdout=True):
        """Return the help message for a given namespace"""
        self._get_help(namespace, print_to_stdout=print_to_stdout, inout='inputs')


class VaspBuilderUpdater(BaseBuilderUpdater):
    WF_ENTRYPOINT = 'vasp.v2.vasp'
    DEFAULT_INPUTSET = DEFAULT_INPUTSET

    def __init__(
        self,
        preset_name: Optional[str] = None,
        builder: Optional[ProcessBuilder] = None,
        root_namespace: Optional[str] = None,
        code: Optional[str] = None,
        verbose: bool = False,
    ):
        """
        Initialise the update object.

        :param builder: The ``ProcessBuilder`` or ``ProcessBuilderNamespace`` to be used for setting
          standared VaspWorkChain inputs.

        :param root_namespace: The namespace to be assumed to be the *root*, e.g. where the input structure
          should be specified.
        """
        super().__init__(preset_name=preset_name, builder=builder, verbose=verbose)
        # Define the root namespace - e.g. the VaspWorkChain namespace where structure should be specified
        if root_namespace is None:
            self.root_namespace = self._builder
        else:
            self.root_namespace = root_namespace

        self.namespace_vasp = self._builder
        # Defult
        self.inputset = None
        self.code = self.preset.default_code if code is None else code

    @property
    def reference_structure(self) -> orm.StructureData:
        """Reference structure used for setting kpoints and other stuff"""
        return self.root_namespace.structure

    def clear(self) -> None:
        """Clear the nodes set to the namespace"""
        self.namespace_vasp.parameters = None
        self.namespace_vasp.options = None
        self.namespace_vasp.settings = None
        self.namespace_vasp.kpoints = None
        self.namespace_vasp.potential_family = None
        self.namespace_vasp.potential_mapping = None

        self.root_namespace.structure = None
        self.root_namespace.metadata.label = None

    def apply_preset(self, initial_structure, code=None, label=None, overrides=None) -> 'VaspBuilderUpdater':
        """
        Apply the preset
        """
        if code is None:
            code = self.code
            logging.info(f'Using code {code}')
        self.use_inputset(
            initial_structure,
            set_name=self.inputset,
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
        structure,
        set_name=None,
        overrides=None,
        apply_preset=False,
        code=None,
        structure_node_name='structure',
    ) -> 'VaspBuilderUpdater':
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

        inset = VASPInputSet(set_name, overrides=overrides_, verbose=self.verbose)
        self.namespace_vasp.parameters = orm.Dict(dict={'incar': inset.get_input_dict(structure)})
        self.namespace_vasp.potential_family = orm.Str(inset.get_potcar_family())
        self.namespace_vasp.potential_mapping = orm.Dict(dict=inset.get_pp_mapping(structure))
        self.namespace_vasp.kpoints_spacing = orm.Float(inset.get_kpoints_spacing())
        setattr(self.root_namespace, structure_node_name, structure)
        self.inputset = inset
        return self

    def set_kspacing(self, kspacing: float) -> 'VaspBuilderUpdater':
        self.namespace_vasp.kpoints_spacing = orm.Float(kspacing)
        if self.namespace_vasp.kpoints:
            del self.namespace_vasp.kpoints
        return self

    update_kspacing = set_kspacing

    @property
    def parameters(self) -> Union[orm.Dict, None]:
        """Return the parameters node"""
        return self.namespace_vasp.parameters

    @property
    def settings(self) -> Union[orm.Dict, None]:
        """Return the wrapped settings node"""
        return self.namespace_vasp.settings

    @property
    def options(self) -> Union[orm.Dict, None]:
        """Return the wrapped options node"""
        return self.namespace_vasp.options

    def set_code(self, code: Union[str, orm.Code, None] = None) -> 'VaspBuilderUpdater':
        if code is None:
            code = self.preset.default_code
        if isinstance(code, str):
            code = orm.load_code(code)

        self.namespace_vasp.code = code
        return self

    def update_code(self, code: Union[str, orm.Code]):
        warn('update_code is deprecated, use set_code instead', DeprecationWarning)
        return self.set_code(code)

    def set_incar(self, *args, **kwargs) -> 'VaspBuilderUpdater':
        """Update incar tags"""
        if self.namespace_vasp.parameters is None:
            self.namespace_vasp.parameters = orm.Dict(dict={'incar': {}})
        content = dict(*args, **kwargs)
        node = update_dict_node(self.namespace_vasp.parameters, content, 'incar')
        self.namespace_vasp.parameters = node
        return self

    def update_incar(self, *args, **kwargs) -> 'VaspBuilderUpdater':
        warn('update_incar is deprecated, use set_incar instead', DeprecationWarning)
        return self.set_incar(*args, **kwargs)

    def set_options(
        self, *args, code: Optional[str] = None, apply_preset: bool = False, **kwargs
    ) -> 'VaspBuilderUpdater':
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

    def update_options(self, *args, **kwargs) -> 'VaspBuilderUpdater':
        warn('update_options is deprecated, use set_options instead', DeprecationWarning)
        return self.set_options(*args, **kwargs)

    def set_kpoints_mesh(self, mesh: List[int], offset: List[float]) -> 'VaspBuilderUpdater':
        """Use mesh for kpoints"""
        kpoints = orm.KpointsData()
        kpoints.set_cell_from_structure(self.reference_structure)
        kpoints.set_kpoints_mesh(mesh, offset)
        self.namespace_vasp.kpoints = kpoints
        try:
            del self.namespace_vasp.kpoints_spacing
        except KeyError:
            pass
        return self

    def update_kpoints_mesh(self, mesh: List[int], offset: List[float]) -> 'VaspBuilderUpdater':
        warn('update_kpoints_mesh is deprecated, use set_kpoints_mesh instead', DeprecationWarning)
        return self.set_kpoints_mesh(mesh, offset)

    def set_settings(
        self, *args, code: Optional[str] = None, apply_preset: bool = False, **kwargs
    ) -> 'VaspBuilderUpdater':
        """Update the settings"""

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

    def update_settings(self, *args, **kwargs):
        warn('update_settings is deprecated, use set_settings instead', DeprecationWarning)
        return self.set_settings(*args, **kwargs)

    def set_label(self, label: Optional[str] = None) -> 'VaspBuilderUpdater':
        """Set the toplevel label, default to the label of the structure"""
        if label is None:
            label = self.root_namespace.structure.label
        self.root_namespace.metadata.label = label
        return self

    def update_label(self, label=None) -> 'VaspBuilderUpdater':
        warn('update_label is deprecated, use set_label instead', DeprecationWarning)
        return self.set_label(label)

    def set_resources(self, *args, **kwargs) -> 'VaspBuilderUpdater':
        """Update resources"""
        if self.namespace_vasp.options is None:
            raise RuntimeError('Please set the options before setting resources')
        resources = dict(self.namespace_vasp.options['resources'])
        resources.update(*args, **kwargs)
        self.namespace_vasp.options = update_dict_node(self.namespace_vasp.options, resources, 'resources')
        return self

    def update_resources(self, *args, **kwargs) -> 'VaspBuilderUpdater':
        warn('update_resources is deprecated, use set_resources instead', DeprecationWarning)
        return self.set_resources(*args, **kwargs)

    def _set_options(
        self, option_class, option_name: str, target_namespace: Union[ProcessBuilder, ProcessBuilderNamespace], **kwargs
    ):
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
        """Return the reference structure"""
        return self.namespace_vasp.initial_structure

    def apply_preset(
        self,
        structure_init,
        structure_final,
        code=None,
        label=None,
        interpolate=True,
        nimages=5,
        **kwargs,
    ) -> 'VaspNEBUpdater':
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
        set_name=None,
        overrides=None,
        apply_preset=False,
        code=None,
    ) -> 'VaspNEBUpdater':
        super().use_inputset(
            structure=initial_structure,
            set_name=set_name,
            overrides=overrides,
            apply_preset=apply_preset,
            code=code,
            structure_node_name='initial_structure',
        )
        return self

    def set_label(self, label: Optional[str] = None) -> 'VaspNEBUpdater':
        """Set the toplevel label, default to the label of the structure"""
        if label is None:
            label = self.root_namespace.initial_structure.label
        self.root_namespace.metadata.label = label
        return self

    def set_final_structure(self, final_structure: orm.StructureData) -> 'VaspNEBUpdater':
        self.namespace_vasp.final_structure = final_structure
        return self

    def set_neb_images(self, images: Union[list, dict, AttributeDict]) -> 'VaspNEBUpdater':
        """Set the NEB images"""

        if isinstance(images, list):
            output = {f'image_{i:02d}': image for i, image in enumerate(images)}
        elif isinstance(images, (dict, AttributeDict)):
            output = images
        self.namespace_vasp.neb_images = output
        return self

    def set_interpolated_images(self, nimages: int) -> 'VaspNEBUpdater':
        """
        Interpolate images and set as inputs structures

        This requires the initial and final structure to be set already.
        This function also update the final image with PBC issue fixed.
        """
        from .transform import neb_interpolate

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

    def view_images(self):
        """
        Visualize the images using ASE
        """
        from ase.visualize import view

        view(
            map(
                lambda x: x.get_ase(),
                [self.builder.initial_structure, *self.builder.neb_images.values(), self.builder.final_structure],
            )
        )


class VaspRelaxUpdater(VaspBuilderUpdater):
    """
    An updater for VaspRelaxWorkChain
    """

    WF_ENTRYPOINT = 'vasp.v2.relax'

    def __init__(
        self,
        preset_name: Optional[str] = None,
        builder: Optional[ProcessBuilder] = None,
        override_vasp_namespace: Optional[ProcessBuilderNamespace] = None,
        namespace_relax: Optional[ProcessBuilderNamespace] = None,
        code: Optional[str] = None,
    ):
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

    def apply_preset(
        self,
        structure: orm.StructureData,
        code: Optional[str] = None,
        label: Optional[str] = None,
        **kwargs,
    ) -> 'VaspRelaxUpdater':
        out = super().apply_preset(structure, code, label, **kwargs)
        self.set_relax_settings()
        return out

    def set_relax_settings(self, **kwargs) -> 'VaspRelaxUpdater':
        """Set/update RelaxOptions controlling the operation of the workchain"""
        self._set_options(RelaxOptions, 'relax_settings', self.namespace_relax, **kwargs)
        return self

    update_relax_settings = set_relax_settings

    def clear_relax_settings(self) -> 'VaspRelaxUpdater':
        """Reset any existing relax options"""
        self.namespace_relax.relax_settings = RelaxOptions().aiida_dict()
        return self

    def clear(self) -> 'VaspRelaxUpdater':
        super().clear()
        self.clear_relax_settings()
        return self


class VaspConvUpdater(VaspBuilderUpdater):
    """Update for VaspConvergenceWorkChain"""

    WF_ENTRYPOINT = 'vasp.v2.converge'

    def apply_preset(self, initial_structure, code=None, label=None, **kwargs) -> VaspBuilderUpdater:
        super().apply_preset(initial_structure, code, label, **kwargs)
        self.set_conv_settings()
        return self

    def set_conv_settings(self, **kwargs) -> 'VaspConvUpdater':
        """
        Use the supplied convergence settings
        """
        from ..converge import ConvOptions

        self._set_options(ConvOptions, 'conv_settings', self.builder, **kwargs)
        return self


class VaspBandUpdater(VaspBuilderUpdater):
    """Updater for VaspBandsWorkChain"""

    WF_ENTRYPOINT = 'vasp.v2.bands'

    def __init__(self, preset_name=None, builder=None, override_vasp_namespace=None, code=None):
        super().__init__(preset_name=preset_name, builder=builder, code=code, root_namespace=builder)
        # The primary VASP namespace is under builder.vasp
        if override_vasp_namespace is None:
            self.namespace_vasp = self.builder.scf
        else:
            self.namespace_vasp = override_vasp_namespace

    def get_relax_updater(self):
        """
        Return the relax updater for this band structure calculation

        The relax updater can be used to populate the `.relax` namespace which will
        trigger the relaxation of the structure.
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
        self, structure: orm.StructureData, run_relax: bool = False, label=None, **kwargs
    ) -> 'VaspBandUpdater':
        super().apply_preset(structure, label=label, **kwargs)

        # Specify the relaxation and NAC namespace
        if run_relax:
            relax_upd = self.get_relax_updater()
            relax_upd.apply_preset(structure, label=label, **kwargs)
        self.set_band_settings()
        return self

    def set_band_settings(self, **kwargs) -> 'VaspBandUpdater':
        from aiida_vasp.workchains.v2.bands import BandOptions

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
    """Check if there is anything specified under a PortNamespace"""
    return any(map(bool, port_namespace.values()))


def update_dict_node(
    node: orm.Dict,
    content: dict,
    namespace: Optional[str] = None,
    reuse_if_possible: bool = True,
) -> orm.Dict:
    """
    Update a Dict node with the content
    Optionally update an item of the Dict node.
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


def builder_to_dict(builder: ProcessBuilder, unpack: bool = True) -> dict:
    """
    Convert a builder to a dictionary and unpack certain nodes.

    When unpacked, the resulting dictionary cannot be used for `submit`/`run`.

    The primary useage of the resulting dictionary is for pretty printing.
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
