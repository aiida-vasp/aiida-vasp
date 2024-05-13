"""
Pipline for in place modification of builders

## Basic usage

```python
builder = VaspRelaxWorkChain.get_builder()
upd = VaspRelaxBuilder(builder)
upd.use_input_set(structure, ..., ...)
upd.use_code(...., ...)
upd.set_wallclock_seconds(...).set_num_machines(...)
```

## Using configuration dictionary

Instead of setting up the builder interactively, the updater can be initialised using
a dictionary.

```python
config_relax = {
    'overrides': {'encut': 520, 'ediff': 1e-6, 'ispin': 1, 'ncore': 2, 'kpar': 8},
    'code': 'vasp-6.3.0-std@mn',
    'options': {'max_wallclock_seconds': 3600 * 24 },
    'resources': {'tot_num_mpiprocs': 48 * 16, 'num_machines': 16},
    'inputset': 'UCLHSE06RelaxSet',
    'relax_settings': {'force_cutoff': 0.02},
    'kspacing': 0.07
}

upd = VaspRelaxWorkChain.init_from_config(structure, config_relax)
```
"""

from pprint import pprint
from typing import List, Union
from warnings import warn

from aiida import orm
from aiida.common.extendeddicts import AttributeDict
from aiida.engine.processes.builder import ProcessBuilder

from ..inputset.vaspsets import VASPInputSet
from ..relax import RelaxOptions
from .dictwrap import DictWrapper

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


class BuilderUpdater:
    """Base class for builder updater"""

    def __init__(self, builder: ProcessBuilder):
        """Instantiate a pipline"""
        self._builder = builder

    @property
    def builder(self):
        """The builder to be used for launching the calculation"""
        return self._builder

    def show_builder(self):
        """Print stuff defined in the builder"""
        pprint(builder_to_dict(self.builder, unpack=True))


DEFAULT_SET = 'UCLRelaxSet'


class VaspBuilderUpdater(BuilderUpdater):
    WF_ENTRYPOINT = 'vaspu.vasp'
    DEFAULT_SET = DEFAULT_SET

    def __init__(self, builder, root_namespace=None):
        """
        Initialise the update object.

        :param builder: The ``ProcessBuilder`` or ``ProcessBuilderNamespace`` to be used for setting
          standared VaspWorkChain inputs.

        :param root_namespace: The namespace to be assumed to be the *root*, e.g. where the input structure
          should be specified.
        """
        super().__init__(builder)
        if root_namespace is None:
            self.root_namespace = builder
        else:
            self.root_namespace = root_namespace

        self.namespace_vasp = builder
        self.parameters_wrapped = None
        self.options_wrapped = None
        self.settings_wrapped = None

    @property
    def reference_structure(self):
        """Reference structure used for setting kpoints and other stuff"""
        return self.root_namespace.structure

    @property
    def builder(self):
        """The builder to be used for launching the calculation"""
        return self.root_namespace

    def use_inputset(self, structure, set_name='UCLRelaxSet', overrides=None):
        inset = VASPInputSet(set_name, structure, overrides=overrides)
        self.namespace_vasp.parameters = orm.Dict(dict={'incar': inset.get_input_dict()})
        self.namespace_vasp.potential_family = orm.Str('PBE.54')
        self.namespace_vasp.potential_mapping = orm.Dict(dict=inset.get_pp_mapping())
        self.namespace_vasp.kpoints_spacing = orm.Float(0.05)
        self.root_namespace.structure = structure

        return self

    def _initialise_parameters_wrapper(self, force=False):
        """Initialise DictWrapper for tracking INCAR tags"""
        if self.parameters_wrapped is None or force:
            self.parameters_wrapped = DictWrapper(self.namespace_vasp.parameters, self.namespace_vasp, 'parameters')

    def _initialise_options_wrapper(self, force=False):
        """Initialise DictWrapper for tracking INCAR tags"""
        if self.options_wrapped is None or force:
            self.options_wrapped = DictWrapper(self.namespace_vasp.options, self.namespace_vasp, 'options')

    def _initialise_settings_wrapper(self, force=False):
        """Initialise DictWrapper for tracking INCAR tags"""
        if self.settings_wrapped is None or force:
            self.settings_wrapped = DictWrapper(self.namespace_vasp.settings, self.namespace_vasp, 'settings')

    def set_kspacing(self, kspacing: float):
        self.namespace_vasp.kpoints_spacing = orm.Float(kspacing)
        if self.namespace_vasp.kpoints:
            del self.namespace_vasp.kpoints
        return self

    update_kspacing = set_kspacing

    @property
    def incar(self):
        """Return the INCAR dictionary"""
        return dict(self.namespace_vasp.parameters['incar'])

    @property
    def settings(self):
        """Return the wrapped settings dictionary"""
        self._initialise_settings_wrapper()
        return self.settings_wrapped

    @property
    def parameters(self):
        """Return the wrapped parameters dictionary"""
        self._initialise_parameters_wrapper()
        return self.parameters_wrapped

    @property
    def options(self):
        """Return the wrapped options dictionary"""
        self._initialise_options_wrapper()
        return self.options_wrapped

    def set_code(self, code: Union[str, orm.Code]):
        if isinstance(code, str):
            self.namespace_vasp.code = orm.Code.get_from_string(code)
        else:
            self.namespace_vasp.code = code
        return self

    update_code = set_code

    def clear_incar(self):
        """Clear existing settings"""
        if self.namespace_vasp.parameters:
            del self.namespace_vasp.parameters
        self.parameters_wrapped = None
        return self

    def update_incar(self, *args, **kwargs):
        """Update incar tags"""
        if self.namespace_vasp.parameters is None:
            self.namespace_vasp.parameters = orm.Dict(dict={'incar': {}})

        self._initialise_parameters_wrapper()
        # Make a copy of the incar for modification
        incar = dict(self.parameters_wrapped['incar'])
        incar.update(*args, **kwargs)
        self.parameters_wrapped['incar'] = incar
        return self

    def set_default_options(self, **override):
        options = None

        # Try to use a sensible default from code's computer/scheduler type
        if self.namespace_vasp.code:
            computer = self.namespace_vasp.code.computer
            for key in OPTIONS_TEMPLATES:
                if key in computer.label.upper():
                    new = dict(OPTIONS_TEMPLATES[key])
                    new.update(override)
                    options = orm.Dict(dict=new)
                    break
            if options is None:
                for key in OPTIONS_TEMPLATES:
                    if key in computer.scheduler_type.upper():
                        new = dict(OPTIONS_TEMPLATES[key])
                        new.update(override)
                        options = orm.Dict(dict=new)
                        break

        # Use the very default settings
        if options is None:
            warn('Using default options template - adjustment needed for the target computer')
            options = orm.Dict(
                dict={
                    'resources': {
                        'tot_num_mpiprocs': 1,
                    },
                    'max_wallclock_seconds': 3600,
                    'import_sys_environment': False,
                    **override,
                }
            )

        self.namespace_vasp.options = options
        self._initialise_options_wrapper()
        return self

    def set_kpoints_mesh(self, mesh: List[int], offset: List[float]):
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

    update_kpoints_mesh = set_kpoints_mesh

    def update_options(self, *args, **kwargs):
        """Update the options"""
        if self.options_wrapped is None:
            self.set_default_options()
        self.options_wrapped.update(*args, **kwargs)
        return self

    def clear_options(self):
        if self.namespace_vasp.options:
            del self.namespace_vasp.options
        self.options_wrapped = None
        return self

    def update_settings(self, *args, **kwargs):
        """Update the settings"""
        if self.namespace_vasp.settings is None:
            self.namespace_vasp.settings = orm.Dict(dict={})
        self._initialise_settings_wrapper()
        self.settings_wrapped.update(*args, **kwargs)
        return self

    def clear_settings(self):
        """Clear existing settings"""
        if self.namespace_vasp.settings:
            del self.namespace_vasp.settings
        self.settings_wrapped = None
        return self

    def set_label(self, label=None):
        """Set the toplevel label, default to the label of the structure"""
        if label is None:
            label = self.root_namespace.structure.label
        self.root_namespace.metadata.label = label
        return self

    update_label = set_label

    def update_resources(self, *args, **kwargs):
        """Update resources"""
        if self.options_wrapped is None:
            self.set_default_options()
        resources = dict(self.options_wrapped['resources'])
        resources.update(*args, **kwargs)
        self.options_wrapped['resources'] = resources
        return self

    set_resources = update_resources

    def update_from_config(self, structure: orm.StructureData, config: dict):
        """
        Setup from a configuration dictionary

        The configuration dictionary should contain the following items:
          - inputset
          - overrides
          - code
          - options
          - resources
          - settings
          - kspacing
          - kmesh
          - potential_mapping
        """

        self.use_inputset(
            structure,
            config.get('inputset', self.DEFAULT_SET),
            overrides=config.get('overrides', {}),
        )
        self.set_code(orm.Code.get_from_string(config['code']))

        # Kpoint spacing/mesh
        if 'kspacing' in config:
            self.set_kspacing(config['kspacing'])
        if 'kmesh' in config:
            self.set_kpoints_mesh(*config['kmesh'])

        # Potential mapping
        if 'potential_mapping' in config:
            self.builder.potential_mapping = orm.Dict(dict=config['potential_mapping'])

        self.set_default_options(**config.get('options', {}))
        self.update_resources(**config.get('resources', {}))
        if 'settings' in config:
            self.update_settings(**config['settings'])
        self.set_label(f'{structure.label}')
        return self

    @classmethod
    def init_from_config(cls, structure, config):
        """Initialise from a configuration dictionary"""
        from aiida.plugins import WorkflowFactory

        upd = cls(WorkflowFactory(cls.WF_ENTRYPOINT).get_builder())
        return upd.update_from_config(structure, config)


class VaspNEBUpdater(VaspBuilderUpdater):
    WF_ENTRYPOINT = 'vasp.neb'

    def __init__(self, builder):
        """Initialise the builder updater"""
        super().__init__(builder)

    @property
    def reference_structure(self):
        """Return the reference structure"""
        return self.namespace_vasp.initial_structure

    def use_inputset(self, initial_structure, set_name='UCLRelaxSet', overrides=None):
        inset = VASPInputSet(set_name, initial_structure, overrides=overrides)
        self.namespace_vasp.parameters = orm.Dict(dict={'incar': inset.get_input_dict()})
        self.namespace_vasp.potential_family = orm.Str('PBE.54')
        self.namespace_vasp.potential_mapping = orm.Dict(dict=inset.get_pp_mapping())
        self.namespace_vasp.kpoints_spacing = orm.Float(0.05)
        self.namespace_vasp.initial_structure = initial_structure
        return self

    def set_final_structure(self, final_structure):
        self.namespace_vasp.final_structure = final_structure
        return self

    def set_neb_images(self, images):
        """Set the NEB images"""

        if isinstance(images, list):
            output = {f'image_{i:02d}': image for i, image in enumerate(images)}
        elif isinstance(images, (dict, AttributeDict)):
            output = images
        self.namespace_vasp.neb_images = output
        return self

    def set_interpolated_images(self, nimages):
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
        # Update the final image
        self.set_final_structure = interpolated['image_final']


class VaspRelaxUpdater(VaspBuilderUpdater):
    """
    An updater for VaspRelaxWorkChain
    """

    WF_ENTRYPOINT = 'vaspu.relax'

    def __init__(self, builder, override_vasp_namespace=None, namespace_relax=None):
        super().__init__(builder, root_namespace=builder)
        # The primary VASP namespace is under builder.vasp
        if override_vasp_namespace is None:
            self.namespace_vasp = builder.vasp
        else:
            self.namespace_vasp = override_vasp_namespace

        if namespace_relax is None:
            self.namespace_relax = builder
        else:
            self.namespace_relax = namespace_relax

    def update_relax_settings(self, **kwargs):
        """Set/update RelaxOptions controlling the operation of the workchain"""

        if self.namespace_relax.relax_settings is None:
            current_options = RelaxOptions()
        else:
            current_options = RelaxOptions(**self.namespace_relax.relax_settings.get_dict())
        for key, value in kwargs.items():
            setattr(current_options, key, value)
        self.namespace_relax.relax_settings = current_options.to_aiida_dict()
        return self

    def clear_relax_settings(self):
        """Reset any existing relax options"""
        self.namespace_relax.relax_settings = RelaxOptions().to_aiida_dict()
        return self

    def update_from_config(self, structure: orm.StructureData, config: dict):
        super().update_from_config(structure, config)
        self.update_relax_settings(**config.get('relax_settings', {}))
        return self


def builder_to_dict(builder, unpack=True):
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


class VaspConvUpdater(VaspBuilderUpdater):
    """Update for VaspConvergenceWorkChain"""

    WF_ENTRYPOINT = 'vaspu.converge'

    def update_from_config(self, structure: orm.StructureData, config: dict):
        """Update from a configuration dictionary"""
        super().update_from_config(structure, config)
        self.use_conv_settings(**config.get('conv_settings', {}))
        return self

    def use_conv_settings(self, **kwargs):
        """
        Use the supplied convergence settings
        """
        from ..converge import ConvOptions

        opts = ConvOptions(**kwargs)
        self.builder.conv_settings = opts.to_aiida_dict()


class VaspBandUpdater(VaspBuilderUpdater):
    """Updater for VaspBandsWorkChain"""

    WF_ENTRYPOINT = 'vaspu.bands'

    def __init__(self, builder, namespace_vasp=None):
        # The primary VASP namespace is under builder.vasp
        if namespace_vasp is None:
            super().__init__(builder.scf, root_namespace=builder)
        else:
            self.namespace_vasp = namespace_vasp

    def update_from_config(self, structure: orm.StructureData, config: dict):
        """
        Update the builder from a configuration dictionary.

        The dictionary must has a ``scf`` key holding the configurations for singlepoint
        calculations.
        The ``relax`` key is optional.
        """
        super().update_from_config(structure, config['scf'])

        # Specify the relaxation and NAC namespace
        if 'relax' in config:
            relax_upd = VaspRelaxUpdater(
                self.root_namespace,
                namespace_relax=self.root_namespace.relax,
                override_vasp_namespace=self.root_namespace.relax.vasp,
            )
            relax_upd.update_from_config(structure, config['relax'])
        return self

    def set_label(self, label=None):
        """Set the toplevel label, default to the label of the structure"""
        super().set_label(label)
        if label:
            if is_specified(self.root_namespace.relax):
                self.root_namespace.relax.metadata.label = label
            if is_specified(self.root_namespace.scf):
                self.root_namespace.scf.metadata.label = label
        return self


class VaspHybridBandUpdater(VaspBandUpdater):
    """Updater for VaspHybridBandsWorkChain"""

    WF_ENTRYPOINT = 'vaspu.hybrid_bands'

    def update_from_config(self, structure: orm.StructureData, config: dict):
        super().update_from_config(structure, config)

        self.builder.symprec = orm.Float(config['symprec'])
        self.builder.kpoints_per_split = orm.Int(config['kpoints_per_split'])
        return self


class VaspAutoPhononUpdater(VaspBuilderUpdater):
    """Updater for VaspAutoPhononWorkChain"""

    WF_ENTRYPOINT = 'vaspu.phonopy'

    def __init__(self, builder: ProcessBuilder):
        """Initialise with an existing ProcessBuilder for VaspAutoPhononWorkChain"""
        super().__init__(builder.singlepoint, root_namespace=builder)

    def set_phonon_settings(self, options):
        """
        Update the phonon-related options

        example::

          {
            'primitive_matrix': 'auto',
            'supercell_matrix': [2, 2, 2],    # Supercell matrix
            'mesh': 30,                       # Mesh for DOS/PDOS/thermal properties
          }


        """
        self.root_namespace.phonon_settings = orm.Dict(options)
        return self

    def update_from_config(self, structure: orm.StructureData, config: dict):
        """
        Update the builder from a configuration dictionary.

        The dictionary must has a ``singlepoint`` key holding the configurations for singlepoint
        calculations, and a ``phonon_options`` for Phonopy options to be used.
        The ``relax`` and ``nac`` keys are optional.
        """

        super().update_from_config(structure, config['singlepoint'])

        # Specify the relaxation and NAC namespace
        if 'relax' in config:
            relax_upd = VaspRelaxUpdater(
                self.root_namespace,
                namespace_relax=self.root_namespace.relax,
                override_vasp_namespace=self.root_namespace.relax.vasp,
            )
            relax_upd.update_from_config(structure, config['relax'])

        if 'nac' in config:
            nac_upd = VaspBuilderUpdater(self.root_namespace.nac, root_namespace=self.root_namespace)
            nac_upd.update_from_config(structure, config['nac'])

        # Update the phonon settings
        self.set_phonon_settings(config['phonon_settings'])
        return self

    def set_kpoints_mesh(self, mesh, offset) -> None:
        """Use mesh for kpoints"""
        kpoints = orm.KpointsData()
        # Use the reference supercell structure
        kpoints.set_cell_from_structure(self.reference_structure)
        kpoints.set_kpoints_mesh(mesh, offset)
        self.namespace_vasp.kpoints = kpoints
        if self.namespace_vasp.kpoints_spacing:
            del self.namespace_vasp.kpoints_spacing
        return self

    def _get_singlepoint_supercell(self) -> orm.StructureData:
        """Obtain the supercell for the singlepoint calculation"""
        import numpy as np
        from ase.build import make_supercell

        ref = self.root_namespace.structure.get_ase()

        # The sueprcell matrix should be a vector or a matrix
        mat = np.array(self.root_namespace.phonon_settings['supercell_matrix'])
        if mat.size == 3:
            mat = np.diag(mat)

        # Convention of phonopy - the supercell matrix is the transpose of that would be used
        # for ase
        return orm.StructureData(ase=make_supercell(ref, mat.T))

    def show_builder(self):
        """Print stuff defined in the builder"""
        pprint(builder_to_dict(self.root_namespace, unpack=True))


def is_specified(port_namespace):
    """Check if there is anything specified under a PortNamespace"""
    return any(map(bool, port_namespace.values()))
