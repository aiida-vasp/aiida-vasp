"""
VASP Convergence WorkChain for AiiDA
====================================

This module provides a workchain for performing convergence tests in VASP calculations using the AiiDA workflow engine.
It automates the process of running multiple VASP calculations with varying plane-wave cutoff energies and
k-point spacings, collects the results, and provides utilities for analyzing and plotting the convergence data.
The pandas package is needed for generated DataFrames of the convergence test results.


Main Features
-------------

- Defines the `VaspConvergenceWorkChain`, which orchestrates a series of VASP calculations to test convergence with
  respect to:
    - Plane-wave cutoff energy (`ENCUT`)
    - K-point spacing
- Collects and summarizes output data such as total energy, maximum force, maximum stress, and magnetization.
- Provides helper functions to extract convergence data as pandas DataFrames and to plot the results.
- Includes a utility to generate a builder for the convergence workchain with user-specified protocols and options.

Inputs
------

- Structure and calculation parameters as required by the underlying VASP workchain.
- `conv_settings`: A dictionary specifying the convergence test ranges and steps for cutoff energy and k-point spacing.

Outputs
-------

- `cutoff_conv_data`: Summary of results for cutoff energy convergence.
- `kpoints_conv_data`: Summary of results for k-point spacing convergence.

Usage
-----

This module is intended to be used as part of an AiiDA plugin for VASP. Users can launch the convergence workchain
by providing the required inputs and analyze the results using the provided utility functions.


"""

from __future__ import annotations

from typing import Any

import numpy as np
from aiida import orm
from aiida.engine import ProcessSpec, WorkChain, append_, calcfunction
from aiida.plugins import WorkflowFactory

from aiida_vasp.utils.extended_dicts import update_nested_dict_node
from aiida_vasp.utils.opthold import ConvOptions

from .common.builder_updater import VaspBuilderUpdater
from .mixins import WithBuilderUpdater

# pylint:disable=no-member,unused-argument,no-self-argument,import-outside-toplevel


class VaspConvergenceWorkChain(WorkChain, WithBuilderUpdater):
    """
    A workchain to perform convergence tests.

    The inputs are essentially the same as for ``VaspWorChain`` but instead of launching
    a single calculation it launches many calculations with different kpoint spacing
    and the cut off energy.

    A ``conv_setting`` input controls the range of cut off energies and kpoint spacings.
    The available options are:

    - cutoff_start
    - cutoff_stop
    - cutoff_step
    - kspacing_start
    - kspacing_stop
    - kspacing_step
    - cutoff_kconv : cut-off energy for the kpoints convergence tests.
    - kspacing_cutconv : the kpoint spacing to be used for cut-off energy convergence tests.

    The the output data are collected and stored in two ``Dict`` output nodes.
    """

    _sub_workchain_string = 'vasp.v2.vasp'
    _sub_workchain = WorkflowFactory(_sub_workchain_string)
    ENERGY_KEY = 'energy_extrapolated'
    option_class = ConvOptions

    @classmethod
    def define(cls, spec: ProcessSpec) -> None:
        super().define(spec)

        spec.expose_inputs(cls._sub_workchain)
        spec.input(
            'conv_settings',
            help=ConvOptions.aiida_description(),
            validator=ConvOptions.aiida_validate,
            serializer=ConvOptions.aiida_serialize,
            valid_type=orm.Dict,
        )
        spec.outline(cls.setup, cls.launch_conv_calcs, cls.analyse)

        spec.exit_code(
            401,
            'ERROR_SUBWORKFLOW_ERRORED',
            message='At leaste one of the launched sub-workchain has failed',
        )
        spec.output('kpoints_conv_data', required=False)
        spec.output('cutoff_conv_data', required=False)

    def setup(self) -> None:
        """Setup the convergence workflow"""
        settings = self.inputs.conv_settings.get_dict()
        self.ctx.settings = settings

        # Planewave cut off energies
        start = settings['cutoff_start']
        stop = settings['cutoff_stop']
        if start < stop:
            cutoff_list = [start]
            cut = start
            # Ensure start and stop are always included
            while True:
                cut += settings['cutoff_step']
                if cut < stop:
                    cutoff_list.append(cut)
                else:
                    cutoff_list.append(stop)
                    break
        else:
            # Start is equal or larger than stop - signalling no need to do the test
            cutoff_list = []

        # Same treatment for kspacing
        start = settings['kspacing_start']
        stop = settings['kspacing_stop']

        if start > stop:
            spacing = start
            kspacing_list = [spacing]
            while True:
                spacing -= settings['kspacing_step']
                if spacing > settings['kspacing_stop']:
                    kspacing_list.append(spacing)
                else:
                    kspacing_list.append(settings['kspacing_stop'])
                    break
        else:
            kspacing_list = []

        self.ctx.cutoff_list = cutoff_list
        self.ctx.kspacing_list = kspacing_list

    def launch_conv_calcs(self) -> None:
        """
        Setup and launch the convergence calculations
        """
        if not self.ctx.cutoff_list:
            cut_k = 400  # Default if not supplied
        else:
            cut_k = min(self.ctx.cutoff_list)
        if not self.ctx.kspacing_list:
            k_cut = 0.06  # Default if not supplied
        else:
            k_cut = min(self.ctx.kspacing_list)

        cutoff_for_kconv = self.ctx.settings.get('cutoff_kconv', cut_k)
        kspacing_for_cutoffconv = orm.Float(self.ctx.settings.get('kspacing_cutconv', k_cut))

        # Launch cut off energy tests
        inputs = self.exposed_inputs(self._sub_workchain)
        inputs.kpoints_spacing = kspacing_for_cutoffconv
        original_label = inputs.metadata.get('label', '')
        for cut in self.ctx.cutoff_list:
            new_param = update_nested_dict_node(inputs.parameters, {'incar': {'encut': cut}})
            inputs.parameters = new_param
            if original_label:
                inputs.metadata.label = original_label + f' CUTCONV {cut:.2f}'
            else:
                inputs.metadata.label = f'CUTCONV {cut:.2f}'

            running = self.submit(self._sub_workchain, **inputs)
            self.report(f'Submitted {running} with cut off energy {cut:.1f} eV.')
            self.to_context(cutoff_conv_workchains=append_(running))

        # Launch kpoints convergence tests
        new_param = update_nested_dict_node(inputs.parameters, {'incar': {'encut': cutoff_for_kconv}})
        for kspacing in self.ctx.kspacing_list:
            inputs.parameters = new_param
            inputs.kpoints_spacing = kspacing
            if original_label:
                inputs.metadata.label = original_label + f' KCONV {kspacing:.3f}'
            else:
                inputs.metadata.label = f'KCONV {kspacing:.3f}'

            running = self.submit(self._sub_workchain, **inputs)
            self.report(f'Submitted {running} with kpoints spacing {kspacing:.3f}.')
            self.to_context(kpoints_conv_workchains=append_(running))

    def analyse(self) -> None:
        """
        Analyse the output of the calculations.
        Collect data to be plotted/analysed against the cut off energy and kpoints spacing
        """

        def get_maximum(forces: np.ndarray | None) -> float | None:
            if forces is None:
                return None
            norm = np.linalg.norm(forces, axis=1)
            return np.amax(norm)

        def collect_data(workchain: orm.WorkChainNode, energy_key: str) -> dict[str, Any]:
            """Collect the data from workchain output"""
            output = workchain.outputs.misc.get_dict()
            data = {}
            data['maximum_force'] = get_maximum(output.get('forces'))
            # Extract the magnetization
            magnetization = output.get('magnetization')
            if magnetization:
                data['magnetization'] = magnetization[0]
            data['maximum_stress'] = get_maximum(output.get('stress'))
            data['energy'] = output['total_energies'][energy_key]
            return data

        def unpack(name: str, input_data: dict[Any, Any]) -> dict[str, list[Any]]:
            """Unpack a dict with numberical keys"""
            output_dict = {name: []}
            for key, data in input_data.items():
                output_dict[name].append(key)
                # Append values to the corresponding lists
                for key_, value in data.items():
                    if key_ not in output_dict:
                        output_dict[key_] = []
                    output_dict[key_].append(value)
            return output_dict

        exit_code = None

        cutoff_data = {}
        cutoff_miscs = {}
        energy_key = None
        if 'cutoff_conv_workchains' in self.ctx:
            for iwork, workchain in enumerate(self.ctx.cutoff_conv_workchains):
                if workchain.exit_status != 0:
                    exit_code = self.exit_codes.ERROR_SUBWORKFLOW_ERRORED
                    self.report(f'Skipping workchain {workchain} with exit status {workchain.exit_status} ')
                    continue

                # Setup the energy key from the first workchain
                if not energy_key:
                    energy_key = next(iter(workchain.outputs.misc.get_dict()['total_energies'].keys()))

                cutoff = workchain.inputs.parameters['incar']['encut']
                cutoff_data[cutoff] = collect_data(workchain, energy_key)
                cutoff_data[cutoff]['mesh'] = workchain.called[0].inputs.kpoints.get_kpoints_mesh()[0]
                cutoff_miscs[f'worchain_{iwork}'] = workchain.outputs.misc

        kspacing_data = {}
        kspacing_miscs = {}

        if 'kpoints_conv_workchains' in self.ctx:
            for iwork, workchain in enumerate(self.ctx.kpoints_conv_workchains):
                if workchain.exit_status != 0:
                    exit_code = self.exit_codes.ERROR_SUBWORKFLOW_ERRORED
                    self.report(f'Skipping Workchain {workchain} with exit status {workchain.exit_status} ')
                    continue

                # Setup the energy key from the first workchain
                if not energy_key:
                    energy_key = next(iter(workchain.outputs.misc.get_dict()['total_energies'].values()))

                spacing = float(workchain.inputs.kpoints_spacing)
                kspacing_data[spacing] = collect_data(workchain, energy_key)
                kspacing_data[spacing]['mesh'] = workchain.called[0].inputs.kpoints.get_kpoints_mesh()[0]
                kspacing_miscs[f'worchain_{iwork}'] = workchain.outputs.misc
                cutoff = workchain.inputs.parameters['incar']['encut']
                kspacing_data[spacing]['cutoff_energy'] = cutoff

        # Calcfunction to link with the calculation output to the summary data node
        @calcfunction
        def create_links_kconv(**miscs: Any) -> orm.Dict:
            """Alias calcfunction to link summary node with miscs"""
            return orm.Dict(dict=unpack('kpoints_spacing', kspacing_data))

        @calcfunction
        def create_links_cutconv(**miscs: Any) -> orm.Dict:
            """Alias calcfunction to link summary node with miscs"""
            return orm.Dict(dict=unpack('cutoff_energy', cutoff_data))

        if kspacing_data:
            self.out('kpoints_conv_data', create_links_kconv(**kspacing_miscs))
        if cutoff_data:
            self.out('cutoff_conv_data', create_links_cutconv(**cutoff_miscs))

        return exit_code

    @staticmethod
    def get_conv_data(conv_work: orm.WorkChainNode, plot: bool = False, **plot_kwargs: Any) -> tuple[Any, Any]:
        """
        Convenient method for extracting convergence data

        Args:
        conv_work (orm.WorkChainNode): Convergence workflow node

        Returns:
        A tuple of cut-off convergence and k-point convergence result dataframe
        """
        cdf, kdf = get_conv_data(conv_work)
        if plot is True:
            plot_conv_data(cdf, kdf, **plot_kwargs)
        return cdf, kdf


def get_conv_data(conv_work: orm.WorkChainNode) -> tuple[Any, Any]:
    """
    Convenient method for extracting convergence data

    Args:
      conv_work (orm.WorkChainNode): Convergence workflow node

    Returns:
      A tuple of cut-off convergence and k-point convergence result data frame
    """
    import pandas as pd  # noqa: PLC0415

    if 'cutoff_conv_data' in conv_work.outputs:
        cutdf = pd.DataFrame(conv_work.outputs.cutoff_conv_data.get_dict())
        cutdf['energy_per_atom'] = cutdf['energy'] / len(conv_work.inputs.structure.sites)
        cutdf['dE_per_atom'] = cutdf['energy_per_atom'] - cutdf['energy_per_atom'].iloc[-1]
    else:
        cutdf = None

    if 'kpoints_conv_data' in conv_work.outputs:
        kdf = pd.DataFrame(conv_work.outputs.kpoints_conv_data.get_dict())
        kdf['energy_per_atom'] = kdf['energy'] / len(conv_work.inputs.structure.sites)
        kdf['dE_per_atom'] = kdf['energy_per_atom'] - kdf['energy_per_atom'].iloc[-1]
    else:
        kdf = None
    return cutdf, kdf


def plot_conv_data(cdf: Any, kdf: Any, **kwargs: Any) -> list[Any]:
    """
    Make two combined plots for the convergence test results.
    """

    import matplotlib.pyplot as plt  # noqa: PLC0415

    # Create a subplot
    figs = []
    if cdf is not None:
        fig, axs = plt.subplots(3, 1, sharex=True, **kwargs)
        figs.append(fig)
        axs[0].plot(cdf.cutoff_energy, cdf.dE_per_atom, '-x')
        axs[0].set_ylabel('dE (eV / atom)')
        i = 0
        if 'maximum_force' in cdf.columns:
            i += 1
            axs[i].plot(cdf.cutoff_energy, cdf.maximum_force, '-x')
            axs[i].set_ylabel(r'$F_{max}$ (eV$\AA^{-1}$)')
        if 'maximum_stress' in cdf.columns:
            i += 1
            axs[i].plot(cdf.cutoff_energy, cdf.maximum_stress, '-x')
            axs[i].set_ylabel(r'$S_{max}$ (kBar)')
        axs[i].set_xlabel('Cut-off energy (eV)')
        fig.tight_layout()

    if kdf is not None:
        fig, axs = plt.subplots(3, 1, sharex=True, **kwargs)
        figs.append(fig)
        axs[0].plot(kdf.kpoints_spacing, kdf.dE_per_atom, '-x')
        axs[0].set_ylabel('dE (eV / atom)')
        i = 0
        if 'maximum_force' in kdf.columns:
            i += 1
            axs[i].plot(kdf.kpoints_spacing, kdf.maximum_force, '-x')
            axs[i].set_ylabel(r'$F_{max}$ (eV$\AA^{-1}$)')
        if 'maximum_stress' in kdf.columns:
            i += 1
            axs[i].plot(kdf.kpoints_spacing, kdf.maximum_stress, '-x')
            axs[i].set_ylabel(r'$S_{max}$ (kBar)')
        axs[i].set_xticks(kdf.kpoints_spacing)
        axs[i].set_xticklabels(
            [f'{row.kpoints_spacing:.3f}\n{row.mesh}' for _, row in kdf.iterrows()],
            rotation=45,
        )
        axs[i].set_xlabel('K-pointing spacing (mesh)')
        fig.tight_layout()

    return figs


def get_convergence_builder(structure: orm.StructureData, config: dict[str, Any]) -> VaspBuilderUpdater:
    """
    Short cut for getting an VaspBuilderUpdater ready to use

    :structure StructureData: The input structure node.
    :config dict: Configuration dictionary specifying the protocol.

    The following files are used from the configuration: ``code``, ``inputset``, ``conv``, ``options``, ``resources``.
    """

    conv_builder = VaspConvergenceWorkChain.get_builder()

    upd = VaspBuilderUpdater(conv_builder)
    upd.use_inputset(
        structure,
        config.get('inputset', VaspBuilderUpdater.DEFAULT_INPUTSET),
        overrides=config.get('overrides', {}),
    )
    upd.set_code(orm.load_code(config['code']))

    upd.set_default_options(**config.get('options', {}))
    upd.update_resources(**config.get('resources', {}))
    upd.set_label(f'{structure.label} CONV')

    # Convergence specific options
    conv = ConvOptions(**config.get('conv', {}))
    upd.builder.conv_settings = conv.aiida_dict()
    return upd
