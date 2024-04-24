"""
Redesigned convergence testing workchain.

Make it more simple and easy to use.

Inputs
- just like a normal VaspWorkChain
- run with different cut off energies
- run with different k spacing
- summarise the results

No added wrapper etc.
"""

from aiida import orm
from aiida.common.utils import classproperty
from aiida.engine import WorkChain, append_, calcfunction
from aiida.plugins import WorkflowFactory

from .common import nested_update_dict_node
from .common.opthold import FloatOption, OptionContainer

#pylint:disable=no-member,unused-argument,no-self-argument,import-outside-toplevel


class VaspConvergenceWorkChain(WorkChain):
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

    @classmethod
    def define(cls, spec):
        super().define(spec)

        spec.expose_inputs(cls._sub_workchain)
        spec.input('conv_settings', help='Settings of the workchain', valid_type=orm.Dict)
        spec.outline(cls.setup, cls.launch_conv_calcs, cls.analyse)

        spec.exit_code(
            401,
            'ERROR_SUBWORKFLOW_ERRORED',
            message='At leaste one of the launched sub-workchain has failed',
        )
        spec.output('kpoints_conv_data', required=False)
        spec.output('cutoff_conv_data', required=False)

    def setup(self):
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

    def launch_conv_calcs(self):
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
            new_param = nested_update_dict_node(inputs.parameters, {'incar': {'encut': cut}})
            inputs.parameters = new_param
            if original_label:
                inputs.metadata.label = original_label + f" CUTCONV {cut:.2f}"
            else:
                inputs.metadata.label = f"CUTCONV {cut:.2f}"

            running = self.submit(self._sub_workchain, **inputs)
            self.report(f"Submitted {running} with cut off energy {cut:.1f} eV.")
            self.to_context(cutoff_conv_workchains=append_(running))

        # Launch kpoints convergence tests
        new_param = nested_update_dict_node(inputs.parameters, {'incar': {'encut': cutoff_for_kconv}})
        for kspacing in self.ctx.kspacing_list:
            inputs.parameters = new_param
            inputs.kpoints_spacing = kspacing
            if original_label:
                inputs.metadata.label = original_label + f" KCONV {kspacing:.3f}"
            else:
                inputs.metadata.label = f"KCONV {kspacing:.3f}"

            running = self.submit(self._sub_workchain, **inputs)
            self.report(f"Submitted {running} with kpoints spacing {kspacing:.3f}.")
            self.to_context(kpoints_conv_workchains=append_(running))

    def analyse(self):
        """
        Analyse the output of the calculations.
        Collect data to be plotted/analysed against the cut off energy and kpoints spacing
        """

        def collect_data(workchain):
            """Collect the data from workchain output"""
            output = workchain.outputs.misc.get_dict()
            data = {}
            data['maximum_force'] = output.get('maximum_force')
            # Extract the magnetization
            magnetization = output.get('magnetization')
            if magnetization:
                data['magnetization'] = magnetization[0]
            data['maximum_stress'] = output.get('maximum_stress', None)
            data['energy'] = list(output['total_energies'].values())[0]
            return data

        def unpack(name, input_data):
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
        if 'cutoff_conv_workchains' in self.ctx:
            for iwork, workchain in enumerate(self.ctx.cutoff_conv_workchains):
                if workchain.exit_status != 0:
                    exit_code = self.exit_codes.ERROR_SUBWORKFLOW_ERRORED
                    self.report(f"Skipping workchain {workchain} with exit status {workchain.exit_status} ")
                    continue

                cutoff = workchain.inputs.parameters['incar']['encut']
                cutoff_data[cutoff] = collect_data(workchain)
                cutoff_data[cutoff]['mesh'] = workchain.called[0].inputs.kpoints.get_kpoints_mesh()[0]
                cutoff_miscs[f"worchain_{iwork}"] = workchain.outputs.misc

        kspacing_data = {}
        kspacing_miscs = {}

        if 'kpoints_conv_workchains' in self.ctx:
            for iwork, workchain in enumerate(self.ctx.kpoints_conv_workchains):
                if workchain.exit_status != 0:
                    exit_code = self.exit_codes.ERROR_SUBWORKFLOW_ERRORED
                    self.report(f"Skipping Workchain {workchain} with exit status {workchain.exit_status} ")
                    continue

                spacing = float(workchain.inputs.kpoints_spacing)
                kspacing_data[spacing] = collect_data(workchain)
                kspacing_data[spacing]['mesh'] = workchain.called[0].inputs.kpoints.get_kpoints_mesh()[0]
                kspacing_miscs[f"worchain_{iwork}"] = workchain.outputs.misc
                cutoff = workchain.inputs.parameters['incar']['encut']
                kspacing_data[spacing]['cutoff_energy'] = cutoff

        # Calcfunction to link with the calculation output to the summary data node
        @calcfunction
        def create_links_kconv(**miscs):
            """Alias calcfunction to link summary node with miscs"""
            return orm.Dict(dict=unpack('kpoints_spacing', kspacing_data))

        @calcfunction
        def create_links_cutconv(**miscs):
            """Alias calcfunction to link summary node with miscs"""
            return orm.Dict(dict=unpack('cutoff_energy', cutoff_data))

        if kspacing_data:
            self.out('kpoints_conv_data', create_links_kconv(**kspacing_miscs))
        if cutoff_data:
            self.out('cutoff_conv_data', create_links_cutconv(**cutoff_miscs))

        return exit_code

    @classproperty
    def option_class(cls):
        return ConvOptions

    @staticmethod
    def get_conv_data(conv_work, plot=False, **plot_kwargs):
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


class ConvOptions(OptionContainer):
    """Template for the Dict node controlling the workchain behaviour"""

    cutoff_start = FloatOption('The starting cut-off energy', 300)
    cutoff_stop = FloatOption('The Final cut-off energy', 700)
    cutoff_step = FloatOption('Step size of the cut-off energy', 50)
    kspacing_start = FloatOption('The starting kspacing', 0.07)
    kspacing_stop = FloatOption('The final kspacing', 0.02)
    kspacing_step = FloatOption('Step size of the cut-off energy', 0.01)
    cutoff_kconv = FloatOption('The cut-off energy used for kpoints convergence tests', 450)
    kspacing_cutconv = FloatOption('The kpoints spacing used for cut-off energy convergence tests', 0.07)


def get_conv_data(conv_work):
    """
    Convenient method for extracting convergence data

    Args:
      conv_work (orm.WorkChainNode): Convergence workflow node

    Returns:
      A tuple of cut-off convergence and k-point convergence result data frame
    """
    import pandas as pd

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


def plot_conv_data(cdf, kdf, **kwargs):
    """
    Make two combined plots for the convergence test results.
    """
    import matplotlib.pyplot as plt

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
            [f"{row.kpoints_spacing:.3f}\n{row.mesh}" for _, row in kdf.iterrows()],
            rotation=45,
        )
        axs[i].set_xlabel('K-pointing spacing (mesh)')
        fig.tight_layout()

    return figs


def get_convergence_builder(structure, config):
    """
    Short cut for getting an VaspBuilderUpdater ready to use

    :structure StructureData: The input structure node.
    :config dict: Configuration dictionary specifying the protocol.

    The following files are used from the configuration: ``code``, ``inputset``, ``conv``, ``options``, ``resources``.
    """

    from .common.builder_updater import VaspBuilderUpdater

    conv_builder = VaspConvergenceWorkChain.get_builder()

    upd = VaspBuilderUpdater(conv_builder)
    upd.use_inputset(
        structure,
        config.get('inputset', VaspBuilderUpdater.DEFAULT_SET),
        overrides=config.get('overrides', {}),
    )
    upd.set_code(orm.load_code(config['code']))

    upd.set_default_options(**config.get('options', {}))
    upd.update_resources(**config.get('resources', {}))
    upd.set_label(f"{structure.label} CONV")

    # Convergence specific options
    conv = ConvOptions(**config.get('conv', {}))
    upd.builder.conv_settings = conv.to_aiida_dict()
    return upd
