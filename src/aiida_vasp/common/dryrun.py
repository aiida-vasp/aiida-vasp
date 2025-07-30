"""
Module to provide dryrun functionality.
"""

from __future__ import annotations

import shutil
from math import ceil, gcd
from typing import Optional
from warnings import warn

import numpy as np
from aiida.common.folders import SubmitTestFolder
from aiida.engine.daemon.execmanager import upload_calculation
from aiida.engine.processes.builder import ProcessBuilder
from aiida.orm import Dict, KpointsData
from aiida.transports.plugins.local import LocalTransport

from aiida_vasp.assistant.parameters import ParametersMassage
from aiida_vasp.calcs.vasp import VaspCalculation
from aiida_vasp.commands.dryrun_vasp import dryrun_vasp as _vasp_dryrun
from aiida_vasp.data.potcar import PotcarData


class JobScheme:
    """
    A class representing the scheme of the jobs.
    """

    def __init__(
        self,
        n_kpoints: int,
        n_procs: int,
        n_nodes: Optional[int] = None,
        cpus_per_node: Optional[int] = None,
        npw: Optional[int] = None,
        nbands: Optional[int] = None,
        ncore_within_node: bool = True,
        ncore_strategy: str = 'maximise',
        wf_size_limit: float = 1000,
    ) -> None:
        """
        Instantiate a JobScheme object.

        :param n_kpoints: Number of kpoints.
        :param n_procs: Number of processes.
        :param n_nodes: Number of nodes.
        :param cpus_per_node: Number of CPUs per node.
        :param npw: Number of plane waves.
        :param nbands: Number of bands.
        :param ncore_within_node: If True, limit plane-wave parallelisation to within each node.
        :param ncore_strategy: Strategy for optimising NCORE, choose from 'maximise' and 'balance'.
        :param wf_size_limit: Limit of the ideal wavefunction size per process in MB.
        """
        self.n_kpoints = n_kpoints
        self.n_procs = n_procs
        self.n_nodes = n_nodes
        self.cpus_per_node = cpus_per_node
        self.npw = npw
        self.nbands = nbands
        self.ncore_within_node = ncore_within_node
        self.ncore_strategy = ncore_strategy
        self.wf_size_limit = wf_size_limit

        self.n_kgroup = None  # KPOINT groups
        self.n_bgroup = None  # Band groups
        self.n_pgroup = None  # Plane wave groups

        self.kpar = None  # Value for the KPAR
        self.npar = None
        self.ncore = None  # Value for the ncore
        self.new_nbands = nbands  # Value for the new nbands
        self.nbands_amplification = None  # Amplification factor for the NBAND round up
        self.ncore_balance = None  # NCORE/NPAR balance factor

        self.solve_kpar()
        self.solve_ncore()

    @classmethod
    def from_dryrun(cls, dryrun_outcome: dict, n_procs: int, **kwargs) -> JobScheme:
        """
        Construct from dryrun results.

        :param dryrun_outcome: The outcome from a dryrun.
        :param n_procs: Number of processes.
        :param kwargs: Additional keyword arguments.

        :returns: A `JobScheme` object
        """
        kwargs['n_kpoints'] = dryrun_outcome.get('num_kpoints')
        kwargs['nbands'] = dryrun_outcome.get('num_bands')
        kwargs['npw'] = dryrun_outcome.get('num_plane_waves')
        kwargs['n_procs'] = n_procs
        return cls(**kwargs)

    def solve_kpar(self) -> int:
        """
        Solve for the optimum strategy for KPAR.

        :returns: The optimized KPAR value.
        """
        kpar = gcd(self.n_kpoints, self.n_procs)
        self.kpar = kpar
        # If we did not set nbands or npw, we cannot adjust KAR to avoid memory issues
        if any(map(lambda x: x is None, [self.nbands, self.npw])):
            warn(
                'Cannot limit KAR for memory requirement without supplying both NBANDS and NPW',
                UserWarning,
            )
            return kpar

        # Reduce the KPAR
        if self.size_wavefunction_per_proc > self.wf_size_limit:
            for candidate in factors(kpar):
                self.kpar = candidate
                if self.size_wavefunction_per_proc < self.wf_size_limit:
                    kpar = candidate
                    break
        if self.size_wavefunction_per_proc > self.wf_size_limit:
            warn(
                ('Expected wavefunction size per process {} MB is larger than the limit {} MB').format(
                    self.size_wavefunction_per_proc, self.wf_size_limit
                ),
                UserWarning,
            )
        return kpar

    @property
    def nk_per_group(self) -> int:
        """Number of kpoints per group."""
        return self.n_kpoints // self.kpar

    @property
    def procs_per_kgroup(self) -> int:
        """Number of processes per kpoint group."""
        return self.n_procs // self.kpar

    def solve_ncore(self) -> int:
        """
        Solve for NCORE.

        :returns: The optimized NCORE value.
        """
        # Cannot solve if no nbands provided or does not know how many cpus per node
        if self.nbands is None:
            return
        if self.ncore_within_node and (self.cpus_per_node is None):
            return

        combs = []
        for ncore in factors(self.procs_per_kgroup):
            if ncore > 12:
                continue
            # Only consider ncore that is a multiple of the cpus per node
            if self.ncore_within_node and self.cpus_per_node % ncore != 0:
                continue
            npar = self.procs_per_kgroup // ncore
            new_nbands = ceil(self.nbands / npar) * npar
            factor = new_nbands / self.nbands
            combs.append((ncore, factor, abs(ncore / npar - 1), new_nbands))

        combs = list(filter(lambda x: x[1] < 1.2, combs))
        if self.ncore_strategy == 'balance':
            combs.sort(key=lambda x: x[2])
        elif self.ncore_strategy == 'maximise':
            combs.sort(key=lambda x: x[0], reverse=True)
        else:
            raise RuntimeError(f'NCORE strategy: <{self.ncore_strategy}> is invalid')

        # Take the maximum ncore
        ncore, factor, balance, new_nbands = combs[0]

        self.ncore = ncore
        self.npar = self.procs_per_kgroup // ncore
        self.nbands_amplification = factor
        self.new_nbands = new_nbands
        self.ncore_balance = balance
        return ncore

    @property
    def size_wavefunction(self) -> float:
        """Memory requirement for the wavefunction in MB."""
        return self.n_kpoints * self.new_nbands * self.npw * 16 / 1048576

    @property
    def size_wavefunction_per_proc(self) -> float:
        """Memory requirement for the wavefunction per process."""
        return self.size_wavefunction / self.procs_per_kgroup


def factors(num: int) -> list[int]:
    """
    Return all factors of a number in descending order, including the number itself.

    :param num: The number to factor.

    :returns: A list of factors.
    """
    result = [num]
    for i in range(num // 2 + 1, 0, -1):
        if num % i == 0:
            result.append(i)
    return result


def dryrun_vasp(
    input_dict: dict | ProcessBuilder,
    vasp_exe: str = 'vasp_std',
    timeout: int = 10,
    work_dir: str | None = None,
    keep: bool = False,
) -> dict:
    """
    Perform a dryrun for a VASP calculation, return obtained information.

    :param input_dict: The input dictionary/builder for `VaspCalculation`.
    :param vasp_exe: The VASP executable to be used.
    :param timeout: Timeout for the underlying VASP process in seconds.
    :param work_dir: Working directory, if not supplied, will use a temporary directory.
    :param keep: Whether to keep the dryrun output.

    :returns: A dictionary of the dry run results parsed from OUTCAR.
    """
    # Deal with passing an process builder
    if isinstance(input_dict, ProcessBuilder):
        try:
            output_node = prepare_inputs(input_dict)
        except Exception as error:
            raise error
    else:
        try:
            output_node = prepare_inputs(input_dict)
        except Exception as error:
            raise error

    folder = output_node.dry_run_info['folder']
    outcome = _vasp_dryrun(folder, vasp_exe=vasp_exe, timeout=timeout, work_dir=work_dir, keep=keep)
    if not keep:
        shutil.rmtree(folder)

    return outcome


def get_jobscheme(input_dict: dict, nprocs: int, vasp_exe: str = 'vasp_std', **kwargs) -> JobScheme:
    """
    Perform a dryrun for the input and work out the best parallelisation strategy.

    :param input_dict: Inputs of the VaspCalculation.
    :param nprocs: Target number of processes to be used.
    :param vasp_exe: The executable of local VASP program to be used.
    :param kwargs: Additional keyword arguments to be passed to `JobScheme`.

    :returns: A `JobScheme` object.
    """
    dryout = dryrun_vasp(input_dict, vasp_exe)
    scheme = JobScheme.from_dryrun(dryout, nprocs, **kwargs)
    return scheme


def prepare_inputs(inputs: dict) -> VaspCalculation:
    """
    Prepare inputs for VASP calculation.

    :param inputs: The inputs to prepare.

    :returns: The prepared inputs.
    """
    inputs = dict(inputs)
    inputs['metadata'] = dict(inputs['metadata'])
    inputs['metadata']['store_provenance'] = False
    inputs['metadata']['dry_run'] = True
    vasp = VaspCalculation(inputs=inputs)

    with LocalTransport() as transport:
        with SubmitTestFolder() as folder:
            calc_info = vasp.presubmit(folder)
            transport.chdir(folder.abspath)
            upload_calculation(
                vasp.node,
                transport,
                calc_info,
                folder,
                inputs=vasp.inputs,
                dry_run=True,
            )
            vasp.node.dry_run_info = {
                'folder': folder.abspath,
                'script_filename': vasp.node.get_option('submit_script_filename'),
            }
    return vasp.node


def dryrun_relax_builder(builder: ProcessBuilder, **kwargs) -> dict:
    """
    Dry run a relaxation workchain builder.

    :param builder: The builder to dry run.
    :param kwargs: Additional keyword arguments.

    :returns: The results of the dry run.
    """
    vasp_builder = VaspCalculation.get_builder()

    # Setup the builder for the bare calculation

    # Convert into CalcJob inputs
    vasp_builder.code = builder.vasp.code
    pdict = builder.vasp.parameters.get_dict()
    parameters_massager = ParametersMassage(pdict, None)

    vasp_builder.parameters = Dict(dict=parameters_massager.parameters.incar)
    if 'dynamics' in parameters_massager.parameters:
        vasp_builder.dynamics = Dict(dict=parameters_massager.parameters.dynamics)

    if builder.vasp.kpoints is not None:
        vasp_builder.kpoints = builder.vasp.kpoints
    else:
        vasp_builder.kpoints = KpointsData()
        vasp_builder.kpoints.set_cell_from_structure(builder.structure)
        vasp_builder.kpoints.set_kpoints_mesh_from_density(builder.vasp.kpoints_spacing.value * np.pi * 2)
    vasp_builder.metadata.options = builder.vasp.options.get_dict()  # pylint: disable=no-member
    vasp_builder.potential = PotcarData.get_potcars_from_structure(
        builder.structure,
        builder.vasp.potential_family.value,
        builder.vasp.potential_mapping.get_dict(),
    )
    vasp_builder.structure = builder.structure

    return dryrun_vasp(vasp_builder, **kwargs)


def dryrun_vaspu_builder(builder: ProcessBuilder, **kwargs) -> dict:
    """
    Dry run a vaspu.vasp workchain builder.

    :param builder: The builder to dry run.
    :param kwargs: Additional keyword arguments.

    :returns: The results of the dry run.
    """
    pdict = builder.parameters.get_dict()
    vasp_builder = VaspCalculation.get_builder()

    parameters_massager = ParametersMassage(pdict, None)

    vasp_builder.parameters = Dict(dict=parameters_massager.parameters.incar)
    if 'dynamics' in parameters_massager.parameters:
        vasp_builder.dynamics = Dict(dict=parameters_massager.parameters.dynamics)

    # Setup the builder for the bare calculation
    vasp_builder.code = builder.code
    if builder.kpoints is not None:
        vasp_builder.kpoints = builder.kpoints
    else:
        vasp_builder.kpoints = KpointsData()
        vasp_builder.kpoints.set_cell_from_structure(builder.structure)
        vasp_builder.kpoints.set_kpoints_mesh_from_density(builder.kpoints_spacing.value * np.pi * 2)
    vasp_builder.metadata.options = builder.options.get_dict()  # pylint: disable=no-member
    vasp_builder.potential = PotcarData.get_potcars_from_structure(
        builder.structure,
        builder.potential_family.value,
        builder.potential_mapping.get_dict(),
    )
    vasp_builder.structure = builder.structure

    return dryrun_vasp(vasp_builder, **kwargs)
