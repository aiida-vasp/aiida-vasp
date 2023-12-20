.. _vasp_workchain:

==============
VASP workchain
==============

This is the base `VASP`_ workchain, it is designed so that one can perform any single DFT `VASP`_ run. The behavior of the calculation, i.e. exactly what kind of run is performed, control values, etc. are controlled via the ``parameters`` inputs.

The :py:class:`VaspWorkChain<aiida_vasp.workchains.vasp.VaspWorkChain>` is the basis that is used to build the rest of the more specialized workchains, e.g. `RelaxWorkChain` to deal with the structural relaxation.

Required inputs
^^^^^^^^^^^^^^^

The :py:class:`VaspWorkChain<aiida_vasp.workchains.vasp.VaspWorkChain>` requires a number of inputs, these comprise the minimum set of information to run a `VASP`_ calculation from `AiiDA`_.

* ``code``, type: :py:class:`InstalledCode<aiida.orm.nodes.data.code.installed.InstalledCode>`. Describes the VASP executable and holds a reference to the :py:class:`Computer<aiida.orm.computers.Computer>` instance on which it lives.
* ``structure``, type: :py:class:`StructureData<aiida.orm.nodes.data.structure.StructureData>` or :py:class:`CifData<aiida.orm.nodes.data.cif.CifData>`. Describes the structure on which VASP is to be run.
* ``kpoints``, type: :py:class:`KpointsData<aiida.orm.nodes.data.array.kpoints.KpointsData>`. The kpoints mesh or path.
* ``potential_family``, type: :py:class:`Str<aiida.orm.nodes.data.str.Str>`. The name given to a set of uploaded POTCAR files.
* ``potential_mapping``, type: :py:class:`Dict<aiida.orm.nodes.data.dict.Dict>`. Dictionary containing an entry for at least every kind name in the ``structure`` input with the full name of the POTCAR from the ``potential_family``. Example: ``{'In1': 'In_d', 'In2': 'In_h'}``.
* ``parameters``, type: :py:class:`Dict<aiida.orm.nodes.data.dict.Dict>`. Dictionary with the parameters for the calculation. Please consult the documentation on how parameters are handled (:ref `parameters`) for details, particularly the section pertaining to the ``VaspWorkChain``.
* ``options``, type: :py:class:`Dict<aiida.orm.nodes.data.dict.Dict>`. Dictionary containing at least the keys ``resources``. More information about the options is available in the `AiiDA documentation`_.

Extra inputs
^^^^^^^^^^^^

The :py:class:`VaspWorkChain<aiida_vasp.workchains.vasp.VaspWorkChain>` can take other inputs that allow for higher control of the workchain itself.

* ``settings``, type: :py:class:`Dict<aiida.orm.nodes.data.dict.Dict>`. Dictionary containing parameters not related to `VASP`_ itself, e.g. parser settings, selective dynamics, etc.
* ``wavecar``, type: :py:class:`WavefunData<aiida_vasp.data.wavefun.WavefunData>`. It contains the wavefunctions of the Kohn-Sham equation as stored in the `WAVECAR`_ file. It can be used to restart a calculation in a very efficient manner.
* ``chgcar```, type: :py:class:`ChargedensityData<aiida_vasp.data.chargedensity.ChargedensityData>`. It contains the charge density and the PAW one-center occupancies and can be used for restarting VASP calculation, as stored in the `CHGCAR`_ file.
* ``site_magnetization``, type: :py:class:`Dict<aiida.orm.nodes.data.dict.Dict>`. Dictionary containing the site dependent magnetization, that can be used to restart the calculation. It currently it is only tested for the collinear case.
* ``restart_folder``, type: :py:class:`RemoteData<aiida.orm.nodes.data.remote.base.RemoteData>`. This is a folder of a previous calculation that can be used as a parent or to restart the calculation.
* ``max_iterations``, type: :py:class:`Int<aiida.orm.nodes.data.int.Int>`, default: 5. How many iterations the restart will be attempted before resulting in failure.
* ``clean_workdir``, type: :py:class:`Bool<aiida.orm.nodes.data.bool.Bool>`, default: True. Whether or not the remote folder of the calculation will be deleted after the end of the calculation.
* ``verbose``, type: :py:class:`Bool<aiida.orm.nodes.data.bool.Bool>`, default: False. Whether or not extra information is displayed during the workchain execution.
* ``dynamics.positions_dof``, type: :py:class:`List<aiida.orm.nodes.data.list.List>`. It controls the selective dynamics of the ions when performing relaxations.

.. _vasp_workchain_outputs:

Required outputs
^^^^^^^^^^^^^^^^

A successful :py:class:`VaspWorkChain<aiida_vasp.workchains.vasp.VaspWorkChain>` would result in the following outputs always being produced

* ``misc``, type: :py:class:`Dict<aiida.orm.nodes.data.dict.Dict>`. Dictionary containing the output parameters containing smaller quantities that do not depend on system size.

Extra outputs
^^^^^^^^^^^^^

Depending on the input variables passed as inputs one or more of the following outputs can be produced

* ``structure``, type: :py:class:`StructureData<aiida.orm.nodes.data.structure.StructureData>`. Output structure from the simulation.
* ``kpoints``, type: :py:class:`KpointsData<aiida.orm.nodes.data.array.kpoints.KpointsData>`. Output k-points mesh.
* ``trajectory``, type: :py:class:`TrajectoryData<aiida.orm.nodes.data.array.trajectory.TrajectoryData>`. Trajectory of the atomic positions.
* ``chgcar``, type: :py:class:`ChargedensityData<aiida_vasp.data.chargedensity.ChargedensityData>`. It contains the charge density and the PAW one-center occupancies and can be used for restarting VASP calculation, as stored in the `CHGCAR`_ file.
* ``wavecar``, type: :py:class:`WavefunData<aiida_vasp.data.wavefun.WavefunData>`. It contains the wavefunctions of the Kohn-Sham equation as stored in the `WAVECAR`_ file.
* ``bands``, type: :py:class:`BandsData<aiida.orm.nodes.data.array.bands.BandsData>`. The output band structure.
* ``forces``, type: :py:class:`ArrayData<aiida.orm.nodes.data.array.array.ArrayData>`. The output forces of the calculation.
* ``stress``, type: :py:class:`ArrayData<aiida.orm.nodes.data.array.array.ArrayData>`. The output stress of the calculation.
* ``dos``, type: :py:class:`ArrayData<aiida.orm.nodes.data.array.array.ArrayData>`. The output density of states of the calculation.
* ``energies``, type: :py:class:`ArrayData<aiida.orm.nodes.data.array.array.ArrayData>`. The output total energies.
* ``projectors``, type: :py:class:`ArrayData<aiida.orm.nodes.data.array.array.ArrayData>`. The output projectors of decomposition.
* ``dielectrics``, type: :py:class:`ArrayData<aiida.orm.nodes.data.array.array.ArrayData>`. The output dielectric functions.
* ``dynmat``, type: :py:class:`ArrayData<aiida.orm.nodes.data.array.array.ArrayData>`. The output dynamical matrix.
* ``charge_density``, type: :py:class:`ArrayData<aiida.orm.nodes.data.array.array.ArrayData>`. The output charge density.
* ``magnetization_density``, type: :py:class:`ArrayData<aiida.orm.nodes.data.array.array.ArrayData>`. The output magnetization density.
* ``site_magnetization``, type: :py:class:`Dict<aiida.orm.nodes.data.dict.Dict>`. Dictionary containing the site dependent magnetization.

Restarting calculations
^^^^^^^^^^^^^^^^^^^^^^^

The main difference between a :py:class:`VaspWorkChain<aiida_vasp.workchains.vasp.VaspWorkChain>` and a  :py:class:`VaspCalculation<aiida_vasp.calcs.vasp.VaspCalculation>` is that the former implements a basic logic of restarting failed or unfinished calculations.
The framework of :py:class:`~aiida.engine.processes.workchains.restart.BaseRestartWorkChain` is used with a set of predefined handlers to fix some (but not all) common pitfalls,
such as restarting an ionic relaxation that has run out of the wall time and electronic convergence issues.

Once a calculation is finished, the :py:class:`CalculationNode<aiida.orm.nodes.process.calculation.calculation.CalculationNode>` is inspected by a series of :py:func:`process_handler<aiida.engine.processes.workchains.utils.process_handler>`,
which are executed in the order of descending priority.
Each handler may be tied to a specific list of ``exit_code`` that the calculation may have.
If any problems are found, and the restart can be performed, a :py:class:`ProcessHandlerReport<aiida.engine.processes.workchains.utils.ProcessHandlerReport>` would be returned and added to a list.
If the ``break`` attribute of the report is set to ``True`` the handling process would be terminated.
Afterwards, the last report is inspected. If it has an none-zero ``exit_code`` the, then the workchain will be aborted with that ``exit_code`` returned, this corresponds to the case where the error cannot be corrected automatically.
Otherwise, it is assumed that calculation should be restarted with the revised inputs.

The flow chart below illustrates how it works coupled with the emission of the :py:class:`ProcessHandlerReport<aiida.engine.processes.workchains.utils.ProcessHandlerReport>` from the handlers:

.. image:: process-handler.png

For more information, please see the docstring of :py:class:`~aiida.engine.processes.workchains.restart.BaseRestartWorkChain`.

One should note that the handlers included here are not intended to give a comprehensive coverage of all of possible errors from VASP,
but instead we focus on improving the robustness by performing simple corrections that would be the right things to do in most times.

New handlers may be registered by adding the method to :py:class:`VaspWorkChain<aiida_vasp.workchains.vasp.VaspWorkChain>` with the :py:func:`process_handler<aiida.engine.processes.workchains.utils.process_handler>` decorator.
Alternatively, one can also extended the :py:class:`VaspWorkChain<aiida_vasp.workchains.vasp.VaspWorkChain>` by sub-classing and add more handlers there.


.. _AiiDA: https://www.aiida.net
.. _VASP: https://www.vasp.at
.. _AiiDA documentation: http://aiida-core.readthedocs.io/en/latest/
.. _Workchain: https://aiida.readthedocs.io/projects/aiida-core/en/latest/concepts/workflows.html#work-chains
.. _WAVECAR: https://www.vasp.at/wiki/index.php/WAVECAR
.. _CHGCAR: https://www.vasp.at/wiki/index.php/CHGCAR
