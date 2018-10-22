How to relax a structure using vasp.relax
=========================================

Running a single relaxation
----------------------------

An example containing necessary steps. Can be found in `example run_relax file`_. This includes creating sample input data and connecting it to the inputs of the workchain and finally how to submit the workchain The only steps not covered are the setup of AiiDA, which is documented in the `AiiDA docs`_

The example relaxes the ion positions only. Other degrees of freedom are available, see the detailed `input reference <howto/relax_wc/inputs>`

.. _example run_relax file: https://github.com/aiidateam/aiida-vasp/blob/develop/examples/run_relax.py
.. _AiiDA docs: https://aiida-core.readthedocs.io/en/stable/work/index.html

Writing workflows that include a relaxation step
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Unless the workflow is solely concerned with providing a more robust or more specialized relaxation step, it is not recommended to derive from `vasp.relax`. Instead composing should be used; your workchain should set up the inputs and submit the relaxation workchain in a step function which will commit the results of the relaxation workchain to it's context. Further steps can then use the `output_structure_relaxed` result.

.. _howto/relax_wc/inputs:
Reference: vasp.relax inputs
----------------------------

Degrees of freedom
^^^^^^^^^^^^^^^^^^

Through its input parameters, VASP can be configured to utilize three degrees of freedom for relaxations: ion positions, cell volume and cell shape. Some, but not all combinations are allowed, read more about that in the `official VASP wiki - ISIF tag page`_. Other posibilities are also posible, but not covered here and typically demands a dedicated VASP version with hard coded changes to the source code.

`vasp.relax` allows to switch each degree of freedom on / off independently, setting the ISIF and IBRION tags accordingly. Each of these inputs is optional and by default only the ion positions are relaxed.

 * `relax.positions`, type: `Bool`, default: True (on)
 * `relax.shape`, type: `Bool`, default: False (off)
 * `relax.volume`, type: `Bool`, default: False (off)

.. _official VASP wiki - ISIF tag page: https://cms.mpi.univie.ac.at/wiki/index.php/ISIF

k-points
^^^^^^^^

A valid KpointData has to be supplied.

Local minima
^^^^^^^^^^^^^^^^

The `vasp.relax` workchain is able to check for convergence on any of the available degrees of freedom by running a fresh relaxation from the output structure of the previous run. This can sometimes lead to further relaxation, if the previous run got stuck in a local charge density minimum. The new calculation starts from scratch with a randomized charge density, but with the last obtained positions. This is done iteratively until the target property does not change more than a given tolerance. Currently, external check on force, stress and energy is not implemented, but will be available as an option in the future.

Keep in mind there is no guarantee that the new run will overcome the barriers of a local minimum. More in-depth workchains could be developed to do that, by deriving from this workchain or using it as a building block. This feature is switched off by default.

All of these inputs are optional

 * `relax`, type: `Bool`, default: False (do not perform relaxations)
 * `steps`, type: `Int`, default: 60 (the number of ionic positions updates to perform)
 * `positions`, type: `Bool`, default: True (if relax is True, perform relaxations of the atomic positions)
 * `shape`, type: `Bool`, default: False (if relax is True, perform relaxation of the cell shape)
 * `volume`, type: `Bool`, default: False (if relax is True, perform relaxation of the cell volume)
 * `convergence_on`, type: `Bool`, default: False (do not check or run additional relaxations)
 * `convergence_absolute`, type: `Bool`, default: False (with respect to the previous relaxation)
   - False: relative tolerances are used (relative convergence)
   - True: absolute tolerances are used (native VASP units)
 * `convergence_max_iterations`, type: `Int`, default: 5 (execute a maximum of five relaxation runs)
 * `convergence_shape_lengths`, type: `Float`, default: 0.1 (allow a maximum of 10 % change of the L2 norm for the unitcell vectors from the previous relaxation)
 * `convergence_shape_angles`, type: `Float`, default: 0.1 (allow a maximum of 10 % change of the unitcell angles from the previous relaxation)
 * `convergence_volume`, type: `Float`, default: 0.01 (allow a maximum of 1 % change of the unitcell volume from the previous relaxation)
 * `convergence_positions`, type: `Float`, default: 0.01 (allow a maximum of 1 % displacement (L2 norm) of the positions from the previous relaxation)

 Exposed from `vasp.vasp`
 ^^^^^^^^^^^^^^^^^^^^^^^^

 The following inputs will be passed straight through to a `vasp.vasp` workchain. See the `reference there <howto/base_wf/reference>`_ for how to use them.

 * `code`
 * `potential_family`
 * `potential_mapping`
 * `options`
 * `max_iterations` -> `max_iterations` on `vasp.vasp`
 * `clean_workdir` -> `clean_workdir` on `vasp.vasp`

Reference: vasp.relax outputs
-----------------------------

The following output nodes are created upon successful completion:

 * `output_parameters`, type: `ParameterData`, the default scalar / fixed dimension output properties of the final relaxation calculation. These might be used to do further sanity checks or take decisions for later stages of a more complex workchain.
 * `output_structure_relaxed`, type: `StructureData`, The relaxed structure.
