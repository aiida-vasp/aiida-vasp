How to relax a structure using vasp.relax
=========================================

Running a single relaxation
----------------------------

An example containing necessary steps. Can be found in `example base_wf file`_. This includes creating sample input data tand connecting it to the inputs of the workchain and finally how to submit the workchain The only steps not covered are the setup of AiiDA, which is documented in the `AiiDA docs`_

The example relaxes the ion positions only. Other degrees of freedom are available, see the detailed `input reference <howto/relax_wf/inputs>`

.. _example base_wf file: https://github.com/aiidateam/aiida-vasp/blob/develop/examples/run_relax_wf.py
.. _AiiDA docs: https://aiida-core.readthedocs.io/en/stable/work/index.html

Writing workflows that include a relaxation step
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Unless the workflow is solely concerned with providing a more robust or more specialized relaxation step, it is not recommended to derive from `vasp.relax`. Instead composing should be used; your workchain should set up the inputs and submit the relaxation workchain in a step function which will commit the results of the relaxation workchain to it's context. Further steps can then use the `relaxed_structure` result.

.. _howto/relax_wf/inputs:
Reference: vasp.relax inputs
----------------------------

Degrees of freedom
^^^^^^^^^^^^^^^^^^

Through its input parameters, VASP can be configured to utilize three degrees of freedoms for relaxations: ion positions, cell volume and cell shape. Some, but not all combinations are allowed, read more about that in the `official VASP wiki - ISIF tag page`_. Other posibilities are also posible, but not covered here and typically demands a dedicated VASP version with changes in to the standard source code embedded.

`vasp.relax` allows to switch each degree of freedom on / off independently, setting the ISIF and IBRION tags accordingly. Each of these inputs is optional and by default only the ion positions are relaxed.

 * `relax.positions`, type: `Bool`, default: True (on)
 * `relax.shape`, type: `Bool`, default: False (off)
 * `relax.volume`, type: `Bool`, default: False (off)

.. _official VASP wiki - ISIF tag page: https://cms.mpi.univie.ac.at/wiki/index.php/ISIF

k-points
^^^^^^^^

There are two ways to set the k-point grid on which to relax the structure: by passing a `KpointsData` node with a mesh, or to pass a maximum distance from k-point to k-point, from which a gamma centered mesh will be constructed. One and only one of the two must be passed.

 * `kpoints.mesh`, type: `KpointsData`
 * `kpoints.distance`, type : `Float`, units: `1/Angstrom`

Local minima
^^^^^^^^^^^^^^^^

The `vasp.relax` workchain is able to check for convergence on any of the available degrees of freedom by running a fresh relaxation from the output structure of the previous run. This can sometimes lead to further relaxation, if the previous run got stuck in a local charge density minimum. The new calculation starts from scratch with a randomized charge density, but with the last obtained positions. This is done iteratively until the target property does not change more than a given tolerance. Currently, external check on force, stress and energy is not implemented, but will be available as an option in the near future.

Keep in mind there is no guarantee that the new run will overcome the barriers of a local minimum. More in-depth workchains could be developed to do that, by deriving from this workchain or using it as a building block. This feature is switched off by default.

All of these inputs are optional

 * `convergence.on`, type: `Bool`, default: False (do not check or run additional relaxations)
 * `convergence.absolute`, type: `Bool`, default: False (with respect to the previous relaxation)
   - False: relative tolerances are used (relative convergence)
   - True: absolute tolerances are used (native VASP units)
 * `convergence.max_iterations`, type: `Int`, default: 5 (execute a maximum of five relaxation runs)
 * `convergence.shape.lengths`, type: `Float`, default: 0.1 (allow a maximum of 10 % change of the L2 norm for the unitcell vectors from the previous relaxation)
 * `convergence.shape.angles`, type: `Float`, default: 0.1 (allow a maximum of 10 % change of the unitcell angles from the previous relaxation)
 * `convergence.volume`, type: `Float`, default: 0.01 (allow a maximum of 1 % change of the unitcell volume from the previous relaxation)
 * `convergence.positions`, type: `Float`, default: 0.01 (allow a maximum of 1 % displacement (L2 norm) of the positions from the previous relaxation)

 Exposed from `vasp.base`
 ^^^^^^^^^^^^^^^^^^^^^^^^

 The following inputs will be passed straight through to a `vasp.base` workchain. See the `reference there <howto/base_wf/reference>`_ for how to use them.

 * `code`
 * `structure`
 * `potcar_family`
 * `potcar_mapping`
 * `options`
 * `restart.max_iterations` -> `max_iterations` on `vasp.base`
 * `restart.clean_workdir` -> `clean_workdir` on `vasp.base`

Reference: vasp.relax outputs
-----------------------------

The following output nodes are created upon successful completion:

 * `output_parameters`, type: `ParameterData`, the default scalar / fixed dimension output properties of the final relaxation calculation. These might be used to do further sanity checks or take decisions for later stages of a more complex workchain.
 * `relaxed_structrue`, type: `StructureData`, The relaxed structure.
