(simple-calculation)=
# Running a simple VASP calculation

In this example, we will run a simple VASP calculation using the aiida-vasp plugin.
The `VaspWorkChain` will be used to run a calculation on a sample structure.

It is assumed that the [installation][#installation] has been completed, the `InstalledCode` for the VASP executable and the `PBE.54` pseduopotentials has been [configured][#configuration].

The following steps can be run with a ipython shell launched by the `verdi shell` command.

## Setting up the inputs
First, we create the sample structure using the `ase` package (assuming the `ase` package is installed):
```python
from ase.build import bulk
from aiida import orm
si = bulk('Si', 'diamond', 5.4)
si_structure = orm.StructureData(ase=si)
```

Note that the `orm.StructureData` class is used to store the structure data in the AiiDA database.

Next, we create necessary input node and set the required parameters:

```python
from aiida.plugins import WorkflowFactory
builder = WorkflowFactory('vasp.v2.vasp').get_builder()

builder.parameters = orm.Dict(
    dict={'incar':
    {
        'encut': 300,
        'ismear': 0,
        'ispin': 1,
    }
    }
)
builder.options = orm.Dict(dict={
    'resources': {'num_machines': 1,
    'tot_num_mpiprocs': 1}}
    )
builder.structure = si_structure
builder.code = orm.load_code('vasp@localhost')   # Assuming VASP is installed locally and the code is configured
builder.potential_family = 'PBE.54'     # Note that we can use python type - they are converted into AiiDA types automatically
builder.potential_mapping = {'Si': 'Si'}
builder.kpoints_spacing = 0.05
builder.metadata.label = 'test-calculation'
builder
```

It is a good practice to inspect the `ProcessBuilder` before submitting the calculation.

## Run and inspect the results

Finally, the calculation can be run with the `run_get_node` method:

```python
from aiida.engine import run_get_node
out_dict, process_node = run_get_node(builder)
```

If everything runs fine, you should see something like

```
Report: [1615|VaspWorkChain|run_process]: launching VaspCalculation<1618> iteration #1
Report: [1615|VaspWorkChain|results]: work chain completed after 1 iterations
Report: [1615|VaspWorkChain|on_terminated]: cleaned remote folders of calculations: 1618
```

The `out_dict` is dictionary containing the output nodes of the `VaspWorkChain` and the `process_node` is a `WorkChainNode` representing the work chain that has been executed.

```python
>>> out_dict
{'remote_folder': <RemoteData: uuid: c9351bfb-237c-4a80-9710-8f3136f44e4a (pk: 1619)>,
 'retrieved': <FolderData: uuid: 8b84eeb2-b58f-474d-9d33-eb8a0f48d61d (pk: 1620)>,
 'misc': <Dict: uuid: a5c71be6-0a6d-4ef3-b2f2-ba91fbc27a67 (pk: 1621)>,
 'dielectrics': <ArrayData: uuid: e24207ea-cf2e-4459-b2b4-66feef172aa9 (pk: 1622)>}
>>> process_node
<WorkChainNode: uuid: a84d5dfb-03be-4b49-86d0-d81b3b9f31e4 (pk: 1615) (aiida.workflows:vasp.v2.vasp)>
```

:::{note}
The actual `uuid` and `pk` displayed on your compute will be different.
:::

The total energy, forces and stress can be found in the `misc` output node:

```python
>>> out_dict['misc']['total_energies']
{'energy_extrapolated': -10.81676131,
 'energy_extrapolated_electronic': -10.81676131}
>>> out_dict['misc']['forces']
{'final': [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]]}
>>> out_dict['misc']['stress']
{'final': [[35.63282046, 0.0, 0.0],
  [0.0, 35.63282046, 0.0],
  [0.0, 0.0, 35.63282046]]}
```

The forces are zero as the structure has high symmetry.
The stresses are not zero since we did not put the equilibrium lattice constant of Si in the first place.
Geometry optimization is needed to fully relax the structure.

We also can access the same output through the `WorkChainNode` using the `outputs` attribute:

```python
>>> process_node.outputs.misc['run_stats']
{'nsw': 0,
 'nelm': 60,
 'nbands': 8,
 'finished': True,
 'ionic_converged': None,
 'contains_nelm_breach': False,
 'electronic_converged': True,
 'last_iteration_index': [1, 10],
 'consistent_nelm_breach': False}
```
