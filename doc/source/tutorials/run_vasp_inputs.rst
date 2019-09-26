.. _tutorial2:

=================================
2. Test launch a VASP calculation 
=================================

In this part we will simply launch a VASP calculation of silicon using the standard PBE silicon potential that you now should have sitting in your database.

1. Please fetch the ``run_vasp_lean.py`` example::

     wget https://github.com/aiida-vasp/aiida-vasp/raw/develop/examples/run_vasp_lean.py

2. First, make sure your AiiDA daemon is running::

     (aiida) max@workhorse:~$ verdi daemon status
     Profile: generic
     Daemon is running as PID 32330 since 2019-09-25 18:34:15
     Active workers [1]:
     PID    MEM %    CPU %  started
     -----  -------  -------  -------------------
     32334    1.068        0  2019-09-25 18:34:15
     Use verdi daemon [incr | decr] [num] to increase / decrease the amount of workers

3. Have a look at the ``run_vasp_lean.py`` file. We will go trough this together.

4. Submit the VASP calculation by executing::

     python run_vasp_lean.py

5. Check the progress::

     (aiida) max@workhorse:~$ verdi process list -a
     PK  Created    Process label    Process State     Process status
     ----  ---------  ---------------  ----------------  ----------------
     34  6h ago     VaspWorkChain    ⏹ Finished [0]
     35  6h ago     VaspCalculation  ⏹ Finished [0]
     
     Total results: 5
     
     Info: last time an entry changed state: 6h ago (at 15:00:40 on 2019-09-25)
     

6. When the calculation is finished, please have a look at the report of the workchain::

     (aiida) max@workhorse:~$ verdi process report 34
     2019-09-25 14:59:34 [5  | REPORT]: [34|VaspWorkChain|run_calculation]: launching VaspCalculation<35> iteration #1
     2019-09-25 15:00:36 [6  | REPORT]: [34|VaspWorkChain|_handle_succesfull]: VaspCalculation<35> completed successfully
     2019-09-25 15:00:36 [7  | REPORT]: [34|VaspWorkChain|results]: VaspWorkChain<34> completed after 1 iterations
     2019-09-25 15:00:36 [8  | REPORT]: [34|VaspWorkChain|results]: attaching the node Dict<38> as 'misc'
     2019-09-25 15:00:36 [9  | REPORT]: [34|VaspWorkChain|results]: attaching the node RemoteData<36> as 'remote_folder'
     2019-09-25 15:00:36 [10 | REPORT]: [34|VaspWorkChain|results]: attaching the node FolderData<37> as 'retrieved'
     2019-09-25 15:00:39 [11 | REPORT]: [34|VaspWorkChain|on_terminated]: cleaned remote folders of calculations: 35

7. What is stored on the VaspWorkChain::

     (aiida) max@workhorse:~$ verdi process show 34
     Property       Value
     -------------  ------------------------------------
     type           WorkChainNode
     pk             34
     uuid           8d513090-b45d-4b94-8810-b1db4f932a74
     label
     description
     ctime          2019-09-25 14:59:27.028678+00:00
     mtime          2019-09-25 15:00:37.091972+00:00
     process state  Finished
     exit status    0
     computer       [2] saga
     
     Called by      PK  Type
     -----------  ----  -------------
     CALL           33  WorkChainNode
     
     Inputs               PK  Type
     -----------------  ----  -------------
     clean_workdir        31  Bool
     code                  2  Code
     kpoints              23  KpointsData
     max_iterations       30  Int
     options              27  Dict
     parameters           24  Dict
     potential_family     25  Str
     potential_mapping    26  Dict
     settings             28  Dict
     structure            22  StructureData
     verbose              29  Bool
     
     Outputs          PK  Type
     -------------  ----  ----------
     misc             38  Dict
     remote_folder    36  RemoteData
     retrieved        37  FolderData
     
     Called      PK  Type
     --------  ----  -----------
     CALL        35  CalcJobNode
     
     Log messages
     ---------------------------------------------
     There are 7 log messages for this calculation
     Run 'verdi process report 34' to see them

8. In the outputs there is a ``misc``. This is a container for properties that does not depend on system size. Typically, total energies, band gaps, maximum forces, different tensors etc. Let us see what it contains::

     (aiida) max@workhorse:~$ verdi data dict show 38
     {
     "maximum_force": 0.0,
     "maximum_stress": 18.17613392,
     "symmetries": {
        "num_space_group_operations": {
            "dynamic": [
                48
            ],
            "static": [
                48
            ]
        },
        "point_group": {
            "dynamic": [
                "O_h"
            ],
            "static": [
                "O_h"
            ]
        },
        "primitive_translations": [
            1
        ]
     },
     "total_energies": {
     "energy_no_entropy": -10.79608481
     }
     }

9. You can control the output by adjusting the parameters related to the :ref:`parser`.
