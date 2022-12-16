.. _tutorial_fcc_si_step3:

=======================================
3. Equilibrium volume of FCC Si, take 2
=======================================

Here we will continue where the :ref:`previous tutorial<tutorial_fcc_si_step2>` left of.
In particular we will now create a workchain that execute the :ref:`vasp_workchain` (basically,
this executes `VASP`_) for each volume and stores it. In the end we will even use some Python
functions to locate the minimum and store it. We hope that you after this tutorial will
finally start to appreciate the convenience of `AiiDA-VASP`_. As before we assume you have
completed the steps of the :ref:`previous tutorial<tutorial_fcc_si_step2>`.

In order to make it as easy as possible, we have created a cookie cutter recipe for creating
`AiiDA-VASP`_ workchains. This is located in ``cookiecutters/workchain`` on the root of the
repository. If you have not cloned it, you can download only the cookiecutter folder using::

  %/$ svn checkout https://github.com/aiida-vasp/aiida-vasp/branches/develop/cookiecutters

What we want to do is to create a workchain that calculates some equation of state, in this
case over volume. For simplicity we will call the workchain ``Eos``. It is always recommended to
dissect the operations we will perform into parts before writing the workchain. In this way, it
is possibly to possibly reuse larger parts of the workchain at a later point. Also, it is easier
to understand what is going on from just inspecting the definitions of the inputs, steps and outputs
of the workchain. Please, have a look at :ref:`workchains` and references therein before continuing.

#. Make sure you have ``cookiecutter`` installed::

     %/$ pip install cookiecutter

#. Now we simply execute the cookiecutter using the ``workchains`` folder as the template::

     %/$ cookiecutter cookiecutters/workchains
     workchain_name [testchain]: eos
     workchain_docstring [This workchain is a testchain]: This workchain will accept a dictionary of structures and extract the total energies for each structure. The data is saved and the energy minimum is calculated and stored.
     next_workchain_to_call [vasp.vasp]:

#. A new folder ``eos`` was created which contains the ``eos.py`` workchain file. Inspect it.
   Maybe the most important part of a workchain is the ``define`` function, which, as its name
   implies, defines the workchain. For this particular case it looks like:

   .. literalinclude:: ../../../tutorials/eos_base.py
      :start-after: classmethod
      :end-before: def initialize

   We have defined some inputs:
   
   - All inputs defined in the ``_next_workchain``, which in this case is the :ref:`vasp_workchain`.
     We can include those very easily by using ``expose_inputs``. You can of course expose
     inputs from other workchains as long as `AiiDA`_ can find them.

   - A dictionary ``structures`` containing as many structures you want to run. Notice that
     the data type of each entry is :py:class:`StructureData<aiida.orm.StructureData>`.

   - Some generic exit codes.

   Then we have defined some outputs:

   - We attach all outputs from ``_next_workchain``, i.e. the :ref:`vasp_workchain` to this
     workchain (basically a link is created to avoid storing things twice). Similar to
     the inputs, ``expose_outputs`` can be leveraged.

   Finally, there is the ``outline`` which tells the workchain which functions to run. In
   this case we have:

   - A function that initializes what we need, ``initialize``.

   - A while loop that iterates over all the structures.

   - A ``init_next_workchain`` function that initializes what we need for the next workchain
     at the current iteration (here we make sure the specific structure is set).

   - A ``run_next_workchain`` which basically executes :ref:`vasp_workchain`.

   - A ``verify_next_workchain`` which verifies that there is a valid :ref:`vasp_workchain`
     and inherits any exit codes present from :ref:`vasp_workchain`.

   - A ``finalize`` which stores the results.

#. The workchain is not completely ready. We need to extract the total energies
   and specify the general workchain further. Please make the following changes
   to the generated workchain:

   .. literalinclude:: ../../../tutorials/eos.py
      :diff: ../../../tutorials/eos_base.py

   and save it as ``eos.py``. Or you could also download it::

     %/$ wget https://github.com/aiida-vasp/aiida-vasp/raw/develop/tutorials/eos.py

   The majority of changes were related to being more specific, except two things:

   - The necessity of decorating the function that generates the output array containing the
     volume and total energies in a ``calcfunction``. This is to preserve data provenance,
     otherwise we would not have known how the data was collected from each of the underlying
     workchains.

   - The inclusion of a ``calcfunction`` that interpolates the calculated data to find
     and ever better estimate of the volume at the energy minima. The example uses a
     cubic fit, which is certainly not very physical and should not be used in production.
     It is only to show how simply it is to leverage the power of Python, NumPy and SciPy.
     Again, this was decorated with a ``calcfunction`` in order to make sure `AiiDA`_ can
     preserve the data provenance.

#. Next, download the launch script that is tailored to launch the workchain we have now developed::

     %/$ wget https://github.com/aiida-vasp/aiida-vasp/raw/develop/tutorials/run_fcc_si_workchain.py

#. Change the ``options`` and ``code_string`` as you did :ref:`previously<tutorial_fcc_si_step1>`.

   .. warning:: Make sure you have activated your `AiiDA`_ virtual environment and
		that the `AiiDA`_ daemon is running before continuing.

#. Now we need to make sure the daemon can pick up the workchain. We can do this by
   making sure the daemon sees the directory where ``eos.py`` and ``run_fcc_si_workchain.py`` is
   located. The simplest approach is to add the following, to your virtual environment ``activate``
   script (assuming you do not use Conda)::

     $ echo "export PYTHONPATH=$PYTHONPATH:<yourdirectory>" >> ~/env/aiida-vasp/bin/activate

   assuming ``<yourdirectory>`` is the directory containing the ``eos.py`` and
   ``run_fcc_si_workchain.py`` files. The location of the ``activate`` is assumed from the
   previous steps in the tutorial. If you use Conda, please do::

     % echo "export PYTHONPATH=$PYTHONPATH:<yourdirectory>" >> $CONDA_PREFIX/etc/conda/activate.d/env_vars.sh
		
#. Submit the workchain by running the call script::

     %/$ python run_fcc_si_workchain.py

#. After a while we check the status::

     %/$ verdi process list -a
     PK  Created    Process label                 Process State      Process status
     ------  ---------  ----------------------------  -----------------  -----------------------------------------------------------
     103573  14m ago    EosWorkChain                  ⏹ Finished [0]
     103574  14m ago    VaspWorkChain                 ⏹ Finished [0]
     103575  14m ago    VaspCalculation               ⏹ Finished [0]
     103579  13m ago    VaspWorkChain                 ⏹ Finished [0]
     103580  13m ago    VaspCalculation               ⏹ Finished [0]
     103584  12m ago    VaspWorkChain                 ⏹ Finished [0]
     103585  12m ago    VaspCalculation               ⏹ Finished [0]
     103589  11m ago    VaspWorkChain                 ⏹ Finished [0]
     103590  11m ago    VaspCalculation               ⏹ Finished [0]
     103594  10m ago    VaspWorkChain                 ⏹ Finished [0]
     103595  10m ago    VaspCalculation               ⏹ Finished [0]
     103599  9m ago     VaspWorkChain                 ⏹ Finished [0]
     103600  9m ago     VaspCalculation               ⏹ Finished [0]
     103604  8m ago     VaspWorkChain                 ⏹ Finished [0]
     103605  8m ago     VaspCalculation               ⏹ Finished [0]
     103609  7m ago     VaspWorkChain                 ⏹ Finished [0]
     103610  7m ago     VaspCalculation               ⏹ Finished [0]
     103614  5m ago     VaspWorkChain                 ⏹ Finished [0]
     103615  5m ago     VaspCalculation               ⏹ Finished [0]
     103620  4m ago     store_total_energies          ⏹ Finished [0]
     103622  4m ago     locate_minimum                ⏹ Finished [0]

   As you can see, seven :ref:`vasp_calculation` was performed, one for each supplied volume.
   Also, there is a separate entry for the storrage of the total energies, which also performs
   a sort. In this was, all generated results are trackable and we have preserved data
   provenance.

#. Let us have a look at the output of ``EosWorkChain``::

     %/$ verdi process show 103573
     Property       Value
     -------------  ------------------------------------
     type           WorkChainNode
     pk             103573
     uuid           1055b7ac-e02e-4510-9ed8-62177a9c2bd1
     label
     description
     ctime          2019-10-08 12:50:46.739221+00:00
     mtime          2019-10-08 13:00:47.530855+00:00
     process state  Finished
     exit status    0
     computer       [6] mycluster
     
     Inputs              PK      Type
     ------------------  ------  -------------
     structures
         silicon_at_4_3  103563  StructureData
         silicon_at_4_2  103562  StructureData
         silicon_at_4_1  103561  StructureData
         silicon_at_4_0  103560  StructureData
         silicon_at_3_9  103559  StructureData
         silicon_at_3_8  103558  StructureData
         silicon_at_3_7  103557  StructureData
         silicon_at_3_6  103556  StructureData
         silicon_at_3_5  103555  StructureData
     clean_workdir       103572  Bool
     code                101271  Code
     kpoints             103564  KpointsData
     max_iterations      103571  Int
     options             103568  Dict
     parameters          103565  Dict
     potential_family    103566  Str
     potential_mapping   103567  Dict
     settings            103569  Dict
     verbose             103570  Bool
     
     Outputs          PK  Type
     -----------  ------  ---------
     eos          103621  ArrayData
     eos_minimum  103623  Dict
     
     Called        PK  Type
     --------  ------  ----------------
     CALL      103622  CalcFunctionNode
     CALL      103620  CalcFunctionNode
     CALL      103614  WorkChainNode
     CALL      103609  WorkChainNode
     CALL      103604  WorkChainNode
     CALL      103599  WorkChainNode
     CALL      103594  WorkChainNode
     CALL      103589  WorkChainNode
     CALL      103584  WorkChainNode
     CALL      103579  WorkChainNode
     CALL      103574  WorkChainNode
     
     Log messages
     ---------------------------------------------
     There are 9 log messages for this calculation
     Run 'verdi process report 103573' to see them

     
#. Inspect the total energies versus volume::

     %/$ verdi data array show 103621
     {
     "eos": [
        [
            10.71875,
            -4.42341934
        ],
        [
            11.664,
            -4.66006377
        ],
        [
            12.66325,
            -4.79595549
        ],
        [
            13.718,
            -4.86303425
        ],
        [
            14.82975,
            -4.87588342
        ],
        [
            16.0,
            -4.8481406
        ],
        [
            17.23025,
            -4.78451894
        ],
        [
            18.522,
            -4.69228806
        ],
        [
            19.87675,
            -4.58122037
        ]
      ]
     }

#. And the located minimum::

     %/$ verdi data dict show 103623
     {
         "energy": -4.8769539175,
         "volume": 14.55935661641
     }

That concludes this tutorial. We hope at this point you have now realized
that `AiiDA-VASP`_ seems somewhat usefull and that you would like to continue to
learn more, maybe even start to write your own :ref:`workflows` or :ref:`workchains`.
   
.. _Gnuplot: http://gnuplot.info/
.. _AiiDA: https://www.aiida.net
.. _FCC Si: https://cms.mpi.univie.ac.at/wiki/index.php/Fcc_Si
.. _VASP: https://www.vasp.at
.. _AiiDA-VASP: https://github.com/aiida-vasp/aiida-vasp
