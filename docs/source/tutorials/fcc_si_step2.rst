.. _tutorial_fcc_si_step2:

===============================
2. Equilibrium volume of FCC Si
===============================

Here we will follow the second part of the `FCC Si`_ example in the `VASP`_ tutorial,
i.e. perform the calculations of the total energies at different volumes. Please make
sure you have completed the :ref:`previous tutorial<tutorial_fcc_si_step1>` before continuing.

#. When the standard tutorial is completed we will now do the same in `AiiDA-VASP`_ using
   different strategies to you get a feel for how you can structurize a simple workflow
   like this. Let us now fetch the `AiiDA-VASP`_ run file for this example::

     %/$ wget https://github.com/aiida-vasp/aiida-vasp/raw/develop/tutorials/run_fcc_si_multiple_volumes.py

#. Inspect the file, which has the following changes compared to the file used in the :ref:`previous tutorial<tutorial_fcc_si_step1>`

   .. literalinclude:: ../../../tutorials/run_fcc_si_multiple_volumes.py
      :diff: ../../../tutorials/run_fcc_si_one_volume.py

   As you can see, we did only minimal changes, all which should be self explanatory.

#. Change the ``options`` and ``code_string`` as you did in the :ref:`previous tutorial<tutorial_fcc_si_step1>`.

   .. warning:: Make sure you have activated your `AiiDA`_ virtual environment and
      that the `AiiDA`_ daemon is running before continuing.

#. Save and execute the resulting run script by issuing::

     %/$ python run_fcc_si_multiple_volumes.py

#. Let us check the progress::

     %/$ verdi process list -a
     PK  Created    Process label                 Process State      Process status
     ------  ---------  ----------------------------  -----------------  -----------------------------------------------------------
     101383  22s ago    VerifyWorkChain               ⏵ Waiting          Waiting for child processes: 101393
     101392  22s ago    VerifyWorkChain               ⏵ Waiting          Waiting for child processes: 101403
     101393  21s ago    VaspWorkChain                 ⏵ Waiting          Waiting for child processes: 101413
     101402  21s ago    VerifyWorkChain               ⏵ Waiting          Waiting for child processes: 101423
     101403  21s ago    VaspWorkChain                 ⏵ Waiting          Waiting for child processes: 101424
     101412  21s ago    VerifyWorkChain               ⏵ Waiting          Waiting for child processes: 101438
     101413  21s ago    VaspCalculation               ⏵ Waiting          Waiting for transport task: upload
     101422  21s ago    VerifyWorkChain               ⏵ Waiting          Waiting for child processes: 101444
     101423  21s ago    VaspWorkChain                 ⏵ Waiting          Waiting for child processes: 101454
     101424  20s ago    VaspCalculation               ⏵ Waiting          Waiting for transport task: upload
     101433  20s ago    VerifyWorkChain               ⏵ Waiting          Waiting for child processes: 101465
     101443  20s ago    VerifyWorkChain               ⏵ Waiting          Waiting for child processes: 101467
     101438  20s ago    VaspWorkChain                 ⏵ Waiting          Waiting for child processes: 101464
     101444  20s ago    VaspWorkChain                 ⏵ Waiting          Waiting for child processes: 101466
     101453  20s ago    VerifyWorkChain               ⏵ Waiting          Waiting for child processes: 101468
     101454  19s ago    VaspCalculation               ⏵ Waiting          Waiting for transport task: upload
     101463  19s ago    VerifyWorkChain               ⏵ Waiting          Waiting for child processes: 101469
     101465  19s ago    VaspWorkChain                 ⏵ Waiting          Waiting for child processes: 101471
     101464  19s ago    VaspCalculation               ⏵ Waiting          Waiting for transport task: upload
     101466  18s ago    VaspCalculation               ⏵ Waiting          Waiting for transport task: upload
     101467  18s ago    VaspWorkChain                 ⏵ Waiting          Waiting for child processes: 101472
     101468  18s ago    VaspWorkChain                 ⏵ Waiting          Waiting for child processes: 101470
     101469  17s ago    VaspWorkChain                 ⏵ Waiting          Waiting for child processes: 101473
     101470  17s ago    VaspCalculation               ⏵ Waiting          Waiting for transport task: upload
     101471  17s ago    VaspCalculation               ⏵ Waiting          Waiting for transport task: submit
     101472  16s ago    VaspCalculation               ⏵ Waiting          Waiting for transport task: upload
     101473  16s ago    VaspCalculation               ⏵ Waiting          Waiting for transport task: submit

   All calculations are submitted to the cluster in parallel, which is pretty convenient and
   certainly an improvement compared to the regular tutorial, where you do the calculations
   sequentially. Not a problem for a few calculations, but maybe not ideal if you want to execute
   a few thousand calculations.

#. After a while, we execut ``verdi process list -a`` and get::

     %/$ verdi process list -a
     PK  Created    Process label                 Process State      Process status
     ------  ---------  ----------------------------  -----------------  -----------------------------------------------------------
     101383  2m ago     VerifyWorkChain               ⏹ Finished [0]
     101392  2m ago     VerifyWorkChain               ⏹ Finished [0]
     101393  2m ago     VaspWorkChain                 ⏹ Finished [0]
     101402  2m ago     VerifyWorkChain               ⏹ Finished [0]
     101403  2m ago     VaspWorkChain                 ⏹ Finished [0]
     101412  2m ago     VerifyWorkChain               ⏹ Finished [0]
     101413  2m ago     VaspCalculation               ⏹ Finished [0]
     101422  2m ago     VerifyWorkChain               ⏹ Finished [0]
     101423  2m ago     VaspWorkChain                 ⏹ Finished [0]
     101424  2m ago     VaspCalculation               ⏹ Finished [0]
     101433  2m ago     VerifyWorkChain               ⏹ Finished [0]
     101443  2m ago     VerifyWorkChain               ⏹ Finished [0]
     101438  2m ago     VaspWorkChain                 ⏹ Finished [0]
     101444  2m ago     VaspWorkChain                 ⏹ Finished [0]
     101453  2m ago     VerifyWorkChain               ⏹ Finished [0]
     101454  2m ago     VaspCalculation               ⏹ Finished [0]
     101463  2m ago     VerifyWorkChain               ⏹ Finished [0]
     101465  2m ago     VaspWorkChain                 ⏹ Finished [0]
     101464  2m ago     VaspCalculation               ⏹ Finished [0]
     101466  2m ago     VaspCalculation               ⏹ Finished [0]
     101467  2m ago     VaspWorkChain                 ⏹ Finished [0]
     101468  2m ago     VaspWorkChain                 ⏹ Finished [0]
     101469  2m ago     VaspWorkChain                 ⏹ Finished [0]
     101470  2m ago     VaspCalculation               ⏹ Finished [0]
     101471  2m ago     VaspCalculation               ⏹ Finished [0]
     101472  2m ago     VaspCalculation               ⏹ Finished [0]
     101473  2m ago     VaspCalculation               ⏹ Finished [0]

   All processes are in a finished state and we can extract the total energies for each step.
   However, it should be obvious that extracting the total energies from ``misc`` from each
   step manually seems to be a waste of time. One could query (AiiDA has shortcuts for this), but
   one would still need to look up the values in some way. Maybe it is easier to be able to access
   them directly when all the calculations are complete? Let us do another modification to the
   run script above, namely:

   .. literalinclude:: ../../../tutorials/run_fcc_si_multiple_volumes_eos.py
      :diff: ../../../tutorials/run_fcc_si_multiple_volumes.py

Save the new file as ``run_fcc_si_multiple_volumes_eos.py`` or fetch it with::

     wget https://github.com/aiida-vasp/aiida-vasp/raw/develop/tutorials/run_fcc_si_multiple_volumes_eos.py

   This file should enable automatic extraction of total energies versus the different volumes.

#. Execute ``run_fcc_si_multiple_volumes_eos.py``::

     %/$ python run_fcc_si_multiple_volumes_eos.py

   Give the script time to complete. Contrary to before you will now see the output
   in the terminal (the similar output we previously got with ``verdi process report`` for
   a given process). More on this later. A file ``eos`` should be generated.

#. Inspect the ``eos`` file that was generated::

     %/$ more eos
     3.5 -4.42341934
     3.6 -4.66006377
     3.7 -4.79595549
     3.8 -4.86303425
     3.9 -4.87588342
     4.0 -4.8481406
     4.1 -4.78451894
     4.2 -4.69228806
     4.3 -4.58122037

#. Plot the data with your favorite plotting tool, for instance `Gnuplot`_::

     %/$ gnuplot
     gnuplot> plot "eos" with lp

   which should give you something similar that is shown in the `FCC Si`_ tutorial.

We have now completed the `FCC Si`_ tutorial in `AiiDA-VASP`_. However, maybe not in the
best way. As you might already have realized, we had to change ``submit`` to ``run`` in order
to make sure we got the output of each executed workchain. If we would have used ``submit``
we would only have access to a link to some output that we did not know was occupied or not at the
point of inspection (unless we waited until we know the execution of the workchain
was complete). In doing so, all the volume steps now run sequentially and thus it takes
far longer for all the volume calculations to be finished (assuming your cluster started
each volume calculations at the time of submission).

In the next tutorial we will instead develop a workchain to perform the different
volume calculations (in fact, more general, namely different structure calculations)
and call that from the run script. In doing so, you will realize that we can keep
the benefits of automatic extraction of the data and parallel execution. In addition,
we get the benefit of storing the result in the database.

.. _Gnuplot: http://gnuplot.info/
.. _AiiDA: https://www.aiida.net
.. _FCC Si: https://cms.mpi.univie.ac.at/wiki/index.php/Fcc_Si
.. _VASP: https://www.vasp.at
.. _AiiDA-VASP: https://github.com/aiida-vasp/aiida-vasp
