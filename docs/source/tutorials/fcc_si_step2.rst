.. _tutorial_fcc_si_step2:

=============
2. Iterations
=============

Here we continue from the :ref:`previous tutorial<tutorial_fcc_si_step1>` performing the initial setup and test run for a static volume for the `FCC Si`_
example in the `VASP tutorials`_. Here we will perform the calculations of the total energies at different volumes. Please make
sure you have completed the :ref:`previous tutorial<tutorial_fcc_si_step1>` before continuing.

#. Let us first fetch the `AiiDA-VASP`_ run file for the multiple volumes::

     $ wget https://github.com/aiida-vasp/aiida-vasp/raw/master/tutorials/run_fcc_si_multiple_volumes.py

#. Inspect the file, which has the following changes compared to the file used in the :ref:`previous tutorial<tutorial_fcc_si_step1>`:

   .. literalinclude:: ../../../tutorials/run_fcc_si_multiple_volumes.py
      :diff: ../../../tutorials/run_fcc_si_one_volume.py

   As you can see, we did only minimal changes, all which should be self explanatory. In essence we perform the same
   workchain as before, only that we now iterate over multiple workchains, each getting supplied a different volume of
   the unit cell.

#. Change the ``options`` and ``code_string`` as you did in the :ref:`previous tutorial<tutorial_fcc_si_step1>` to
   suit your setup.

   .. warning:: Make sure you have activated your `AiiDA`_ virtual environment and
      that the `AiiDA`_ daemon is running before continuing.

#. Save and execute the resulting run script by issuing::

     $ python run_fcc_si_multiple_volumes.py

#. Let us check the progress::

     $ verdi process list -a
       PK  Created    Process label    Process State     Process status
     ----  ---------  ---------------  ----------------  -----------------------------------
     1043  37s ago    VaspWorkChain    ⏵ Waiting         Waiting for child processes: 1067
     1054  37s ago    VaspWorkChain    ⏵ Waiting         Waiting for child processes: 1102
     1065  36s ago    VaspWorkChain    ⏵ Waiting         Waiting for child processes: 1132
     1067  36s ago    VaspCalculation  ⏵ Waiting         Waiting for transport task: upload
     1078  36s ago    VaspWorkChain    ⏵ Waiting         Waiting for child processes: 1139
     1089  36s ago    VaspWorkChain    ⏵ Waiting         Waiting for child processes: 1141
     1100  35s ago    VaspWorkChain    ⏵ Waiting         Waiting for child processes: 1143
     1102  35s ago    VaspCalculation  ⏵ Waiting         Waiting for transport task: upload
     1113  35s ago    VaspWorkChain    ⏵ Waiting         Waiting for child processes: 1145
     1124  35s ago    VaspWorkChain    ⏵ Waiting         Waiting for child processes: 1147
     1137  34s ago    VaspWorkChain    ⏵ Waiting         Waiting for child processes: 1149
     1132  34s ago    VaspCalculation  ⏵ Waiting         Waiting for transport task: upload
     1139  33s ago    VaspCalculation  ⏵ Waiting         Waiting for transport task: upload
     1141  32s ago    VaspCalculation  ⏵ Waiting         Waiting for transport task: upload
     1143  31s ago    VaspCalculation  ⏵ Waiting         Waiting for transport task: upload
     1145  30s ago    VaspCalculation  ⏵ Waiting         Waiting for transport task: upload
     1147  29s ago    VaspCalculation  ⏵ Waiting         Waiting for transport task: upload
     1149  27s ago    VaspCalculation  ⏵ Waiting         Waiting for transport task: upload

     Total results: 18

     Report: last time an entry changed state: 26s ago (at 09:24:10 on 2022-12-20)
     Report: Using 9% of the available daemon worker slots.

   All calculations are submitted to the cluster in parallel, which is pretty convenient and
   certainly an improvement compared to the regular tutorial, where you do the calculations
   sequentially. Not a problem for a few calculations, but maybe not ideal if you want to execute
   a few thousand calculations.

#. After a while, we again execute ``verdi process list -a``, hoping our runs are finished and get::

     $ verdi process list -a
       PK  Created    Process label    Process State     Process status
     ----  ---------  ---------------  ----------------  -----------------------------------
     1043  4m ago     VaspWorkChain    ⏹ Finished [0]
     1054  4m ago     VaspWorkChain    ⏹ Finished [0]
     1065  4m ago     VaspWorkChain    ⏹ Finished [0]
     1067  4m ago     VaspCalculation  ⏹ Finished [0]
     1078  4m ago     VaspWorkChain    ⏹ Finished [0]
     1089  4m ago     VaspWorkChain    ⏹ Finished [0]
     1100  4m ago     VaspWorkChain    ⏹ Finished [0]
     1102  4m ago     VaspCalculation  ⏹ Finished [0]
     1113  4m ago     VaspWorkChain    ⏹ Finished [0]
     1124  4m ago     VaspWorkChain    ⏹ Finished [0]
     1137  4m ago     VaspWorkChain    ⏹ Finished [0]
     1132  4m ago     VaspCalculation  ⏹ Finished [0]
     1139  4m ago     VaspCalculation  ⏹ Finished [0]
     1141  3m ago     VaspCalculation  ⏹ Finished [0]
     1143  3m ago     VaspCalculation  ⏹ Finished [0]
     1145  3m ago     VaspCalculation  ⏹ Finished [0]
     1147  3m ago     VaspCalculation  ⏹ Finished [0]
     1149  3m ago     VaspCalculation  ⏹ Finished [0]

     Total results: 18

     Report: last time an entry changed state: 30s ago (at 09:27:33 on 2022-12-20)
     Report: Using 0% of the available daemon worker slots.

   All processes are in a finished state and we can extract the total energies for each step.
   However, it should be obvious that extracting the total energies from ``misc`` from each
   step manually seems a bit inefficient. One could query (`AiiDA`_ has shortcuts for this), but
   one would still need to look up the values in some way. Surely the computer should be able
   to eject what we need, right? Maybe it is easier to be able to access
   them directly when all the calculations are complete? Let us do another modification to the
   run script above, namely:

   .. literalinclude:: ../../../tutorials/run_fcc_si_multiple_volumes_eos.py
      :diff: ../../../tutorials/run_fcc_si_multiple_volumes.py

   Notice that except for the iteration over the call to the ``main`` which eventually runs
   the workchain, we have replaced ``submit`` with ``run``. The reason for this is that we have to
   make sure the VASP run completes before we can extract the total energies. This also has the consequence
   that the VASP calculations are now executed sequentially, instead of in parallel. We will later suggest
   how we can fix this, but for the time being, let us go with this approach. Save the new file as
   ``run_fcc_si_multiple_volumes_eos.py`` or fetch it with::

     wget https://github.com/aiida-vasp/aiida-vasp/raw/master/tutorials/run_fcc_si_multiple_volumes_eos.py

   This file should enable automatic extraction of total energies versus the different volumes.

#. Execute ``run_fcc_si_multiple_volumes_eos.py``::

     $ python run_fcc_si_multiple_volumes_eos.py

   Give the script time to complete. Contrary to before you will now see the output
   in the terminal (the similar output we previously got with ``verdi process report`` for
   a given process). This is because we now use ``run`` instead of ``submit``.
   More on this later. A file ``eos`` should be generated in the folder you executed the
   Python command.

#. Inspect the ``eos`` file that was generated::

     $ more eos
     3.5 -4.42341939
     3.6 -4.66006381
     3.7 -4.79595554
     3.8 -4.86303429
     3.9 -4.87588353
     4.0 -4.8481407
     4.1 -4.78451926
     4.2 -4.69228837
     4.3 -4.58122058

#. Plot the data with your favorite plotting tool, for instance `Gnuplot`_::

     $ gnuplot
     gnuplot> plot "eos" with lp

   which should give you something similar that is shown in the `FCC Si`_ tutorial.

We have now completed the `FCC Si`_ tutorial in `AiiDA-VASP`_. However, maybe not in the
best way. As you have already realized, we had to change ``submit`` to ``run`` in order
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
.. _FCC Si: https://www.vasp.at/wiki/index.php/Fcc_Si
.. _VASP: https://www.vasp.at
.. _VASP tutorials: https://www.vasp.at/wiki/index.php/Category:Tutorials
.. _AiiDA-VASP: https://github.com/aiida-vasp/aiida-vasp
