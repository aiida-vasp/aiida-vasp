.. _tutorial_fcc_si_step4:

=======================
 4. Running in parallel
=======================

In this tutorial, we will look into how to enable parallel execution of
say multiple workchains, in this case for each volume.
Make sure you have completed the steps of the :ref:`previous tutorial<tutorial_fcc_si_step3>`.

Even though the previous implementation of the ``EosWorkChain`` was easy to comprehend and
clearly written, it was not possible to submit all the :ref:`vasp_workchain` in one go. The reason
for this is that we did not perform a simple iteration over the ``self.submit`` and ``self.to_context``
and the iteration included other methods in the ``define`` section composed in a ``_while`` statement.

We will now show how we can modify the previous workchain to be able to calculate the total energies for
each volume in parallel.

#. First modify the previous ``eos.py`` accordingly:

   .. literalinclude:: ../../../tutorials/eos.py
      :diff: ../../../tutorials/eos_parallel.py

   or download it::

     $ wget https://github.com/aiida-vasp/aiida-vasp/raw/master/tutorials/eos_parallel.py

   make sure to place it in the ``PYTHONPATH`` added in the previous tutorial.

#. Replace ``EosWorkChain`` with ``EosParallelWorkChain`` and ``eos`` with ``eos_parallel`` in ``run_fcc_si_workchain.py``.

#. Restart the daemon with ``verdi daemon restart``.

   .. warning::
      Not restarting the daemon is a very common mistake when debugging and updating code.
      Any updates to code related to AiiDA should be accompanied by a daemon restart so that
      it can also pick up the updated code.

#. Submit the workchain by running the call script::

     $ python run_fcc_si_workchain.py

#. Check the status quickly to verify that it indeed now started all the ``VaspWorkChain``::

     $ verdi process list
       PK  Created    Process label         Process State     Process status
     ----  ---------  --------------------  ----------------  --------------------------------------------------------------------------------------
     2239  14s ago    EosParallelWorkChain  ⏵ Waiting         Waiting for child processes: 2240, 2241, 2242, 2245, 2248, 2251, 2254, 2257, 2260
     2240  13s ago    VaspWorkChain         ⏵ Waiting         Waiting for child processes: 2244
     2241  13s ago    VaspWorkChain         ⏵ Waiting         Waiting for child processes: 2247
     2242  12s ago    VaspWorkChain         ⏵ Waiting         Waiting for child processes: 2250
     2244  12s ago    VaspCalculation       ⏵ Waiting         Waiting for transport task: upload
     2245  12s ago    VaspWorkChain         ⏵ Waiting         Waiting for child processes: 2253
     2247  12s ago    VaspCalculation       ⏵ Waiting         Waiting for transport task: upload
     2248  11s ago    VaspWorkChain         ⏵ Waiting         Waiting for child processes: 2256
     2250  11s ago    VaspCalculation       ⏵ Waiting         Waiting for transport task: upload
     2251  11s ago    VaspWorkChain         ⏵ Waiting         Waiting for child processes: 2259
     2253  11s ago    VaspCalculation       ⏵ Waiting         Waiting for transport task: upload
     2254  10s ago    VaspWorkChain         ⏵ Waiting         Waiting for child processes: 2262
     2256  10s ago    VaspCalculation       ⏵ Waiting         Waiting for transport task: upload
     2257  10s ago    VaspWorkChain         ⏵ Waiting         Waiting for child processes: 2264
     2260  9s ago     VaspWorkChain         ⏵ Waiting         Waiting for child processes: 2266
     2259  9s ago     VaspCalculation       ⏵ Waiting         Waiting for transport task: upload
     2262  9s ago     VaspCalculation       ⏵ Waiting         Waiting for transport task: upload
     2264  8s ago     VaspCalculation       ⏵ Waiting         Waiting for transport task: upload
     2266  8s ago     VaspCalculation       ⏵ Waiting         Waiting for transport task: upload

     Total results: 19

     Report: last time an entry changed state: 7s ago (at 20:04:23 on 2022-12-21)
     Report: Using 1% of the available daemon worker slots.

#. Wait a bit and check again::

     $ verdi process list -a
       PK  Created    Process label         Process State     Process status
     ----  ---------  --------------------  ----------------  --------------------------------------------------------------------------------------
     2362  12h ago    EosParallelWorkChain  ⏹ Finished [0]
     2363  12h ago    VaspWorkChain         ⏹ Finished [0]
     2364  12h ago    VaspWorkChain         ⏹ Finished [0]
     2365  12h ago    VaspWorkChain         ⏹ Finished [0]
     2366  12h ago    VaspWorkChain         ⏹ Finished [0]
     2367  12h ago    VaspWorkChain         ⏹ Finished [0]
     2370  12h ago    VaspWorkChain         ⏹ Finished [0]
     2369  12h ago    VaspCalculation       ⏹ Finished [0]
     2371  12h ago    VaspWorkChain         ⏹ Finished [0]
     2373  12h ago    VaspCalculation       ⏹ Finished [0]
     2375  12h ago    VaspCalculation       ⏹ Finished [0]
     2376  12h ago    VaspWorkChain         ⏹ Finished [0]
     2378  12h ago    VaspCalculation       ⏹ Finished [0]
     2379  12h ago    VaspWorkChain         ⏹ Finished [0]
     2381  12h ago    VaspCalculation       ⏹ Finished [0]
     2383  12h ago    VaspCalculation       ⏹ Finished [0]
     2385  12h ago    VaspCalculation       ⏹ Finished [0]
     2387  12h ago    VaspCalculation       ⏹ Finished [0]
     2389  12h ago    VaspCalculation       ⏹ Finished [0]
     2418  12h ago    store_total_energies  ⏹ Finished [0]
     2420  12h ago    locate_minimum        ⏹ Finished [0]

     Total results: 21

     Report: last time an entry changed state: 12h ago (at 20:15:52 on 2022-12-21)
     Report: Using 0% of the available daemon worker slots.

   By running these in parallel, it took, given we used the same computational cluster, and that all jobs started,
   significantly less time overall to complete the workflow. Very often this is a great way to manage many calculations
   simultaneously. You can inspect the same content as in the previous tutorial and conclude it gives the same output.

   .. note::
      Sometimes the cluster administrators see the activity of AiiDA as too intense, either due to many
      SSH connections opening at the same time or because the queue is filling up. One can then consider to
      submit in batches and/or limit the SSH connections. If the latter, please have a look at the
      `connection overload`_ documentation.

.. _Gnuplot: http://gnuplot.info/
.. _AiiDA: https://www.aiida.net
.. _tutorial for writing workflows: https://aiida.readthedocs.io/projects/aiida-core/en/latest/intro/tutorial.html#workflows
.. _AiiDA documentation: https://aiida.readthedocs.io/projects/aiida-core/en/latest/index.html
.. _FCC Si: https://cms.mpi.univie.ac.at/wiki/index.php/Fcc_Si
.. _VASP: https://www.vasp.at
.. _AiiDA-VASP: https://github.com/aiida-vasp/aiida-vasp
.. _connection overload: https://aiida.readthedocs.io/projects/aiida-core/en/latest/howto/run_codes.html#mitigating-connection-overloads
