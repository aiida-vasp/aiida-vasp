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

#. Make sure you have cookiecutter installed::

     %/$ pip install cookiecutter

.. warning:: Make sure you have activated your `AiiDA`_ virtual environment and
	     that the `AiiDA`_ daemon is running before continuing.

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

#. The workchain is not completely ready. We need to extract the total energies and specify
the general workchain further. Please make the following changes to the generated workchain::

.. literalinclude:: ../../../tutorials/eos_base.py
   :diff: ../../../tutorials/eos.py

   Or download it::

     %/$ wget https://github.com/aiida-vasp/aiida-vasp/raw/develop/tutorials/eos.py

#. Download the launch script that is tailored to launch the workchain we have now developed::

     %/$ wget https://github.com/aiida-vasp/aiida-vasp/raw/develop/tutorials/run_fcc_si_workchain.py

#. Change the ``options`` and ``code_string`` as you did :ref:`previously<tutorial_fcc_si_step1>`.

.. warning:: Make sure you have activated your `AiiDA`_ virtual environment and                                                                                                                                                  that the `AiiDA`_ daemon is running before continuing.

#. Submit the workchain by running the call script::

     %/$ python run_fcc_si_workchain.py

.. _Gnuplot: http://gnuplot.info/
.. _AiiDA: https://www.aiida.net
.. _FCC Si: https://cms.mpi.univie.ac.at/wiki/index.php/Fcc_Si
.. _VASP: https://www.vasp.at
.. _AiiDA-VASP: https://github.com/aiida-vasp/aiida-vasp
