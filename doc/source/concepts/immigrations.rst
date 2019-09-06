.. _immigrations:

Immigrations
============

Sometimes it is usefull or necessary to import `VASP`_ calculations that have not been executed in `AiiDA`_ or with the AiiDA-VASP plugin. For these cases we an existing `VASP`_ calculation need to be imported, or immigrated into the AiiDA-VASP realm. To make this process simple, we have designed am particular calculation that perform an immigration called :ref:`immigrator_calculation`.

Working principles
------------------
The immigration process is straightforward and utilizes the underlying calculation framework in `AiiDA`_.

1. First, we take the input files from the previous calculation and construct a new calculation
2. The calculation is executed, for instance through a workchain, but `VASP`_ is never run.
3. The output files are parsed similar to a normal execution of `VASP`_ using AiiDA-VASP.

This strategy is very convenient as it gives full data provenance.

Basic process
-------------

The following will assume that you have

 * `AiiDA`_ and AiiDA-VASP installed and a profile set up
 * added ``<computer>`` to your profile
 * added a `VASP`_ ``<code>`` (does not have to be runnable in principle) associated with an existing computer ``<computer>`` to your profile
 * a working understanding of how to run calculations in `AiiDA`_
 * a directory where `VASP`_ have previously been executed that you want to immigrate.

Assuming you ran `VASP`_ in ``/scratch/<user>/<runfolder>/`` on the remote (or local) computer ``<computer>``, the following python snippet will help you integrate that run into your `AiiDA`_ profile database::

   from aiida.orm import Code
   from aiida.plugins import CalculationFactory
   from aiida.engine import run

   # the folder you want to import
   calculation_folder = '/scratch/<user>/<runfolder>/'
   # first, get the code which the calculation will use
   code = Code.get_from_string('<code>@<computer>')  # use the name of your code and computer

   # then comes the information which can not be read from input files
   resources = {'num_machines': 1, 'num_mpiprocs_per_machine': 16}  # whatever is appropriate for you

   process, inputs = CalculationFactory('vasp.vasp').immigrant(code, calculation_folder, resources=resources)
   result = run(process, **inputs)

After running the above snippet, ``result`` will contain the output node dictionary of the created calculation.

Parsing additional information
------------------------------

In order to parse additional quantities, first, please see :ref:`parsing`. For the :ref:`immigrator_calculation`, this follows the same mechanism as for normal calculations. Here we take as an example, that the user want to also parse and include the electronic band structure (not enabled by default). We then make this available by passing the ``settings`` kwarg to the ``immigrant`` method::

   settings_dict = {'parser_settings': {'add_bands': True}}
   process, inputs = CalculationFactory('vasp.vasp').immigrant(..., settings=settings_dict)

And running the job as described above should produce a ``result`` containing the electronic band structure.

Missing POTCAR file
-------------------

Requires:

 * You have the POTCAR potentials to be used uploaded to your profile (see :ref:`potentials`)

In case you have removed the POTCAR files (e.g. to save space storing many runs), you need to specify which potcars to choose from the ones uploaded to your database::

   potential_mapping = {'Si': 'Si'}
   potential_family = 'PBE.54'
   process, inputs = CalculationFactory('vasp.vasp').immigrant(..., potential_family=potential_family, potential_mapping=potential_mapping)

This will use the chosen Si potential. It will fail if the POSCAR does not specify element names.

Additional input files
----------------------

For CHGCAR and/or WAVECAR it can not always be inferred from the INCAR, whether they are input or output files. If they are not automatically recognized as inputs, you may pass additional arguments to the `immigrant` method which support those inputs::

   process, inputs = CalculationFactory('vasp.vasp').immigrant(..., use_chgcar=True, use_wavecar=True)

which will create input nodes from these files without checking the INCAR file.

.. _AiiDA: https://www.aiida.net
.. _VASP: https://www.vasp.at
