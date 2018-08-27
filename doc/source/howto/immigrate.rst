How to import non-AiiDA VASP runs
=================================

Basic process
-------------

The following will assume that you have

 * AiiDA and AiiDA-VASP installed and a profile set up
 * added ``<computer>`` to your profile
 * added a VASP code (does not have to be runnable in principle) with computer ``<computer>`` to your profile
 * a working understanding of how to run calculations in AiiDA

See :ref:`main-quickstart` for an overview of setting up AiiDA-VASP and :ref:`main-running` for a quick intro to running calculations.

Assuming you ran VASP in ``/scratch/<user>/<runfolder>/`` on the remote (or local) computer ``<computer>``, the following python snippets will help you integrate that run into your AiiDA profile database::

   from aiida.orm import CalculationFactory, Code
   from aiida.work import run

   calculation_folder = '/scratch/<user>/<runfolder>/'
   # first, get the code which the calculation will use
   code = Code.get_from_string('vasp@<computer>')  # use the name of your code

   # then comes the information which can not be read from input files
   resources = {'num_machines': 1, 'num_mpiprocs_per_machine': 16}  # whatever is appropriate for you

   # replace `vasp.vasp` with what best fits your run
   process, inputs = CalculationFactory('vasp.vasp').immigrant(code, calculation_folder, resources=resources)
   result = run(process, **inputs)

After running the above snippet, ``result`` will contain the output node dictionary of the created calculation.

Parsing additional information
------------------------------

In order to parse additional quantities, the same mechanism as for normal calculations is available by passing the ``settings`` kwarg to the ``immigrant`` method::

   settings_dict = {'parser_settings': {'add_bands': True}}
   process, inputs = CalculationFactory('vasp.vasp').immigrant(..., settings=settings_dict)

Missing POTCAR file
-------------------

Requires:

 * You have the POTCAR potentials to be used uploaded to your profile (see :ref:`howto-potcar`)

In case you have removed the POTCAR file (e.g. to save space storing many runs), you need to specify which potcars to chose from the ones uploaded to your DB::

   potcar_spec = {
      'family': 'PBE_54',
      'map': {'Si': 'Si'}
   }
   process, inputs = CalculationFactory('vasp.vasp').immigrant(..., potcar_spec=potcar_spec)

This will use the chosen Si potential, assuming the POSCAR specifies only Si atoms. It will fail if the POSCAR does not specify element names.

Additional input files
----------------------

For CHGCAR and / or WAVECAR it can not always be inferred from INCAR, whether they are input or output files. If they are not automatically recognized as inputs, you may pass additional arguments to the `immigrant` method of vasp calculations which support those inputs::

   process, inputs = CalculationFactory('vasp.vasp').immigrant(..., use_chgcar=True, use_wavecar=True)

Will create input nodes from these files without checking the INCAR file.
