.. _bulk_modulus_script:

=======================
6. Designing a workflow
=======================

This tutorial focuses on how to approach designing, developing, and finally launching a workflow.

We will use as an example, the calculation of the bulk modulus of wurtzite-type SiC.

The tutorial is divided in two sections. In the first part we will prepare the mental state by identifying
the different steps we need to perform in order to calculate the bulk modulus using VASP.
While, in the second part we will make a simple script that depend on
`AiiDA`_ and `AiiDA-VASP`_ to perform the same calculation.

Also, note that we will in this tutorial try to motivate the reader to write most of this and not
just downloading the scripts and executing them. This is good training practice. However, in case you need or want to
there is also a link to the scripts used in this tutorial. But please try not to be tempted to download them
right away. For this reason, they are also mentioned at the end.

1. Flow of work
---------------

.. _workflow_bulk_modulus:

Let us first start by trying to sketch the outline of the steps needed to investigate a property,
phenomena or something else. These steps typically then define the flow of work to be performed, or
the workflow.

Here we take a rather simple example to illustrate how to approach this problem, such that it becomes
easier to write up the workflow in `AiiDA`_ later. Maybe a good mental picture of this would be that
we want to bridge the gap between the two concepts of *scientific workflows* and the *workflows as code*.
The first being how you would approach a problem from the scientific side, which is obviously not something
you can directly do on a computer, while the second is how it can be performed in a computer.

To obtain the bulk modulus we typically need to perform several different calculations.
Broadly, we would like to follow the following steps:

#. First we need to prepare the inputs, or starting point. In this case the crystal structure is the most
   important. Here, this would be the wurtzite-type SiC.

#. Then we need to `Relax` the crystal structure.

#. Then we need to `Wait` until (2) finishes and verify results.

#. Then we need to `Create two structures` at fixed volumes with +/- 1% from the relaxed
   structure obtained at the step (3).

#. Then we need to `Relax the shape` of the structures generated in step (4).

#. Then we need to `Wait` until (5) finishes and verify results.

#. Then we need to `Compute bulk modulus` as a post process by :math:`K \simeq -V_0\frac{\Delta P}{\Delta  V}`,
   where we use the previous outputs as inputs.

Let us now try to perform these steps exclusively using `VASP`_. The steps above are of course general and could
be executed on any code that provides the functionality to honor the sub-steps.

2. Manual calculation
---------------------

1. Preparing input files
^^^^^^^^^^^^^^^^^^^^^^^^

#. We will utilize the following crystal structure for wurtzite-type SiC, for `VASP`_ this means using
   the following ``POSCAR`` file::

     wurtzite-type SiC
       1.0000000000
       3.0920000000   0.0000000000   0.0000000000
      -1.5460000000   2.6777505485   0.0000000000
       0.0000000000   0.0000000000   5.0730000000
     Si    C
	 2     2
     Direct
       0.3333333333   0.6666666667   0.0000000000
       0.6666666667   0.3333333333   0.5000000000
       0.3333333333   0.6666666667   0.3758220000
       0.6666666667   0.3333333333   0.8758220000

#. We also need to prepare a suitable ``KPOINTS`` file, here::

     # Half grid shift along c*
     0
     Gamma
		 6             6             4
       0.000000000   0.000000000   0.500000000

#. Then, we need the ``INCAR`` file to tell `VASP`_ to relax the structure::

     EDIFF = 1e-08
     EDIFFG = -1e-05
     ENCUT = 500
     IALGO = 38
     IBRION = 2
     ISIF = 3
     ISMEAR = 0
     LCHARG = .FALSE.
     LREAL = .FALSE.
     LWAVE = .FALSE.
     NELM = 100
     NSW = 100
     NELMIN = 5
     PREC = Accurate
     SIGMA = 0.01

#. You also need to assemble the potential files for ``SiC`` by concatenating
   the recommended ``PBE.54`` potential files for ``Si`` and ``C`` into a suitable ``POTCAR``.
   If in doubt about these steps, please have a look at `VASP lectures`_, `VASP tutorials`_, `VASP howtos`_,
   `VASP tutorials using notebooks`_ or `VASP videos`_ or ask experienced `VASP`_ users. `AiiDA-VASP`_ is not
   a substitute for not knowing how to use `VASP`_, its intent is to make it more efficient, reproducible to use `VASP`_,
   in addition to open up for new areas of applications.

#. Create a jobfile, or another way to execute `VASP`_ in your environment and continue to the next step.

2. Perform the relaxation and 3. Wait for it to finish
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

#. Execute `VASP`_ using the prepared ``INCAR``, ``POSCAR``, ``KPOINTS`` and ``POTCAR`` from the previous step.

#. Inspect the resulting relaxed structure that will be housed in the created ``CONTCAR`` file. Its content should
   be something like::

     wurtzite-type SiC
	1.00000000000000
	  3.0920838394862775    0.0000000000000000    0.0000000000000000
	 -1.5460419197431388    2.6778231556249570   -0.0000000000000000
	 -0.0000000000000000   -0.0000000000000000    5.0729992295718249
	Si   C
	  2     2
     Direct
       0.3333333332999970  0.6666666667000030 -0.0000493014454436
       0.6666666667000030  0.3333333332999970  0.4999506985545563
       0.3333333332999970  0.6666666667000030  0.3758713014454431
       0.6666666667000030  0.3333333332999970  0.8758713014454431

       0.00000000E+00  0.00000000E+00  0.00000000E+00
       0.00000000E+00  0.00000000E+00  0.00000000E+00
       0.00000000E+00  0.00000000E+00  0.00000000E+00
       0.00000000E+00  0.00000000E+00  0.00000000E+00

4. Create two structures with decreased and increased volume
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

We now need to create two sets of VASP inputs.

#. Create a new ``POSCAR`` using the ``CONTCAR`` from the previous step, but modify it. The 2nd line of
   ``CONTCAR`` should be modified by replacing the scaling factor of 1.0 to a reduced volume. We will reduce it to
   99% of the original volume, meaning that the scaling factor should be set at :math:`0.99^{1/3}` = 0.9966554934125964.
   Save this to ``POSCAR_SMALLER``.

#. Create yet another ``POSCAR`` using the ``CONTCAR`` from the step 1 and 2.
   Again modify the second line to a larger scaling factor, which represent a 1% volume increase,
   meaning that the scaling factor should be set at :math:`1.01^{1/3}` = 1.0033222835420892.
   Save this to ``POSCAR_BIGGER``.

5. Perform the relaxation of modified volumes and 6. Wait for it to finish
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

We will now run `VASP`_ on the two modified structures with a decreased and increased volume, respectively
and inspect the resulting stress. One performing these relaxations, we want to keep the cell volume
static since we have modified this manually.

#. Modify the ``INCAR`` and set ``ISIF=4`` to ensure `VASP`_ does not relax the cell volume.

#. Execute VASP using the previous ``KPOINTS`` and ``POTCAR``, and newly generated ``INCAR`` and ``POSCAR_SMALLER``.
   Inspect the ``vasprun.xml`` file to find something like::

     <varray name="stress" >
      <v>      21.84700539      -0.00000000       0.00000000 </v>
      <v>       0.00000000      21.84700539      -0.00000000 </v>
      <v>       0.00000000       0.00000000      21.84572011 </v>
     </varray>

   for the stress on the last ionic step. Note that the volume of the cell is 41.59 :math:`\unicode{x212B}^3`. Keep these results at hand for later.

#. Execute VASP yet again using the previous ``INCAR``, ``KPOINTS`` and ``POTCAR``, but use now ``POSCAR_BIGGER``.
   Inspect again the new ``vasprun.xml`` file to find something like::

     <varray name="stress" >
      <v>     -20.80399703       0.00000000       0.00000000 </v>
      <v>       0.00000000     -20.80399703       0.00000000 </v>
      <v>       0.00000000       0.00000000     -20.80406790 </v>
     </varray>

   for the stress on the last ionic step. Note that the volume of the cell is 42.43 :math:`\unicode{x212B}^3`. Keep these results at hand for later.

7. Calculate the bulk modulus
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The bulk modulus can now be calculated from these results as::

  $ verdi shell
  In [1]: - ( ( (42.43 + 41.59) / 2  ) * ( (-20.80 - 21.85) / (42.43 - 41.59) ) ) / 10
  Out[1]: 213.3007738095248

We thus obtain the bulk modulus of ~213 GPa for this calculation.

3. Calculations using AiiDA
---------------------------

If there is any intention to perform bulk modulus calculations in a repeatedly and robust manner, the workflow above should be
define more formally. `AiiDA`_ comes into play with the concept of workflows.
Let us now try to perform the same calculation with assistance from `AiiDA`_.

We start by making the calling script, which you can get by running::

  $ wget https://raw.githubusercontent.com/aiida-vasp/aiida-vasp/master/tutorials/run_bulk_modulus.py

and would look something like this:

.. literalinclude:: ../../../tutorials/run_bulk_modulus.py
   :language: python

Notice a few things:

#. We set ``builder.clean_workdir = Bool(False)``. This makes sure that the folder used on the defined
   computer (usually remotely on the cluster in most cases) for the VASP calculations are not removed after a successful calculation
   (which is usually the default). Then one can after a ``VaspCalculation`` is finished execute ``verdi calcjob gotocomputer <PK>``,
   where ``<PK>`` is the id given to the ``VaspCalculation`` process you want to inspect. This takes you directly to the
   submission folder where the input and output files can be directly inspected. This setting is rather useful to debug and
   verify that your input and output files are as expected, regardless of the plugin.

#. When we use ``submit`` the process node is returned when we call ``group.add_nodes()`` passing the relevant node. This does however
   not mean that the node represents a finalized and successful calculation. This process will continue. Hence, we this need
   to explicitly check that the calculations are finished, here in a separate step using the node attribute ``is_terminated`` and
   later also checking that the node does not contain any exit code, which means the calculation was not successful by checking the
   node attribute ``is_finished_ok``. Only when all calculations are in a successful finished state can we proceed to calculate the
   bulk modulus.

#. This bulk modulus script assumes the AiiDA group named ``Bulk_modulus_SiC_test`` already exists. This group is created by::

     $ verdi group create "Bulk_modulus_SiC_test"

   And we can use::

       $ verdi group list

   to inspect the present groups.

#. We only rudimentary search the group to find the nodes which contains the results as outputs we need to calculate the
   bulk modulus. Hence, we can not run this script a second time as we would then pick up the old members
   of the group and print the bulk modulus of that. For some cases, this is good, meaning we reuse previous calculations so that
   we do not have to recalculate them, but for the sake of the example it can be a bit confusing. Hence, we check that the
   group ``Bulk_modulus_SiC_test`` is empty as the beginning of the ``main`` and stop if it is not. One can then, if for instance
   if rerunning the script, delete the group by first getting the id with ``verdi group list``, locating the ``<PK>`` and then
   issuing ``verdi group delete <PK>``, followed by a new ``verdi group create "Bulk_modulus_SiC_test"``. Then the script can be
   launched again.

Modify the ``run_bulk_modulus.py`` script to your setup. If in doubt, go back to previous tutorials. In essence, you need to look at the
code defined, the potential family you have uploaded potentials into and the resources to be used
where the code is installed. When that is done, execute the script::

  $ python run_bulk_modulus.py
  Relaxing starting cell: [[3.092, 0.0, 0.0], [-1.546, 2.6777505485015, 0.0], [0.0, 0.0, 5.073]]
  Waiting for relaxation of starting cell to be done.
  Waiting for relaxation of starting cell to be done.
  Waiting for relaxation of starting cell to be done.
  Waiting for relaxation of starting cell to be done.
  Waiting for relaxation of starting cell to be done.
  Waiting for relaxation of starting cell to be done.
  Waiting for relaxation of starting cell to be done.
  Waiting for relaxation of starting cell to be done.
  Waiting for relaxation of starting cell to be done.
  Waiting for relaxation of starting cell to be done.
  Relaxation positions of scaled cell: Cell([[3.081742345228316, 0.0, 0.0], [-1.540871172614158, 2.668867162801478, 0.0], [0.0, 0.0, 5.056032550657371]])
  Relaxation positions of scaled cell: Cell([[3.102356619252392, 0.0, 0.0], [-1.551178309626196, 2.686719647813093, 0.0], [0.0, 0.0, 5.0898531718508595]])
  Waiting for all remaining relaxations to be done.
  Waiting for all remaining relaxations to be done.
  Waiting for all remaining relaxations to be done.
  Waiting for all remaining relaxations to be done.
  Waiting for all remaining relaxations to be done.
  Waiting for all remaining relaxations to be done.
  Waiting for all remaining relaxations to be done.
  Waiting for all remaining relaxations to be done.
  Waiting for all remaining relaxations to be done.
  Waiting for all remaining relaxations to be done.
  Waiting for all remaining relaxations to be done.
  Bulk modulus: 213.2516947999474 GPa

Depending on your setup it will take some minutes before the bulk modulus is printed, assuming your submitted calculations does
not queue and start immediately. As we see, it is similar to the result we obtained using the manual method above.
Obviously, if we need to redo this calculation, either using different input parameters
and/or other input structures it would be little overhead involved and we could immediately benefit from it during our daily workflows.
However, we can do better.

4. Writing bulk modulus workchain
---------------------------------

In the previous section, the initial relaxation and the two volume restricted relaxations are performed
independently. The resulting nodes are just grouped and we execute a
calculation of the bulk modulus by fetching these nodes later. This
means the workflow defined in the ``main`` method is not very transparent and we have to introduce while
iterations and sleep the code for some time before continuing the iteration. This works, but is
not very elegant. In fact, we would also like our workflow to follow the steps of the flow of work we defined at the top,
also formally. In doing so, it makes it easier to both design new and understand existing workflows. Furthermore,
it makes reuse easier due to the modular design of it all.

The question is thus; can we improve the somewhat large script introduced in the previous
section? The answer is a definitive; yes. `AiiDA`_ has the
concept of a :ref:`workchains` which basically is a container for a
workflow, or a component in a workflow. It is supposed to be *the* modular concept in `AiiDA`_
that gives us enough freedom to dissect and combine workchains into basic and complex workflows.
It is apparent that we would like to write a workchain as modular
and reusable as possible, such that ultimately several :ref:`WorkChains` can
be cherry picked into a bigger composition of some kind of master
piece of a workchain to solve some given problem. And in doing so, keeping as close to the standard and
expectations given in the domain so that others can rapidly relate to and understand what the workflow
and its workchains are doing.

.. note::

   The additional :ref:`workchains` not bundled with `AiiDA-VASP`_ have to be exposed as a
   Python module to be seen from `AiiDA`_ daemon. It should be availble in ``$PYTHONPATH`` of the daemon. The easiest way
   to achieve this is usually to make a small plugin that you can use as a container for the installation and documentation
   of the workchain, or set of workchain. This can then be installed (with ``pip`` for instance) in
   an AiiDA environment where you want to utilize the workchain(s).

So, let us now try to preserve the workflow described at the top more formally in code as well.
The challenge will then be writing a suitable workchain that represents this workflow. It turns out that
the migration from the ``run_bulk_modulus.py`` script above is rather straightforward. In the process we will
also try to develop the workchain in a way which makes it easy to maintain, develop and share. In,
addition to just transferring the run script developed above into a more robust workflow we will at the same time
look into:

#. Using the `AiiDA plugin cookie cutter`_ to prepare the structure of the plugin, which will be a very simple
   one in our case. The workchain what will define our workflow which we will call from a suitable run script.

   .. note::
      `AiiDA`_ is designed to be very modular and by making our single workflow for the bulk modulus is on its own a
      plugin,we make sure it can also be installed easily and correctly by many users. The overhead of doing this
      is highly worth it and sounds more daunting than it is in practice. You will also benefit from this.

#. Writing documentation.

   .. note::
      Having at least basic documentation in place is very important, not only for others, but also for yourself.
      Everyone forget and you forget more quickly than you think. Also notice that it makes sense to be as explicit
      as possible when writing the workchain code. When doing this, it is much easier for new people to inspect and
      understand what the code actually does.

#. Enabling possibilities to share the plugin with others due to easy install.

Okey, let us now start. We have already a pretty good mental picture on what we need to calculate and what to extract from the
calculations. Buckle up, this is going to be a very compact introduction to how you can write plugins and workchains following
good practices.

1. Generate the plugin structure and basics
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

We will now start by utilizing the `AiiDA plugin cookie cutter`_. If you are not familiar with the cookie cutter
principle, it is a way to quickly get baseline structure in place for modular code pieces or extensions, which is
perfect for `AiiDA`_ plugins.

#. Make sure you have your `AiiDA`_ Python environment enabled.

#. Go to some location which will house the folder and `Git`_ repository of the bulk modulus plugin we will develop.

#. Install ``cookiecutter``::

     $ pip install cookiecutter

#. Run the ``cookiecutter`` and answer the questions accordingly::

     $ cookiecutter https://github.com/aiidateam/aiida-plugin-cutter.git
     plugin_name [aiida-diff]: aiida-vasp-bm
     module_name [aiida_vasp_bm]:
     short_description [AiiDA demo plugin that wraps the `diff` executable for computing the difference between two files.]: A workflow for AiiDA-VASP that calculates the bulk modulus.
     entry_point_prefix [vasp_bm]:
     github_user [aiidateam]: <yourgithubusername>
     repo_name [aiida-vasp-bm]:
     contact_email []: <yourcontactemail>
     version [0.1.0a0]:
     author [The AiiDA Team]: The AiiDA-VASP developer team
     aiida_min_version [2.0]: 2.2
     year [2022]: 2023
     Running black on aiida-vasp-bm
     All done! ✨ 🍰 ✨
     14 files left unchanged.

     ### IMPORTANT###
     Register your plugin NOW by making a pull request to the AiiDA Plugin Registry at

     https://aiidateam.github.io/aiida-registry/

   where ``<yourgithubusername>`` is your username on GitHub, assuming you eventually want to publish the workflow there. Also,
   a contact email ``<yourcontactemail>`` should be provided.

   Enter the resulting ``aiida-vasp-bm`` folder and inspect it. The first thing you should do is to read the ``README.md`` file to
   know what the  ``AiiDA plugin cookie cutter`_ provides and how to further configure it. Out of the box, the following is
   provided:

     - A license file. MIT is provided out of the box.

     - The ``pyproject.toml`` file. Contains information about the plugin, dependencies etc.
       This also enabled us to install the plugin with e.g.

     - A basic example calculation (``aiida_vasp_bm/calculations.py::DiffCalculation``) and
       a basic parser (``aiida_vasp_bm/parsers.py::DiffParser``). In addition a few command line examples are given in
       ``aiida_vasp_bm/cli.py``. None of this is relevant for us. We will instead reuse the ``VaspCalculation``
       in the `AiiDA-VASP`_ plugin, the `AiiDA-VASP`_ parsers and its command line infrastructure.

     - Basic tests. Testing your code is important, so a testing framework is also provided to test the example calculation
       and cli in the ``tests`` folder.

     - Documentation. A sample using Sphinx that already can be compiled is found in the ``docs`` folder. In there, the
       web documentation can be built by running ``make html``.

2. Adapting the generated plugin template
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Adapting the ``cookiecutter`` output from the previous step to our specific involves both removing, modifying existing files
and adding the file containing our workchain.

Since what we need is a very simple plugin, basically only containing the workchain housing our workflow, some light
documentation etc. we will now remove what we do not need from the output of the ``cookiecutter`` command in the
previous step. Please remove:

- The ``tests`` folder in the root directory. We would generally advice against not prioritizing testing, but as stated above, testing
  workchains which use external codes is rather complicated and is not the topic of this tutorial. Also remove everything inside
  the ``examples`` folder.

- The ``calculations.py``, ``cli.py``, ``parsers.py`` and ``helpers.py`` file and the ``data`` folder
  from the ``aiida_vasp_bm`` folder.

- The ``conftest.py`` from the root directory.

- Remove the ``testing`` section from the ``[project.optional-dependencies]`` section in ``pyproject.toml``. Also remove the
  ``project.entry-points."aiida.data``, ``project.entry-points."aiida.parsers``, ``project.entry-points."aiida.calculations`` and
  ``project.entry-points."aiida.cmdline.data`` sections. Furthermore, since the `AiiDA-VASP`_ plugin depends on
  `AiiDA`_ we can remove that and add the proper version of `AiiDA-VASP`_ to
  the ``dependencies`` section so that it now reads::

    dependencies = [
      "voluptuous",
      aiida-vasp
    ]

  and the ``project.optional-dependencies`` so that it reads::

    [project.optional-dependencies]
    testing = [
      "aiida-vasp[tests]"
    ]
    pre-commit = [
      "aiida-vasp[pre-commit]"
    ]
    docs = [
      "aiida-vasp[pre-commit]"
    ]

  By changing this, we make sure the dependencies of this plugin follows what is given in the latest `AiiDA-VASP`_ version.
  We also add what is called an entry point for our workchain::

    [project.entry-points."aiida.workflows"]
    "vasp.vasp_bm" = "aiida_vasp_bm.workchains.bulkmodulus:BulkModulusWorkChain"

  which makes it possible later to use the ``WorkflowFactory`` function of `AiiDA`_ to load our workchain
  in the calling script using the ``vasp.bm`` key.

  Furthermore, make sure ``requires-python = ">=3.8"`` and remove ``37`` under the ``testenv:py`` entry,
  which is also something we require in `AiiDA-VASP`_. The ``pyproject.toml`` should now look something like:

  .. rli:: https://raw.githubusercontent.com/aiida-vasp/aiida-vasp-bm/master/aiida_vasp_bm/pyproject.toml

We now need to add the actual workchain to the plugin.

First get the workchain file by running::

  $ wget https://raw.githubusercontent.com/aiida-vasp/aiida-vasp-bm/master/aiida_vasp_bm/workchains/bulkmodulus.py

and put this file in the ``aiida_vasp_bm`` folder. It should look something like:

.. rli:: https://raw.githubusercontent.com/aiida-vasp/aiida-vasp-bm/master/aiida_vasp_bm/workchains/bulkmodulus.py
   :language: python

The workflow is the same as what we implemented in ``run_bulk_modulus.py``, but here, we have increased the
verbosity and tried to be more explicit to comply with a typical implementation of a workchain in `AiiDA`_.
Before continuing, you should read up on `how you write workchains in AiiDA`_. The main functionality for our workchain,
which will calculate the bulk modulus is housed in the class ``BulkModulusWorkChain``, which inherits from
the general ``WorkChain`` class of `AiiDA`_. Consult the `how you write workchains in AiiDA`_ for more details of the structure of this,
but summarized very briefly here: (i) we define what is going to be supplied as inputs, what is going to be done, or executed,
and what is going to be an output in the ``define`` method of the ``BulkModulusWorkChain`` class. As you also noticed
from the ``run_bulk_modulus.py`` script, we eventually called the ``RelaxWorkChain`` class of the `AiiDA-VASP`_ plugin multiple
times. We also need to do this here. The ``spec.expose_inputs()`` configures the ``BulkModulusWorkChain`` to expect the
same inputs as defined in another workchain, in this case, the workchain we will call later, the ``RelaxWorkChain``. In
addition, we specify an output that is required for this particular workchain, the ``bulk_modulus``, which is required to be
a simple float datatype.

.. note::
   Again, remember that `AiiDA`_ have wrappers for the usual Python datatypes, including ``float``, represented as ``Float``.
   This is due to the fact that in `AiiDA`_, the datatypes often need to be passed around and they need to obey data provenance
   principles. As such they are defined as nodes which gives additional functionalities over the basic data types, such as
   the ability to store, lock and link them. Also, as you have maybe already noticed from previous tutorials, `AiiDA`_
   also supplies other dedicated data types. All of these can be used as input or output of workchains.

Finally, we define what is to be done in the ``outline`` section. This is sequential and refer to functions defined on
this class further down in the file. It contains the following functions:

#. ``initialize`` Initialize whatever needs initializing.

#. ``run_relax`` Run a relaxation calculation to fully relax the initial crystal structure.

#. ``create_two_structures`` Create two structures with a +/- 1% change
   in the volume from the relaxed structure obtained at the previous step.

#. ``run_two_volumes`` Run two relaxation calculations to relax the shape and positions of the
   structures created at the previous step at fixed volume.

#. ``calc_bulk_modulus`` Compute the bulk modulus using the
   formula :math:`K \simeq -V_0 \frac{\Delta P}{\Delta V}`, where the
   pressure :math:`P \equiv \mathrm{Tr}(\sigma)/3` follows the VASP
   convention where :math:`\sigma` is the stress tensor.

You can inspect the details of the ``BulkModulusWorkChain`` otherwise and inspect that it largely
takes the same inputs and gives the same outputs as we did in the manual and the ``run_bulk_modulus.py`` way above.
In the end of this section, we now would like to raise a few notes:

.. note::
   Since the workchain definition is sequential, the ``calc_bulk_modulus`` function will not be executed until the previous steps
   are completed.

.. note::
   Data created inside a workchain should be registered so that data provenance can be kept. This can for simple functions
   be performed using the `calcfunction concept`_ and decorator in `AiiDA`_. For this purpose, the
   strained crystal structures are created in ``get_strained_structure`` and the bulk modulus is calculated in
   ``calculate_bulk_modulus``, and its method is decorated by ``@calcfunction``. These functions are called from
   the ``BulkModulusWorkChain`` when needed. The calcfunction concept is a very compact and slim way to pass tracked
   input and output to and from a function, thus honoring data provenance. For more complex functions, it makes sense to
   instead consider porting the function to a dedicated workchain.

.. note::
   We have utilize ``ToContext`` and ``to_context``, which ensures that we can make sure the workflow does not continue until
   what is set in these are done. Please consult the `AiiDA`_ documentation on the `to context`_ topic to learn more.

3. Creating the launch script
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The launch script will define our input structure and launch the ``BulkModulusWorkChain`` workchain which
will contain the ``bulk_modulus`` output node that we can inspect when it is finalized. Obtain the launch script::

  $ wget https://raw.githubusercontent.com/aiida-vasp/aiida-vasp-bm/master/aiida_vasp_bm/examples/run_vasp_bm.py

It should look something like:

.. rli:: https://raw.githubusercontent.com/aiida-vasp/aiida-vasp-bm/master/aiida_vasp_bm/examples/run_vasp_bm.py
   :language: python

Place this call script in the ``examples`` directory of the plugin.

We also have to install the plugin we have created, such that the `AiiDA`_ daemon can pick it up. It is convenient
to install it in what we call editable mode. Doing this, makes it possible to utilize the framework of installation
and the fact that other packages like `AiiDA`_ becomes aware of our new plugin and workchain, while maintaining the
convenience of being able to quickly modify the code and relaunch. Assuming you are not in the root folder of the
``aiida-vasp-bm`` plugin (the folder containing the ``aiida_vasp_bm`` and ``pyproject.toml`` file), execute::

  $ pip install -e .

which will install the plugin in editable mode.

.. warning::
   Please remember to restart the daemon if you change any code that the ``AiiDA`_ daemon should have access to. This typically
   goes for all editable (the ones installed with the ``-e`` option for ``pip``) installations of Python packages. This
   is very easy to miss and during heavy development it sometimes make sense to restart the daemon in the call script
   to not forget. Beware that if you forget this, it can cause a significant time drain on debugging something that already
   might be working just fine.

Again, tailor the launch script to your needs as we now have done a few times. Launch it. Note the id ``<PK>>>`` printed. Monitor its
progress with e.g. ``verdi process report <PK>``. When done, fetch its ``bulk_modulus`` output id ``<<PK'>>`` from
``verdi process show <PK>`` and inspect it::

  $ verdi node attributes <PK'>
  PK: 6700
  {
      "value": 213.25056684995
  }

or check it in Python with e.g.::

  $ verdi shell
  In [1]: node = load_node(<PK>)

  In [2]: node.outputs.bulk_modulus.value
  Out[3]: 213.25056684995

.. note::
   The basic data types, meaning, the ones that are not listed under ``verdi data``. In this case, ``Float`` can
   be inspected with ``verdi node attributes`` from the command line. From Python the entry point for inspecting this
   is the same.

As we see, it yields the same bulk modulus as before. We have now made a new plugin and a workflow that we can
easily share with others and perform bulk modulus in a consistent manner on many different input structures. We can also
now include the ``BulkModulusWorkChain`` in other workflows if we need to calculate the bulk modulus as a sub component of
that.

4. Sharing your workflow
^^^^^^^^^^^^^^^^^^^^^^^^

In order to spread the uptake and also open for improvements and suggestions from the field, we suggest that you consider to
make your workflow publicly accessible. There are mainly two ways:

- Contribute to the `AiiDA-VASP`_ plugin such that your developed workflow eventually is part of the `AiiDA-VASP`_ distrubution.

- Publish it on GitHub or other standardized frameworks that supports publications of code using `Git`_.

For the first option, please open an issue on our `AiiDA-VASP`_ repository and explain what you have done and then we take it from there.
In order to show what you have done, showing the code is usually a good approach anyway, so the last option is still highly useful.

In order to upload your plugin to GitHub, please first make sure you have an account and that your GitHub account is reflected in the
``<yourgithubusername>>`` in the ``pyproject.toml``. Also, make sure to update the author information in the README.md, ``pyproject.toml``
and the documentation folder ``docs``. Please make sure everyone that should be credited is in fact credited. Then create a new
project with the name ``<projectname>`` on GitHub and perform the following steps in the plugin root folder (here ``aiida-vasp-bm`` for this example, but it is
rather unlikely you want to upload this particular plugin) folder::

  $ git init
  $ git add -A
  $ git commit -m "Initial commit"
  $ git branch -M main
  $ git remote add origin git@github.com:<yourgithubusername>/<projectname>.git
  $ git push -u origin main

Your repository, plugin and workflow should now, if you marked it with ``Public`` during the project creation, be available to everyone on
``https://github.com/<yourgithubusername>/aiida-vasp-bm`` (for this particular example).

.. _Git: https://git-scm.com/
.. _AiiDA: https://www.aiida.net
.. _how you write workchains in AiiDA: https://aiida.readthedocs.io/projects/aiida-core/en/latest/howto/write_workflows.html
.. _calcfunction concept: https://aiida.readthedocs.io/projects/aiida-core/en/latest/topics/calculations/concepts.html#topics-calculations-concepts-calcfunctions
.. _VASP: https://www.vasp.at
.. _AiiDA-VASP: https://github.com/aiida-vasp/aiida-vasp
.. _AiiDA plugin cookie cutter: https://github.com/aiidateam/aiida-plugin-cutter
.. _to context: https://aiida.readthedocs.io/projects/aiida-core/en/latest/topics/workflows/usage.html#to-context
.. _VASP lectures: https://www.vasp.at/wiki/index.php/Lectures_and_presentations
.. _VASP tutorials: https://www.vasp.at/wiki/index.php/Category:Tutorials
.. _VASP howtos: https://www.vasp.at/wiki/index.php/Category:Howto
.. _VASP tutorials using notebooks: https://www.vasp.at/tutorials/latest/
.. _VASP videos: https://www.youtube.com/channel/UCBATkNZ7pkAXU9tx7GVhlaw
