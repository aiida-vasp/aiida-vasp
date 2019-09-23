Tips to write scripts using AiiDA-VASP
=======================================

Python interactive shell and ipython
------------------------------------

Python has the interactive shell mode. Using this, python script can
be run one line by one line interactively. This mode is invoked
usually by::

   % python
   Python 3.7.3 | packaged by conda-forge | (default, Jul  1 2019, 21:52:21)
   [GCC 7.3.0] :: Anaconda, Inc. on linux
   Type "help", "copyright", "credits" or "license" for more information.
   >>> print("This is the interactive mode.")
   This is the interactive mode.
   >>>

Instead of using python in this way, the interactive python with
enhanced features is used by ``ipython``::

   % ipython
   Python 3.7.3 | packaged by conda-forge | (default, Jul  1 2019, 21:52:21)
   Type "copyright", "credits" or "license" for more information.

   IPython 5.8.0 -- An enhanced Interactive Python.
   ?         -> Introduction and overview of IPython's features.
   %quickref -> Quick reference.
   help      -> Python's own help system.
   object?   -> Details about 'object', use 'object??' for extra details.

   In [1]: print("This is the interactive mode.")
   This is the interactive mode.

   In [2]:

In the following, it is assumed users use ``ipython`` instead of
normal ``python`` for using the interactive shell mode. The details of
``ipython`` is found at https://ipython.readthedocs.io/en/stable/. But
as explained below, ``ipython`` is used via ``verdi shell`` command
to interact with AiiDA.


``verdi shell``
---------------

If AiiDA core is properly installed and set-up, we are ready to use
``verdi shell`` command. This invokes ``ipython`` with loading AiiDA
profile and some useful classes and methods. So normally when we want
to use AiiDA from the interactive shell mode, ``verdi shell`` command
is used instead of ``ipython``.

::

   % verdi shell
   Python 3.7.3 | packaged by conda-forge | (default, Jul  1 2019, 21:52:21)
   Type "copyright", "credits" or "license" for more information.

   IPython 5.8.0 -- An enhanced Interactive Python.
   ?         -> Introduction and overview of IPython's features.
   %quickref -> Quick reference.
   help      -> Python's own help system.
   object?   -> Details about 'object', use 'object??' for extra details.

   In [1]: Str("PBE.54")
   Out[1]: <Str: uuid: 9f2e8b02-c3e6-491e-ba0a-b3fcd529a35e (unstored) value: PBE.54>

   In [2]:

More details on ``verdi shell`` is found `here
<https://aiida.readthedocs.io/projects/aiida-core/en/latest/working_with_aiida/scripting.html#verdi-shell>`_.

Python Docstring
----------------

In python, method, class, or module can have their Docstrings. The
Docstring explains what it is and we can read it easily by the ``help``
method in the interactive mode. By

::

   In [1]: help(Str)

a command line viewer is opened to show the Docstring.


Tab completion
--------------

In the interactive mode invoked by ``verdi shell``, tab completion is
extremely helpful in writing an AiiDA script to show the like of
attributes of classes and instances. Therefore we must set up
the verdi tab-completion following `this documentation
<https://aiida.readthedocs.io/projects/aiida-core/en/latest/install/configuration.html#verdi-tab-completion>`_.

::

   In [1]: val = Float(0.1)

   In [2]: val.          <-(now hit tab key!)
                 val.add_comment               val.backend
                 val.add_incoming              val.backend_entity
                 val.attributes                val.class_node_type           >
                 val.attributes_items          val.clear_attributes
                 val.attributes_keys           val.clear_extras

``DataFactory`` and ``WorkflowFactory``
----------------------------------------

These are methods that return classes. These method names have
"Factory" in them. This is from the factory method pattern, one of the
popular design patterns in object oriented programming.

In AiiDA, these take one parameter that is a string of
entry_point. For examples,

::

   Dict = DataFactory('dict')
   StructureData = DataFactory('structure')
   KpointsData = DataFactory('array.kpoints')
   Workflow = WorkflowFactory('vasp.vasp')
   Workflow = WorkflowFactory('vasp.relax')

Again, what is returned is a python class. Therefore to use, e.g.,
AiiDA Dict, we have to instantiate it by::

   mydict = Dict(dict={'foo': 'bar'})



``ProcessBuilder``
------------------

For WorkChain and CalcJob, we can take so called ProcessBuilder,
which is easily done by ``get_builder()``::

   MyBuilder = MyWorkflow.get_builder()

Below, how to use the ProcessBuilder is explained shortly. More
details are found at `AiiDA documentation
<https://aiida-core.readthedocs.io/en/latest/working/processes.html#working-processes-builder>`_. On
AiiDA process, it is nice to read `this
<https://aiida-core.readthedocs.io/en/latest/concepts/processes.html>`_
and `this
<https://aiida-core.readthedocs.io/en/latest/working/processes.html>`_
in the AiiDA documentation.

There are two ways to submit a process to AiiDA daemon. They are like
either

::

   from aiida.engine import submit
   submit(MyWorkchain, **inputs)

or

::

   from aiida.engine import submit
   submit(MyBuilder)

``inputs`` is a python dictionary containing parameters of the
process. These parameters are stored in ``MyBuilder`` as the
attributes, i.e.,

::

   MyBuilder.label = "My label"

instead of writing ``inputs['label'] = "My label"``. The advantage of
use of ProcessBuilder is that we can use tab completion on the
interactive mode.
