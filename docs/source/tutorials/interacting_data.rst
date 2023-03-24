.. _interacting_with_data:

====================
6. Where is my data?
====================

So far, we have followed tutorials, step by step and everything might seem very clear. However,
after doing your first `VASP`_ calculation that is not in a tutorial form, you might have
a few questions. One of them is typically, where is my data and how can I interact with it?
In this tutorial we will try to shed some light on this and how you can approach this from `AiiDA`_
and the

As you might have noticed we have now populated the database with the silicon structure
multiple times. Let us try instead to load some of the structures that already are present
to save coding, reuse previous results and save data storage.

In order to illustrate this we will need a finished ``VaspWorkChain`` and we will use the one
from the :ref:`previous tutorial<tutorial_fcc_si_dos>`. You can chose any ``VaspWorkChain`` you
want from your side, but notice that the input and outputs might be slightly different. If you are
confused, please complete :ref:`previous tutorial<tutorial_fcc_si_dos>` and continue to use the
results from that as we do here.

One of the first questions new plugin, but experienced `VASP`_ users ask is, where can I find my original files?
Notice first, that the default configuration of the plugin removes the downloaded files after successful
parsing unless it is explicitly told not to do so, or to keep specific files. Let us first investigate
how we can tailor this. It is rather simple. First, we always try to retrieve what is stored in the
``_ALWAYS_RETRIEVE_LIST`` in ``calcs/vasp.py``, which currently is the ``CONTCAR``, ``OUTCAR``, ``vasprun.xml``, ``vasp_output``
, any ``wannier90*`` file. The default setting ``settings['ALWAYS_STORE'] = True`` make sure these files will always be
kept after a successful parsing step, which also results in a successful ``VaspCalculation``. One can specify
``settings['ALWAYS_STORE'] = False`` to delete these files after a successful parsing step if one only want to
keep the parsed data, which are typically stored as `AiiDA`_ data nodes on the output of the process node. In
addition, it is possible to fine tune this. For instance, one might want to make sure to also download additional
files, like the ``CHGCAR``. This can be done by specifying ``settings['ADDITIONAL_RETRIEVE_LIST'] = ['CHGCAR']`` or
``settings['ADDITIONAL_RETRIEVE_TEMPORARY_LIST'] = ['CHGCAR']``, where the additional entries in the former (latter)
will be stored (deleted) by default after parsing, regardless of the values of ``ALWAYS_STORE``.

With this in mind, let us interact with data and see how this is manifested for our chosen example.

#. Let us first inspect the outputs of our previous ``VaspWorkChain``::

     $ verdi process show 2431
     Property     Value
     -----------  ------------------------------------
     type         VaspWorkChain
     state        Finished [0]
     pk           2431
     uuid         40ce7bd6-cd38-405e-951e-c56251a0cf1b
     label
     description
     ctime        2022-12-22 11:17:25.967623+01:00
     mtime        2022-12-22 11:19:39.816274+01:00

     Inputs               PK  Type
     -----------------  ----  -------------
     clean_workdir      2430  Bool
     code                818  InstalledCode
     kpoints            2422  KpointsData
     max_iterations     2429  Int
     options            2426  Dict
     parameters         2423  Dict
     potential_family   2424  Str
     potential_mapping  2425  Dict
     settings           2427  Dict
     structure          1529  StructureData
     verbose            2428  Bool

     Outputs          PK  Type
     -------------  ----  ----------
     dos            2436  ArrayData
     misc           2437  Dict
     remote_folder  2434  RemoteData
     retrieved      2435  FolderData

     Called          PK  Type
     ------------  ----  ---------------
     iteration_01  2433  VaspCalculation

     Log messages
     ---------------------------------------------
     There are 3 log messages for this calculation
     Run 'verdi process report 2431' to see them

   Pay particular notice to the outputs, especially the ``retrieved``, which is of a ``FolderData`` type.
   Let us now inspect it.

#. We can use the ``verdi node repo`` command. Let us first check what is in the folder::

     $ verdi node repo ls 2435
     CONTCAR
     DOSCAR
     EIGENVAL
     OUTCAR
     _scheduler-stderr.txt
     _scheduler-stdout.txt
     vasp_output
     vasprun.xml

   As we can see, this is the default files listed in ``_ALWAYS_RETREIVE_LIST``. In addition, there are
   the scheduler standard stream files, which is added by `AiiDA`_.

#. Let us have a look at the content of for instance ``CONTCAR``::

     $ verdi node repo cat 2435 CONTCAR
     # Compound: Si. Old comment: silicon_at_
        1.0000000000000000
	  1.9500000000000000    1.9500000000000000    0.0000000000000000
	  0.0000000000000000    1.9500000000000000    1.9500000000000000
	  1.9500000000000000    0.0000000000000000    1.9500000000000000
	Si
	  1
     Direct
       0.0000000000000000  0.0000000000000000  0.0000000000000000

       0.00000000E+00  0.00000000E+00  0.00000000E+00

   If you want, this can be piped to a file and displayed using regular tools::

     $ verdi node repo cat 2435 CONTCAR > /tmp/contcar
     $ more /tmp/contcar
     # Compound: Si. Old comment: silicon_at_
        1.0000000000000000
	  1.9500000000000000    1.9500000000000000    0.0000000000000000
	  0.0000000000000000    1.9500000000000000    1.9500000000000000
	  1.9500000000000000    0.0000000000000000    1.9500000000000000
	Si
	  1
     Direct
       0.0000000000000000  0.0000000000000000  0.0000000000000000

       0.00000000E+00  0.00000000E+00  0.00000000E+00

   So getting to your files requires a bit more typing than what seems comparable to working
   with folders and files in the traditional way, but this is only relevant for simple one off examples.
   Once, the workflow becomes more involved and the nesting of folders much more complicated, the
   typing involved quickly becomes more compact using `AiiDA`_, but of course, the main benefits is
   in everything that comes along with it.

#. Inspecting data, or working with it in general programmatic way is also very easy using the ``verdi shell``, which
   gives you access to an `IPython`_ instance where most of the needed `AiiDA`_ functionality is loaded for you.
   Launch the ``verdi shell``::

     $ verdi shell

   Then we load the node::

     In [1]: node = load_node(2435)

   And inspect the objects residing in the ``retrieved`` folder::

     In [2]: node.base.repository.list_object_names()
     Out[2]:
     ['CONTCAR',
     'OUTCAR',
     '_scheduler-stderr.txt',
     '_scheduler-stdout.txt',
     'vasp_output',
     'vasprun.xml']

   As we can see, as before, this is the default files listed in ``_ALWAYS_RETREIVE_LIST``, in addition to the
   scheduler files.

   .. note::
      For most commands, tab completion is available so you can write ``node.`` and then tab
      complete it to check what methods (with parenthesis) or attributes (no parenthesis) are available on the node.
      Notice however, that most of the useful methods and attributes are now placed into sub-namespaces under ``base``,
      see documentation on `namespace change`_ for more details.

   We can now inspect the content of these files::

     In [3]: node.base.repository.get_object_content('CONTCAR')
     Out[3]: '# Compound: Si. Old comment: silicon_at_\n   1.0000000000000000     \n     1.9500000000000000    1.9500000000000000    0.0000000000000000\n     0.0000000000000000    1.9500000000000000    1.9500000000000000\n     1.9500000000000000    0.0000000000000000    1.9500000000000000\n   Si\n     1\nDirect\n  0.0000000000000000  0.0000000000000000  0.0000000000000000\n\n  0.00000000E+00  0.00000000E+00  0.00000000E+00\n'

   And the content is available as a string. We can also of course dump this to a file::

     In [4]: with open('/tmp/contcar', 'w') as fo:
     ...:     fo.write(node.base.repository.get_object_content('CONTCAR'))
     ...:

   exit the ``verdi shell`` by typing ``exit`` and issue::

     $ more /tmp/contcar

   and there you again see the ``CONTCAR`` from the `VASP`_ calculation.

.. _AiiDA: https://www.aiida.net
.. _density of states for FCC Si: https://www.vasp.at/wiki/index.php/Fcc_Si_DOS
.. _VASP: https://www.vasp.at
.. _AiiDA-VASP: https://github.com/aiida-vasp/aiida-vasp
.. _IPython: https://ipython.org/
.. _namespace change: https://aiida.readthedocs.io/projects/aiida-core/en/latest/reference/_changelog.html#node-namespace-restructuring
