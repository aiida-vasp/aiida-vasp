.. _potentials:

Potentials
==========
As you already know, `VASP`_ relies on potentials represented by the POTCAR files.

AiiDA-VASP takes care of managing your POTCAR files, but you need to obtain them separately and make them available to AiiDA-VASP. These are usually supplied by the `VASP`_ team and is part of the license. You should have recieved a folder (``tar`` archive) containing multiple subfolders (``tar`` archives), each representing a set of POTCAR files intended to be used together. AiiDA-VASP allows you to import only the sets (or even individual potentials) you require, and keep them grouped in potential families. The import process uploads these potentials to you working AiiDA database.


Why do I need to import POTCAR files?
-------------------------------------

`AiiDA`_ does more than prepare calculations and send them to a cluster. The main focus of `AiiDA`_ lies on tracking data provenance. The same goes for AiiDA-VASP. Importing the POTCAR files into your working `AiiDA`_ database yields some advantages:

   * AiiDA-VASP stores a unique hash for each file. This can help users navigate when different potentials have very similar looking headers, but do in fact contain a different potential.
   * POTCAR files uploaded to the database cannot be modified accidentally, thus it is recorded unambiguously, which file was used for which execution of each run.
   * Storing the file's contents rather than a link prevents accidentially breaking the link by moving the file away (or renaming it).

How to import a set of POTCAR files?
------------------------------------

The command line tools for these tasks are written as plugins to `AiiDA`_, and can be called through the `AiiDA`_ command ``verdi``::

   $ (aiida-vasp) verdi data vasp-potcar --help
   Usage: verdi data vasp-potcar [OPTIONS] COMMAND [ARGS]...

      Top level command for handling VASP POTCAR files.

   Options:
     -h, --help  Show this message and exit.

   Commands:
     exportfamily  Export a POTCAR family into a compressed tar...
     listfamilies  List available families of VASP potcar files.
     uploadfamily  Upload a family of VASP potcar files.



::

   $ verdi data vasp-potcar uploadfamily --path=<path> --name=<potential_family> --description=<desc>

Where ``<path>`` is the path to the folder or tar archive containing the POTCAR set. The command expects the folder or archive to look like::

   <path>/
   |
   +- Ac/
   |  +- POTCAR
   |  +- ...
   |
   +- Ag/
   |  +- POTCAR
   |  +- ...
   ...

If it encounters anything different, it will recursively search the given path for subpaths matching this structure and import all the POTCAR files found in that way.

``<potential_family>`` is the label you will use to access the potentials from this set or to specify which potentials you want to use in a particular `VASP`_ run. The meaning of ``<description>`` is self-explanatory.

Custom sets can simply be arranged in a matching folder structure and then imported using the same command.

Uploading a set of potentials
-----------------------------

For this purpose, we can use that the ``uploadfamily`` command by default adds any POTCAR files not yet uploaded to the family of the given ``name``, for example::

   $ verdi data vasp-potcar uploadfamily --path=path/to/Ac --name="PBE_custom" --description="A custom set"
   $ verdi data vasp-potcar uploadfamily --path=other/path/to/Ag --name="PBE_custom"

Note, that the description does not have to be given if the family already exists.

Due to the recursive nature of the search, this also works for combining several small sets of POTCARs in a few commands, without having to arrange them in a different way first.

How to check what potential families are present in the database?
-----------------------------------------------------------------

::

   $ verdi data vasp-potcar listfamilies

How to access uploaded potentials and search?
---------------------------------------------

The data structure used to find and retrieve potentials is called :py:class:`aiida_vasp.data.potcar.PotcarData` and can be accessed through AiiDA's data factory as ``DataFactory('vasp.potcar')``. This class provides shortcuts for simple or frequent cases, for complex queries, please refer to the `AiiDA documentation`_ on querying the provenance graph.

Find potentials by properties
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

More advanced searches, like for ranges of properties etc can be done using the :py:class:`QueryBuilder<aiida.orm.querybuilder.QueryBuilder>` tool, which is part of `AiiDA`_ and documented there.

Use::

   PotcarData.find(<property>=<value>, <property2>=<value2>, ...)

which returns a list of all stored :py:class:`PotcarData<aiida_vasp.data.potcar.PotcarData>` instances fullfilling the criteria. Some important supported ``<property>`` entries are:

   * ``sha512`` - An SHA512 hash of the file contents
   * ``title`` - Title of the potential, typically the title of the POTCAR
   * ``element`` - The chemical element described by this potential
   * ``full_name`` - The name of the containing folder from which it was uploaded. This is used to specify a potential inside a family. Example: ``Zn_sv_GW``
   * ``original_file_name`` - The filename (+ last three directories) from which it was uploaded (May help identifying exactly where it came from).

and for each you supply the ``<value>`` which is relevant for you given search.

Find potentials by a list of elements
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To find one potential for each element in a list of element names, all from the same family::

   mapping = {
      'Ac': 'Ac',
      'Ag': 'Ag_GW'  # or 'Ag_pv', 'Ag_sv_GW', ...
   }
   potcars_for_elements = PotcarData.get_potcars_dict(
      elements=['Ac', 'Ag', ..], <potential_family>, mapping=mapping)

The ``mapping`` dictionary is required to decide which of the variants should be chosen for each element. The mapping can also conveniently be stored in a ``Dict`` node for reuse. The potential family is specified with ``<potential_family>``.

How to pass potentials to a VASP calculation?
---------------------------------------------

For a single `VASP`_ calculation run, you should at the very minimum use the :ref:`vasp_workchain` (although we recommend to use the :ref:`converge_workchain` as the standard entry point), which takes the family as a database-storable string and a dictionary mapping elements to a particular variant for that element::

   from aiida.plugins import DataFactory
   from aiida.common.extendeddicts import AttributeDict
   from aiida.orm import Str

   inputs = AttributeDict()
   inputs.potential_family = Str('<potential_family>')
   inputs.potential_mapping = DataFactory('dict')(dict={'In': 'In_d', 'As': 'As'})

The :ref:`vasp_workchain` takes care of finding the right files and concatenating them for you.

For a more complex workflow, the process may be different, it may for example use heuristics to find a default potential for you.

.. _AiiDA: https://www.aiida.net
.. _VASP: https://www.vasp.at
.. _AiiDA documentation: http://aiida-core.readthedocs.io/en/latest/
