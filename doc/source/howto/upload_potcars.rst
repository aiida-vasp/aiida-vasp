How to import POTCAR files
==========================

Why do I need to import POTCARS?
--------------------------------

AiiDA does more than prepare calculations and send them to a cluster. The main focus of AiiDA lies on tracking data provenance. The same goes for AiiDA-VASP, importing (or "uploading" the POTCAR files into your personal AiiDA database helps with that in the following ways:

   * AiiDA-VASP stores a unique hash for each file. This can help when different potentials have very similar looking headers.
   * POTCAR files uploaded to the DB can not be modified accidentally, thus it is recorded unambiguously, which file was used for which run.
   * Storing the file's contents rather than a link prevents accidentially breaking the link by moving the file away (or renaming it).

How to import a set of POTCAR files?
------------------------------------

::

   $ verdi data vasp-potcar uploadfamily --path=<path> --name=<name> --description=<desc>

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

``<name>`` is the label you will use to access POTCARs from this set or to specify which "family" of potentials you want to use in a particular VASP run. The meaning of ``<description>`` is self-evident.

Custom sets can simply be arranged in a matching folder structure and then imported using the same command.

How to hand-pick a set of POTCAR files?
---------------------------------------

For this purpose, we can use that the ``uploadfamily`` command by default adds any POTCAR files not yet uploaded to the family of the given ``name``, for example::

   $ verdi data vasp-potcar uploadfamily --path=path/to/Ac --name="PBE_custom" --description="A custom set"
   $ verdi data vasp-potcar uploadfamily --path=other/path/to/Ag --name="PBE_custom"

Note, that the description does not have to be given if the family already exists.

Due to the recursive nature of the search, this also works for combining several small sets of POTCARS in a few commands, without having to arrange them in a different way first.

How to Check what POTCAR families I already have?
-------------------------------------------------

::

   $ verdi data vasp-potcar listfamilies

How to Access uploaded POTCAR files?
------------------------------------

The data structure used to find and retrieve POTCAR files is called ``aiida_vasp.data.potcar.PotcarData`` and can be accessed through AiiDA's factory as ``DataFactory('vasp.potcar')``. This class provides shortcuts for simple or frequent cases, for complex queries, please refer to the AiiDA documentation on querying the provenance graph.

Find by exact properties
^^^^^^^^^^^^^^^^^^^^^^^^

More advanced searches, like for ranges of properties etc can be done using the ``QueryBuilder`` tool, which is part of AiiDA and documented there.

Use::

   PotcarData.find(property=value, property2=value2)

, it returns a list of all stored ``PotcarData`` instances fullfilling the criteria. Some important supported properties are:

   * ``md5``, an MD5 hash of the file contents
   * ``title``, Title of the POTCAR file
   * ``element``, The chemical element described by this potential
   * ``full_name``, The name of the containing folder from which it was uploaded. This is used to specify a potential inside a family. Example: ``Zn_sv_GW``
   * ``original_file_name``, The filename (+ last three directories) from which it was uploaded (May help identifying exactly where it came from).

Find by a list of elements
^^^^^^^^^^^^^^^^^^^^^^^^^^

To find one POTCAR for each in a list of element names, all from the same family::

   mapping = {
      'Ac': 'Ac',
      'Ag': 'Ag_GW'  # or 'Ag_pv', 'Ag_sv_GW', ...
   }
   potcars_for_elements = PotcarData.get_potcars_dict(
      elements=['Ac', 'Ag', ..], <family_name>, mapping=mapping)

The ``mapping`` dictionary is required to decide which of the variants the set provides for each element should be chosen. The mapping can also conveniently be stored in a ``ParameterData`` node for reuse.

How to pass POTCAR potentials to a VASP run?
--------------------------------------------

For a single VASP run, you might use the ``vasp.base`` workflow, which takes the family as a database-storable string and a dictionary mapping elements to a particular variant for that element::

   from aiida.orm import WorkflowFactory
   from aiida.common.extendeddicts import AttributeDict
   from aiida.work.db_types import Str

   from aiida_vasp.utils.aiida_utils import get_data_node

   inputs = AttributeDict()
   inputs.potcar_famliy = Str('<name>')
   inputs.potcar_mapping = get_data_node('parameter', dict={'In': 'In_d', 'As': 'As'})

The ``VaspBaseWf`` takes care of finding the right files and concatenating them for you.

For a more complex workflow, the process may be different, it may for example use heuristics to find a default potential for you. Refer to the documentation of that specific workflow in that case.
