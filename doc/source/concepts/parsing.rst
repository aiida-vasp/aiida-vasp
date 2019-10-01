.. _parsing:

=======
Parsing
=======
AiiDA-VASP provides flexible parsing of `VASP`_ output files to store data in the `AiiDA`_ database and repository.

The quantities that can be parsed are now fully customisable. The user interface for configuring the parsing settings takes place in the ``settings['parser_settings']`` dictionary entry. The defualt ``parser_settings`` is presently:

.. literalinclude:: aiida_vasp.parsers.vasp
   :start-after: defaults
   :lines: 42

There are four ways to interact and set the parser properties.

#. Using a boolean::

     settings['parser_settings'] = {'add_<node_name>': True}

  This will set the name of the output to ``<node_name>`` as defined
  in ``aiida_vasp.parsers.settings.NODES``.

(2) list
^^^^^^^^

Which quantities are used for a node can be adjusted by setting::

  'add_<node_name>': ['<quantity1>', '<quantity2>', ... ]'

Where the list contains strings with quantity names that should be used for this node. This will work only for predefined nodes.

(3) dict
^^^^^^^^

The third option is to define completely custom nodes by setting e.g.,

::

  'add_<node_name>': { 'type': '<node_type>', 'quantities': ['<quantity1>', '<quantity2>', ...], 'link_name': '<link_name>'}

This will define a completely new custom node. The ``link_name`` is optional and will be set to the ``<node_name>`` if it is not present.

(4) Interactive way
^^^^^^^^^^^^^^^^^^^

These above three ways to interact with nodes were all mediated by the calculations settings. There is also an interactive way, if you have an instance of ``VaspParser``::

  VaspParser.add_custom_node(node_name, node_definition)

where the format for the ``node_definition`` is as in the previous example with the custom nodes.

NodeComposer
------------

A ``NodeComposer`` has been added to compose output nodes based on a given set of quantities. It can be initialised in two different ways::

  composer = NodeComposer(vasp_parser=...)
  composer = NodeComposer(file_parsers=[...])

where either a ``VaspParser`` object or a list of ``FileParser`` objects are required for the ``NodeComposer`` to get the required quantities for a node. Here is an example for the interface for composing a node from the VasprunParser test::

  composer = NodeComposer(file_parsers=[vasprun_parser])
  data_obj = composer.compose('array.kpoints', quantities=['kpoints'])

For most node types there are some default quantities defined in ``aiida_vasp.parsers.node_composer.NODES_TYPES``, so that the above could be shortened to

::

  data_obj = composer.compose('array.kpoints')


Short-cut properties
--------------------

Since now none of the quantities will actually return an Aiida node anymore, there are some built-in short-cut properties in order to save some typing and also to avoid using the ``get_quantity`` interface directly. An example would be::

  parser = PoscarParser(file_path=path)
  result = parser.structure

where the whole composer part is dealt with internally::

  @property
  def structure(self):
      if self._structure is None:
          composer = NodeComposer(file_parsers=[self])
          self._structure = composer.compose('structure', quantities=['poscar-structure'])
      return self._structure

At the moment existing short-cut properties are:

- ``PoscarParser.structure``
- ``KpParser.kpoints``
- ``OutcarParser.parameter``
- ``DoscarParser.dos``
- ``WavecarParser.wavecar``
- ``ChgcarParser.chgcar``
- ``IncarParser.incar`` (This is an exception in so far that it actually just returns a dictionary)

Prominently missing from this list are ``bands`` from the ``EigenvalParser`` and everything from the ``VasprunParser``. For the former the reason is that ``bands`` require quantities from other FileParsers and for the latter the namespace is already very crowded.

Quantity definitions and alternatives
-------------------------------------

``quantity.is_alternative`` has been replaced by ``quantity.name``, where ``name`` is the identifier of the main quantity this quantity is an alternative to. In general the alternative system now works in the following way:

- The main quantity (the one with the highest priority) has a list with ``alternatives``.
- All other quantities will if their identifier does not equal their ``quantity.name`` add their identifier to the ``quantity.names`` alternative list.
- When a quantity gets requested as part of node, the VaspParser will check, whether that quantity can be parsed. If not it will go through the list of alternatives and check those.

Setting the ``alternatives`` list on the main quantity is only required if a well defined sequence of priority is required by the developer. Otherwise the alternatives lists will be automatically set when loading the quantities.

That the quantities now have a built-in mapping between their quantity identifier and the main quantity they are an alternative to, helped resolving the issue of assigning the correct quantities to a node. E.g. the Fermi level will now always be 'fermi_level' in ``ParameterData`` and not sometimes ``outcar-fermi_level`` if it has been parsed from OUTCAR.

Open issues
-----------

- All the short-cut properties are very similar. There should be a way to generalise them and set them in ``__init__`` or somewhere else.
- The ``NodeComposer`` still depends a lot on the actual input format provided. @DropD suggested a more generic form for the ``NodeComposer`` that would look like

::

  example_data = {
      'kpoints_mesh': {
          'args': ([10, 10, 10]),
          'offset': [1, 1, 1]
      },
      'cell_from_structure': {
          'args': (my_structure)
      }
  }

  def make_node(node_type, data):
      import deepcopy
      new_node = node_type()
      for property, property_args in deepcopy(data).items():
          setter_name = 'set_{}'.format(property)
          setter = getattr(new_node, setter_name)
          args = property_args.pop('args', [])
          kwargs = property_args
          setter(*args, **kwargs)
      return new_node

Initially there is a concern was that this very strict requirement on the format of quantities would make writing new FileParsers more complicated. Or we would have to implement conversion functions that require approximately the same amount of code than the current solution. However, we are already now doing something very similar in requiring a specific format for the quantities. Therefore it would not even be that much of change.

.. _VASP: https://www.vasp.at
.. _AiiDA: https://www.aiida.net
