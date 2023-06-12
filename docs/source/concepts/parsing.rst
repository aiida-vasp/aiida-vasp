.. _parsing:

=======
Parsing
=======
AiiDA-VASP provides flexible parsing of `VASP`_ output files to store data in the `AiiDA`_ database and repository.

The quantities that can be parsed are now fully customizable. The user interface for configuring the parsing settings takes place in the ``settings['parser_settings']`` dictionary entry. The default ``parser_settings`` is presently:

.. warning::
   Notice however, that even though the parser and the `node composer` is configurable, the output check in `AiiDA`_ will complain that your newly added custom node is not detected in the ``spec`` definitions of for instance your :py:class:`VaspCalculation<aiida_vasp.calcs.vasp.VaspCalculation>`. If you do add additional nodes outside the ones defined in the code already, please consider to also add it to the ``spec.output`` section in the :py:class:`VaspCalculation<aiida_vasp.calcs.vasp.VaspCalculation>` class and potentially also to workchains if need be.

.. literalinclude:: ../../../aiida_vasp/parsers/vasp.py
   :start-after: ParserSettings
   :end-before: VaspParser

where ``<node_name>`` is defined as ``add_<node_name>``. The ``<node_name>`` usually correspond to
the parsed ``quantities``, except for ``misc`` which is a container for all system size independent quantities.
This can be seen from the current default node definitions:

.. literalinclude:: ../../../aiida_vasp/parsers/settings.py
   :start-after: NODES
   :end-before: ParserDefinitions

As you can see, the ``<node_name>`` named ``bands`` is composed of three quantities,
the ``eigenvalues``, the ``kpoints`` and the ``occupancies``. You typically need all three
when you analyze the band structure. Similarly it is possible to customize the output stored by
composing different ``quantities``. However, the regular user should not need to utilize these
functions. It is more useful when developing new workchains, where it makes sense to introduce a
new output container which does not already exists.

There are four ways to interact and set the parser properties.

#. Using a boolean::

     settings['parser_settings'] = {'add_<node_name>': True}

  This will enable the parsing of ``<node_name>``. This is the recommended way to interact with the current parser.

#. Specifying a list::

     settings['parser_settings'] = {'add_<node_name>': ['<quantity1>', '<quantity2>', ... ]'}

   Where the list contains strings with quantity names that should be used for this node.
   This will work only for predefined nodes and available ``quantities``. This options is not for regular users.

#. Specifying a dict::

     settings['parser_settings'] = {'add_<node_name>': { 'type': '<node_type>', 'quantities': ['<quantity1>', '<quantity2>', ...], 'link_name': '<link_name>'}}

   This will define a new custom node. The ``link_name`` is optional and will be set to the ``<node_name>`` if it is not present.
   This option is not for regular users.

#. Defining the parsing interactively

   The three ways above were all mediated by modifying ``settings``. There is also an interactive way,
   if you have an instance of :py:class:`VaspParser<aiida_vasp.parsers.vasp.VaspParser>`::

     VaspParser.add_custom_node(node_name, node_definition)

   where the format for the ``node_definition`` is as in the previous example with the custom nodes.

The parser will check for any error detected for the underlying VASP calculation.
This information is stored in the `notification` quantity, which contains a list of error/warnings detected.
By default, a non-zero exit state is returned if any critical error is found.
The default list of critical errors is defined under the `critical_notifications` key inside :py:data:`~aiida_vasp.parsers.vasp.DEFAULT_SETTINGS`.
Additional settings may also be supplied under ``parser_settings`` to modify the behaviors:

* ``ignore_all_errors``: a boolean value to control whether all notifications should be ignored, defaults to ``False``.
* ``critical_notifications``: a dictionary with items like ``'add_<error>': True`` for controlling which errors to be regarded as `critical` or not.

Composing the quantities into an output node
--------------------------------------------

A :py:class:`NodeComposer<aiida_vasp.parsers.node_composer.NodeComposer>` has been added to compose output nodes based on a given set of ``quantities``. It can be initialized in two different ways::

  composer = NodeComposer(vasp_parser=...)
  composer = NodeComposer(content_parsers=[...])

where either a :py:class:`VaspParser<aiida_vasp.parsers.vasp.VaspParser>` object is passed to ``vasp_parser`` or a list of specific object parsers is passed to ``content_parsers``. Either are required for the :py:class:`NodeComposer<aiida_vasp.parsers.node_composer.NodeComposer>` to assemble the required quantities for a node. Here is an example for the interface for composing a node from the ``VasprunParser`` test::

  composer = NodeComposer(content_parsers=[vasprun_parser])
  data_obj = composer.compose('array.kpoints', quantities=['kpoints'])

For most node types there are default quantities defined as::

  .. literalinclude:: ../../../aiida_vasp/parsers/node_composer.py
		      :start-after: NODE_TYPES
		      :end-before: NodeComposer

so that the above could be shortened to::

  data_obj = composer.compose('array.kpoints')

Short-cut properties
--------------------

Since none of the ``quantities`` will actually return an AiiDA output node anymore, there are some built-in short-cut properties
in order to save some typing and also to avoid using the ``get_quantity`` interface directly. An example would be::

  parser = PoscarParser(file_path=path)
  result = parser.structure

where the whole composer part is dealt with internally::

  @property
  def structure(self):
      if self._structure is None:
          composer = NodeComposer(content_parsers=[self])
          self._structure = composer.compose('structure', quantities=['poscar-structure'])
      return self._structure

At the moment existing short-cut properties are:

- ``PoscarParser.structure``
- ``KpointsParser.kpoints``
- ``OutcarParser.parameter``
- ``DoscarParser.dos``
- ``WavecarParser.wavecar``
- ``ChgcarParser.chgcar``
- ``IncarParser.incar`` (This is an exception in so far that it actually just returns a dictionary)

Prominently missing from this list are ``bands`` from the ``EigenvalParser`` and everything from the ``VasprunParser``.
For the former the reason is that ``bands`` require quantities from other object parsers and for the latter the name space is crowded.

Quantity definitions and alternatives
-------------------------------------

For alternatives, ``quantity.name`` can be defined, where ``name`` is the identifier of the main quantity this
quantity is an alternative to. In general the alternative system now works in the following way:

- The main quantity (the one with the highest priority) has a list with ``alternatives``.
- All other quantities will if their identifier does not equal their ``quantity.name`` add their identifier to
  the ``quantity.names`` alternative list.
- When a quantity gets requested as part of node, the ``VaspParser`` will check, whether that quantity can be parsed.
  If not it will go through the list of alternatives and check those.

Setting the ``alternatives`` list on the main quantity is only required if a well defined sequence of priority
is required by the developer. Otherwise the alternatives lists will be automatically set when loading the quantities.

The quantities now have a built-in mapping between their quantity identifier and the main quantity they are an
alternative to, helps resolving issues of assigning the correct quantities to a node.
E.g. the Fermi level will now always be ``fermi_level`` in ``misc`` and not sometimes
``outcar-fermi_level`` if it has been parsed from OUTCAR.

.. _VASP: https://www.vasp.at
.. _AiiDA: https://www.aiida.net
