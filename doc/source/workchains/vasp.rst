.. _vasp_workchain:

==============
VASP workchain
==============

Required inputs
^^^^^^^^^^^^^^^

 The VaspWorkChain requires a number of inputs, these comprise the minimum set of information to run a `VASP`_ calculation from `AiiDA`_.

* ``code``: an `AiiDA`_ :py:class:`aiida.orm.nodes.data.Code` instance, describes the VASP executable and holds a reference to the ``Computer`` instance on which it lives.
* ``structure``: an `AiiDA`_ :py:class:`aiida.orm.nodes.data.StructureData` or :py:class:`aiida.orm.nodes.data.CifData` instance, describes the structure on which VASP is to be run.
* ``kpoints``: an `AiiDA`_ :py:class:`aiida.orm.nodes.data.KpointsData` instance, describing the kpoints mesh or path.
* ``potential_family``: an `AiiDA`_ :py:class:`aiida.orm.nodes.data.Str` instance, the name given to a set of uploaded POTCAR files.
* ``potential_mapping``: an `AiiDA`_ :py:class:`aiida.orm.nodes.data.Dict` instance, containing an entry for at least every kind name in the ``structure`` input with the full name of the POTCAR from the ``potential_family``. Example: ``{'In1': 'In_d', 'In2': 'In_h'}``.
* ``parameters``: an `AiiDA`_ :py:class:`aiida.orm.nodes.data.Dict` instance, containing key/value pairs that is later written to an INCAR file as ``KEY =VALUE``. Keys and values can be lower case and the standard built-in Python data types should be used for values.
* ``options``, an `AiiDA`_ :py:class:`aiida.orm.nodes.data.Dict` instance, containing at least the keys ``resources``. More information about the options is available in the `AiiDA documentation`_.
		       
.. _AiiDA: https://www.aiida.net
.. _VASP: https://www.vasp.at
.. _AiiDA documentation: http://aiida-core.readthedocs.io/en/latest/
.. _Workchain: https://aiida.readthedocs.io/projects/aiida-core/en/latest/concepts/workflows.html#work-chains

