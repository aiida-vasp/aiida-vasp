.. _running_tests:

=============
Running tests
=============

`AiiDA-VASP`_ utilizes `pytest`_ which makes it easy to run and construct new tests. We also utilize the fixture
managers that are present, see the `AiiDA documentation`_.

In order to run all tests execute, in the root folder::

  $ pytest

This is a useful test to see if the plugin works and that everything is installed correctly.

.. note::
   For contributors, one principle rule is in place: always try to construct tests when submitting a PR.

.. note::
   We also utilize ``tox`` in order to make the test environment fully isolated. In order to run tests
   in ``tox``, please issue ``tox`` in the root folder of `AiiDA-VASP`_.

.. _pytest: https://docs.pytest.org/en/latest/
.. _AiiDA-VASP: https://github.com/aiida-vasp/aiida-vasp
.. _AiiDA documentation: https://aiida.readthedocs.io/projects/aiida-core/en/latest/index.html
.. _VASP: https://www.vasp.at/
