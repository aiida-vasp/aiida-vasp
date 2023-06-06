.. _docs_for_docs:

=============
Documentation
=============

We honor the importance of an up to date and useful documentation and would like to
motivate also the readers to contribute. Here we will describe how one
can

Contributing to the documentation
---------------------------------
Please consult the `open issue list`_ labeled documentation and see if there is
already an issue for what you find lacking or not good enough. If not, please open an issue and tag it
with documentation. The community also greatly appreciate contributions to solve issues, by
submitting a pull request with suggested changes. Before doing so, please check the
documentation by building it on your local computer and verifying that the Sphinx generated documentation
in your browser displays what you expect.

Building the documentation locally on your computer
---------------------------------------------------

Make sure your existing AiiDA-VASP virtual environment is enabled.
Install dependencies from the AiiDA-VASP root folder by issuing
the following command::

  $ pip install -e .[docs]

With a temporary web server
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Now, we can build the documentation and setup a temporary web
server for us to browse and check the documentation. From the
AiiDA-VASP root folder, issue the command::

  $ sphinx-autobuild docs/source docs/_build

Now open ``http://127.0.0.1:8000`` in a browser and inspect the documentation.
The webserver can be brought down by killing the ``sphinx-autobuild`` process,
i.e. by pressing ``Ctrl+C`` or similar.

Without a web server
~~~~~~~~~~~~~~~~~~~~

To just build the documentation without setting up a temporary web server execute::

  $ sphinx-build docs/source docs/_build

With a link checker
~~~~~~~~~~~~~~~~~~~

We can also check links in our documentation, but issuing::

  $ sphinx-build -b linkcheck docs/source docs/_build

.. _open issue list: https://github.com/aiida-vasp/aiida-vasp/issues?q=is%3Aissue+is%3Aopen+label%3A%22documentation%22
