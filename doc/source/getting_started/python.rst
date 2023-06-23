.. _python_env:

======================================
2a. Setup a regular Python environment
======================================

We will now set up a standard Python virtual environment called ``aiida-vasp``.

.. include:: ../../../README.rst
   :start-line: 22
   :end-line: 35

Please make sure that your virtual environment is activated at all times when performing
operations that depends on it, or you want to install into it.

All relevant `AiiDA`_ configurations and the repository reside in a dedicated directory.
We would recommend that you put this inside your directory of the virtual environment in
order to keep everything pertaining to this environment in one location. In order to make
sure `AiiDA`_ detects this, we need to set the ``AIIDA_PATH`` environmental variable. This can be
done by exporting it in the active script of your virtual environment::

  $ echo "export AIIDA_PATH=~/env/aiida-vasp" >> ~/env/aiida-vasp/bin/activate

where we have assumed you virtual environment is located in ``~/env/aiida-vasp``.

In addition to this, we would like to enable tab completion, which you can do by::

  $ echo 'eval "$(_VERDI_COMPLETE=source verdi)"' >> ~/env/aiida-vasp/bin/activate
