.. _conda_env:

===============================
2b. Setup the Conda environment
===============================

We assume `Conda`_ is installed and operational and will now set up a `Conda`_ virtual
environment called ``aiida-vasp`` and prepare it::

   % cd ~
   % mkdir aiida-vasp
   % cd aiida-vasp
   % conda create -n aiida-vasp python=3
   % source activate aiida-vasp
   % conda install gcc_linux-64 gxx_linux-64 ipython

The command ``conda activate aiida-vasp`` enables the `Conda`_ environment ``aiida-vasp``.
All the settings, installs or running are always done in this `Conda`_
environment. It is mandatory to activate this every time logging to
this computer. Please also remember to activate this everytime you
need to make changes to the respective environments. You can of course
created multiple `Conda`_ environments depending on your use case or versions
you would like to have installed.

To enable tab completions, the file
``$CONDA_PREFIX/etc/conda/activate.d/env_vars.sh`` needs to be created
with the following content:

::

   #!/bin/sh
   export AIIDA_PATH=~/aiida-vasp
   eval "$(_VERDI_COMPLETE=source-bash verdi)"
   export HOST=`hostname`

and ``$CONDA_PREFIX/etc/conda/deactivate.d/env_vars.sh`` containing::

   #!/bin/sh
   export HOST=`hostname`

In order to enable the recently placed files, we need to reactivate
the Conda environment by issuing the following:

::

   % conda deactivate
   % conda activate aiida-vasp

.. _AiiDA-VASP: https://github.com/aiida-vasp/aiida-vasp
.. _Conda: https://docs.conda.io/en/latest/
.. _AiiDA: https://www.aiida.net
