.. _general_notes:

================================
1. General notes about the setup
================================

Please consider that `AiiDA`_ have prerequisites that needs to be installed and ensured working. Consult `AiiDA prerequisites`_ and follow the instructions therein. Typically, on a standard Linux installation, like a modern `Ubuntu`_, `Fedora`_ or similar distribution, the only thing you need to install on top of a development workstation install of said Linux distributions is `RabbitMQ`_, in addition to installing a `PostgreSQL`_ server and setting up a `PostgreSQL`_ database on the server. Both the `RabbitMQ`_ and `PostgreSQL`_ server can be hosted on a different computer than the one you run `AiiDA`_ on, but we will not cover such installation procedures in this documentation.

Please consider to work on a bleeding distribution of `Fedora`_ or `Ubuntu`_. You can make it work on `CentOS`_ or other long term stable releases, but it requires more knowledge of the users to update system packages.

`AiiDA`_ and `AiiDA-VASP`_ relies heavily on the Python framework and it is thus recommended to install it into a virtual environment that is dedicated to it. If you already use Python chances are you have a working system to manage your virtual environments, maybe even a favorite. Note that `AiiDA`_ as well `AiiDA-VASP`_ are tested with both ``virtualenv``, ``conda`` and ``virtualenvwrapper``. Since Python 2 is soon to be deprecated the documentation is written assuming Python 3.6 or later is present and working. However, `AiiDA`_ and `AiiDA-VASP`_ are compatible with Python 2 in case you prefer (we suggest you not). If that is the case you should be able to translate the procedures that are different to Python 2. As previously :ref:`mentioned <convention>` we will try to give commands valid for both the regular Python ``virtualenv`` and `Conda`_ virtual environments. One reason to use `Conda`_ is that it is easy
to manage and install a `PostgreSQL`_ database as a non-privileged user. However, the user still needs to be a privileged user to install `RabbitMQ`_. We thus believe the choice of using the standard ``virtualenv`` or `Conda`_ is more dependent on other activities you are doing.


When the virtual environment, the required prerequisites for `AiiDA`_ and `AiiDA`_ itself are installed and configured, a `profile`, a `computer` and a `code` needs to be configured in `AiiDA`_.

The documentation gives information on how to go from a clean system to the very end, which is to launch a test calculation using `AiiDA-VASP`_. If there is a difference in the system commands between Linux distros, all commands we give are based on `Fedora`_ 30.

If at any point you need more information, please consult the general `AiiDA-VASP documentation`_. In case you have already configured a `profile` and a `computer` please skip to the step which creates a `code`. It is highly unlikely that you need help with the installation of `AiiDA-VASP`_ and already have a `VASP`_ code configured in `AiiDA`_. In case that is indeed so, please skip this step and launch your first test calculation and then consider to continue on with the tutorial steps to learn more about the plugin.

.. _VASP: https://www.vasp.at
.. _AiiDA-VASP documentation: https://aiida-vasp.readthedocs.io/en/latest/
.. _AiiDA-VASP: https://github.com/aiida-vasp/aiida-vasp
.. _Conda: https://docs.conda.io/en/latest/
.. _CentOS: https://www.centos.org/?
.. _Fedora: https://getfedora.org/
.. _Ubuntu: https://ubuntu.com/
.. _PostgreSQL: https://www.postgresql.org/
.. _RabbitMQ: https://www.rabbitmq.com/
.. _AiiDA: https://www.aiida.net
.. _AiiDA prerequisites: https://aiida-core.readthedocs.io/en/latest/install/prerequisites.html
