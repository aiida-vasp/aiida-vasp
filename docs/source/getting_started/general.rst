.. _general_notes:

================
1. General notes
================

To utilize `AiiDA-VASP`_ you will need to make sure that:

- You have a working `AiiDA`_ version >= 2 installation (see note).
- You have setup a `profile` in `AiiDA`_.
- That you hold a valid `VASP license`_ and that `VASP`_ is installed on some computer, for instance a remote HPC cluster.
- VASP >= 5.4.4 is used. The plugin has been tested with both VASP 5.4.4 and VASP 6 versions.
- You have defined a `computer` where `VASP`_ is installed and that you can SSH to that computer without using a password.

Since `AiiDA`_ is continuously evolving, we do not give details on how to install and configure it here. Please consult
the `AiiDA documentation`_ for details regarding this. In the documentation you will also find details on how to setup a `profile` and a `computer`.
`VASP`_ is licensed software and you need to obtain your own `VASP license`_. If you need to install `VASP`_ yourself or need
to assist someone, for instance HPC maintenance staff, please consult the `VASP wiki`_.

.. note::
   We do have a compatibility release that supports `AiiDA`_ version 1.6.4, but this is not maintained.
   Also, we strongly recommend users to move to an `AiiDA`_ version >= 2.

.. _VASP: https://www.vasp.at
.. _VASP license: https://www.vasp.at/sign_in/registration_form/
.. _VASP wiki: https://www.vasp.at/wiki/index.php/The_VASP_Manual
.. _AiiDA-VASP: https://github.com/aiida-vasp/aiida-vasp
.. _AiiDA: https://www.aiida.net
.. _AiiDA documentation: https://aiida.readthedocs.io/projects/aiida-core/en/latest/index.html
