(general-notes)=

# 1. General notes

To utilize [AiiDA-VASP] you will need to make sure that:

- You have a working [AiiDA] version >= 2 installation (see note).
- You have setup a `profile` in [AiiDA].
- That you hold a valid [VASP license] and that [VASP] is installed on some computer, for instance a remote HPC cluster.
- VASP >= 5.4.4 is used. The plugin has been tested with both VASP 5.4.4 and VASP 6 versions.
- You have defined a `computer` where [VASP] is installed and that you can SSH to that computer without using a password.

Since [AiiDA] is continuously evolving, we do not give details on how to install and configure it here. Please consult
the [AiiDA documentation] for details regarding this. In the documentation you will also find details on how to setup a `profile` and a `computer`.
[VASP] is licensed software and you need to obtain your own [VASP license]. If you need to install [VASP] yourself or need
to assist someone, for instance HPC maintenance staff, please consult the [VASP wiki].

:::{note}
We do have a compatibility release that supports [AiiDA] version 1.6.4, but this is not maintained.
Also, we strongly recommend users to move to an [AiiDA] version >= 2.
:::

[aiida]: https://www.aiida.net
[aiida documentation]: https://aiida.readthedocs.io/projects/aiida-core/en/latest/index.html
[aiida-vasp]: https://github.com/aiida-vasp/aiida-vasp
[vasp]: https://www.vasp.at
[vasp license]: https://www.vasp.at/sign_in/registration_form/
[vasp wiki]: https://www.vasp.at/wiki/index.php/The_VASP_Manual
