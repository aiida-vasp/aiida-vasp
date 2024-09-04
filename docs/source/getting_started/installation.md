(installation)=
# Installation

Here we briefly describe how to install the [aiida-vasp] plugin. The plugin is available on PyPI and can be installed using `pip`. However, it is recommended to install the plugin from the source code if you plan to contribute to its development.

## Prerequisites

Before starting to use this plugin, please make sure the following prerequisites are met:

- A working [AiiDA] version >= 2.3 installation.
- A configured `profile` in AiiDA.
- VASP is installed on some computer, for instance a remote HPC cluster. Ideally also on the computer running this plugin so you can run some quick calculations to explore.
- VASP >= 5.4.4 is used. The plugin has been tested with both VASP 5.4.4 and VASP 6 versions.
- A `computer` where VASP is installed and that you can SSH to that computer without using a password. This can be your local computer.

For details on how to install and configure AiiDA, please consult its own online documentation.
In the documentation you will also find details on how to setup a `profile`, `Computer` and `Code`.

[VASP] is licensed software and you need to obtain your own [VASP license]. If you need to install [VASP] yourself or need to assist someone, for instance HPC maintenance staff, please consult the [VASP wiki].


## Install the plugin

In most cases, the AiiDA is installed into a virtual environment activate the virtual environment associated with it:

::::{tab-set}

:::{tab-item} venv

Assuming the virtual environment is installed in ``~/env/aiida-vasp``, activate it using:
```bash
   $ source ~/env/aiida-vasp/bin/activate
```
:::

:::{tab-item} conda

Assuming the conda environment is named ``aiida-vasp``, activate it using:
```bash
   $ conda activate aiida-vasp

```
:::

::::


The [aiida-vasp] plugin can now be installed using `pip`:

```bash
   $ (aiida-vasp) pip install aiida-vasp
```

However, it is likely that the pypi version is not up to date with the latest development version. In this case, you can install the plugin from the source code using:

```bash
   $ (aiida-vasp) git clone https://github.com/aiida-vasp/aiida-vasp.git
   $ (aiida-vasp) cd aiida-vasp
   $ (aiida-vasp) pip install -e .[pre-commit,testing]
```

This is also the recommended way to install the plugin if you plan to contribute to its development.

Note that in both cases, the dependencies of the plugin will also be installed automatically.

To verify the installation, you can run the following command:

```bash
   $ verdi plugin list aiida.calculations
```

and the printed list should include the  `vasp.vasp` entry.



[aiida-vasp]: https://github.com/aiida-vasp/aiida-vasp
[vasp]: https://www.vasp.at
[pymatgen]: https://pymatgen.org
[aiida-quantumespresso]: https://github.com/aiidateam/aiida-quantumespresso
[aiida-castep]: https://github.com/zhubonan/aiida-castep
[aiida]: https://www.aiida.net
[online documentation]: https://aiida.readthedocs.io/projects/aiida-core/en/latest/index.html
[vasp license]: https://www.vasp.at/sign_in/registration_form/
[vasp wiki]: https://www.vasp.at/wiki/index.php/The_VASP_Manual
