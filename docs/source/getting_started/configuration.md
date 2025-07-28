---
myst:
  substitutions:
    InstalledCode: "{external:py:class}`InstalledCode <aiida.orm.InstalledCode>`"
---

(configuration)=
# Configurations

Here we explain the post-installation configuration steps for [AiiDA-VASP].

## Setting up a {{ InstalledCode }} for your VASP executable

A {{ InstalledCode  }} is a pointer to a VASP executable that is installed on a remote computer. In this example, we assume that the VASP executable is installed on a remote computer `mycluster`. We will now set up a {{ InstalledCode }} for this executable in [AiiDA].

The `verdi code` commands allow one to setup/update/duplicate `InstalledCode` objects in [AiiDA].
Please consult the [AiiDA documentation] for more details.

Before configuring the code in [AiiDA], please make sure it runs and functions as normal on the computer `mycluster`.

To add a code with label `vasp`:

```
% verdi code create core.code.installed
Report: enter ? for help.
Report: enter ! to ignore the default and set no value.
Computer: mycluster
Absolute filepath executable: /cluster/software/vasp/vasp6.3.2/vasp_std
Label: vasp-std
Description: VASP 6.3.2 standard version (complex)
Default `CalcJob` plugin: vasp.vasp
Escape using double quotes [y/N]:
Success: Created InstalledCode<6>
```

The `Absolute filepath executable` is the full path to the VASP executable installed on the remote computer.
Very often one needs to utilize different versions/executables VASP, for instance, running calculations with the gamma only or non-collinear configurations or with additional auxiliary libraries, like BEEF included.
One can add multiple `InstalledCode` objects to AiiDA for different versions of VASP with different labels.

:::{tip}
You can first create a Code for the `vasp_std` executable, then use:

```
verdi code duplicate vasp-std@mycluster vasp-ncl@mycluster
```

to create a new {{ InstalledCode }} for `vasp_ncl`. The previous configurations will be used as default and the only one you need to change is the *Absolute filepath executable* field.
:::

During the end of the setup, the user will be asked to enter the prepend and append text.
The prepend section is any command that should run before the VASP executed.
For most cluster systems, they corresponds to loading the correct modules.

Enter something along the lines of:

```
module purge
module load <myvaspmodule>
```

in the first open section of the file.

The `Absolute filepath executable` can be obtain with for instance `module show <myvaspmodule>` under the `PATH` environment variable.

It is important to check what the name of the actual VASP executable is as this could be tailored by your HPC maintainers.
The default VASP build system yields the `vasp_std`, `vasp_ncl` and `vasp_gam`, which is maybe a good start.
After saving this file, a new file opens, which is for the text to be appended, e.g. what is done after the executable has been executed.
One may want to enter cleanup routines etc, or just keep it empty.


:::{note}
For local VASP installation, just create a computer with `localhost` as the hostname, `core.local` transport and `core.direct` scheduler.
The `localhost` computer should be created already if you used the `verdi presto` command to initialize the profile.
Then create a {{ InstalledCode }}  with the absolute path to the VASP executable and label `vasp`.
:::


The information about the code can be inspected `verdi code show vasp@mycluster`.

```bash
% verdi code show vasp@mycluster
--------------  -------------------------------------
PK              6
UUID            e5de4c5a-ca44-4a3b-a42c-e5f7e1c21cbb
Label           vasp
Description     VASP 6.3.2 standard version
Default plugin  vasp.vasp
Prepend text    module purge
                module load <myvaspmodule>
Append text
--------------  ------------------------------------
```

## Configure pseduopotentials (POTCARs)

To run a VASP calculation, potentials (the POTCAR files) have to be uploaded to
the database. For more details regarding the handling of the potentials, please see [potentials][#potentials].

:::{hint}
If you have `pymatgen` fully configured with POTCAR data, they can be uploaded to AiiDA-VASP.
See [this section](#potcar-from-pymatgen) for the details.
:::

Assuming you already have a valid license you can download the VASP potentials from their portal and save
it to a convenient location.
For the example here, we use `$HOME/myaiida/potpar_PBE.54.tar`.

In [AiiDA-VASP] we refer to a set of potentials as a potential family, which are typically a
specific version of the PAW dataset relaxed with VASP.


Execute the following command to upload the whole potential family to the database:

```
% verdi data vasp.potcar uploadfamily --path=$HOME/myaiida/potpaw_PBE.54.tar --name=PBE.54 --description="PBE potentials version 54"
POTCAR files found: 327. New files uploaded: 327, Added to Family: 327
```

:::{tip}
For historical reasons we used `PBE.54` as the name for the `PBE_54` PAW dataset.
Any name can be used of course, but some of the presets of the workflows expect the naming convention with
a `.` instead of `_`.
:::


The `name` and `description` are not optional and have to be specified.
The `path` could be either an archive, or one could use a folder name.
It is also possible, not to specify path, where the plugin will traverse all folders from the folder in which the command above is executed from and try to upload all potentials it finds to the specified family.


We use the POTCAR parser from [pymatgen] to get the metadata and sometimes this issues a warning if it detects unknown metadata flags in the potentials. You can usually ignore these warnings.


:::{note}
A *potential family* here can have multiple potential for a single element.
This is different from the same terminologies used in [aiida-pseudo] and [aiida-castep] where a *potential family* provides a one-to-one mapping between elements and pseudopotentials.
Hence, a *potential mapping* is also need when running calculations/workflows, to provided the one-to-one  mapping between the symbols and the potentials with the *potential family*.
:::


[aiida-vasp]: https://github.com/aiida-vasp/aiida-vasp
[vasp]: https://www.vasp.at
[pymatgen]: https://pymatgen.org
[aiida-quantumespresso]: https://aiida-pseudo.readthedocs.io/en/latest/design.html#families
[aiida-castep]: https://github.com/zhubonan/aiida-castep
