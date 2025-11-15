(potentials)=

# How to import and use pseudopotentials

In this guide we will go through how to import and use the pseudopotentials (PAW datasets) in aiida-vasp.


As we know, any [VASP] calculation relies on potentials provided by the POTCAR files.

The actual POTCAR in the calculation should be concatenated from multiple POTCARs provided by the dataset in the same order as the symbols appear in the POSCAR.
Incorrectly ordering or misuse of the POTCAR files are common mistakes when deploying VASP calculation *by hand*.

When using aiida-vasp, the POTCAR files are managed automatically.
First, we need to make them available to AiiDA-VASP.
These are usually supplied by the [VASP] team and is part of the license.
Each POTCAR dataset is an archive (`.tar.gz` or `tgz`) containing multiple subfolders,
the name of the sub-folder is typically referred as the **symbol** of the POTCAR file it contains.

:::{hint}
Since PAW data files are all named as  `POTCAR`, the only way to tell them apart is to look at the content.
However, early releases of the VASP PAW dataset contains mistakes such that the **symbol** of the POTCAR
is not what it say in the file! This will disrupt the routines in aiida-vasp to assign the correct POTCAR
for each element.
:::

While aiida-vasp does allow importing only the sets (or even individual potentials) you require,
it is more common to import the entire dataset and keep them grouped as a potential family.


[AiiDA] does more than prepare calculations and send them to a cluster. The main focus of [AiiDA] lies on tracking data provenance.
Importing the POTCAR files into your working [AiiDA] database yields some advantages:

- aiida-vasp stores a unique hash for each file. This can help users navigate when different potentials have very similar looking headers, but do in fact contain a different potential.
- POTCAR files uploaded to the database cannot be modified accidentally, thus it is recorded unambiguously, which file was used for which execution of each run.
- Storing the file's contents rather than a link prevents accidentally breaking the link by moving the file away (or renaming it).

## How to import a set of POTCAR files

The command line tools for these tasks can be called through the `aiida-vasp` command:

```
$ (aiida-vasp) aiida-vasp potcar --help
Usage: aiida-vasp potcar [OPTIONS] COMMAND [ARGS]...

  Top level command for handling VASP POTCAR files.

Options:
  -v, --verbosity [notset|debug|info|report|warning|error|critical]
                                  Set the verbosity of the output.
  -h, --help                      Show this message and exit.

Commands:
  exportfamily              Export a POTCAR family into a compressed tar...
  fix-inconsistent-symbols  Fix inconsistent families
  integrity                 Check the integrity of a POTCAR family
  listfamilies              List available families of VASP potcar files.
  listsymbols               List available symbols in a POTCAR family group.
  migratefamilies           Migrate the type_string associated with the...
  upload-from-pymatgen      Upload a family of VASP potcar files from...
  uploadfamily              Upload a family of VASP potcar files.
```

```
$ aiida-vasp potcar uploadfamily --path=<path> --name=<potential_family> --description=<desc>
```

Where `<path>` is the path to the folder or tar archive containing the POTCAR set. The command expects the folder or archive to look like:

```
<path>/
|
+- Ac/
|  +- POTCAR
|  +- ...
|
+- Ag/
|  +- POTCAR
|  +- ...
...
```

If it encounters anything different, it will recursively search the given path for subpaths matching this structure and import all the POTCAR files found in that way.

`<potential_family>` is the label you will use to access the potentials from this set or to specify which potentials you want to use in a particular [VASP] run. The meaning of `<description>` is self-explanatory.

Custom sets can simply be arranged in a matching folder structure and then imported using the same command.

## Uploading a set of potentials incrementally

For this purpose, we can use that the `uploadfamily` command by default adds any POTCAR files not yet uploaded to the family of the given `name`, for example:

```
$ aiida-vasp potcar uploadfamily --path=path/to/Ac --name="PBE_custom" --description="A custom set"
$ aiida-vasp potcar uploadfamily --path=other/path/to/Ag --name="PBE_custom"
```

Note, that the description does not have to be given if the family already exists.

Due to the recursive nature of the search, this also works for combining several small sets of POTCARs in a few commands, without having to arrange them in a different way first.

## How to check what potential families are present in the database

```
$ aiida-vasp potcar listfamilies
```

## How to access uploaded potentials and search

The data structure used to find and retrieve potentials is called {py:class}`PotcarData<aiida_vasp.data.potcar.PotcarData>` and can be accessed through AiiDA's data factory as `DataFactory('vasp.potcar')`. This class provides shortcuts for simple or frequent cases, for complex queries, please refer to the [AiiDA documentation] on querying the provenance graph.

More advanced searches, like for ranges of properties etc can be done using the {py:class}`QueryBuilder<aiida.orm.querybuilder.QueryBuilder>` tool, which is part of [AiiDA] and documented there.

Use:

```
PotcarData.find(<property>=<value>, <property2>=<value2>, ...)
```

which returns a list of all stored {py:class}`PotcarData<aiida_vasp.data.potcar.PotcarData>` instances fulfilling the criteria. Some important supported `<property>` entries are:

> - `sha512` - An SHA512 hash of the file contents
> - `title` - Title of the potential, typically the title of the POTCAR
> - `element` - The chemical element described by this potential
> - `full_name` - The name of the containing folder from which it was uploaded. This is used to specify a potential inside a family. Example: `Zn_sv_GW`
> - `original_file_name` - The filename (+ last three directories) from which it was uploaded (May help identifying exactly where it came from).

and for each you supply the `<value>` which is relevant for you given search.

To find one potential for each element in a list of element names, all from the same family:

```
mapping = {
   'Ac': 'Ac',
   'Ag': 'Ag_GW'  # or 'Ag_pv', 'Ag_sv_GW', ...
}
potcars_for_elements = PotcarData.get_potcars_dict(
   elements=['Ac', 'Ag', ..], <potential_family>, mapping=mapping)
```

The `mapping` dictionary is required to decide which of the variants should be chosen for each element. The mapping can also conveniently be stored in a {py:class}`Dict<aiida.orm.nodes.data.dict.Dict>` node for reuse. The potential family is specified with `<potential_family>`.

## How to pass potentials to a VASP calculation

For a single [VASP] calculation run, you should at the very minimum use the VaspWorkChain, which takes the family as a database-storable string and a dictionary mapping elements to a particular variant for that element:

```
from aiida.plugins import DataFactory
from aiida.common.extendeddicts import AttributeDict
from aiida.orm import Str

inputs = AttributeDict()
inputs.potential_family = Str('<potential_family>')
inputs.potential_mapping = DataFactory('dict')(dict={'In': 'In_d', 'As': 'As'})
```

The VaspWorkChain takes care of finding the right files and concatenating them for you.

For a more complex workflow, the process may be different, it may for example use heuristics to find a default potential for you.

[aiida]: https://www.aiida.net
[aiida documentation]: http://aiida-core.readthedocs.io/en/latest/
[vasp]: https://www.vasp.at
