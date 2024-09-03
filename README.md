# AiiDA VASP plugin

| Release |  [![PyPI](https://img.shields.io/pypi/v/aiida-vasp)](https://pypi.org/project/aiida-vasp/) | [![PyPI](https://img.shields.io/pypi/status/aiida-vasp )](https://pypi.org/project/aiida-vasp/)|
|:--------|:------ |:----|
| Build   | [![Coverage](https://codecov.io/gh/espenfl/aiida-vasp/branch/master/graph/badge.svg)](https://codecov.io/gh/espenfl/aiida-vasp)| [![Docs](https://readthedocs.org/projects/aiida-vasp-plugin/badge/?version=latest)](http://aiida-vasp-plugin.readthedocs.io/en/latest/?badge=latest)  |
| Stats  | [![Dowloads](https://img.shields.io/pypi/dm/aiida-vasp)](https://pypi.org/project/aiida-vasp/) | [![Commits]( https://img.shields.io/github/commit-activity/m/aiida-vasp/aiida-vasp)](https://github.com/aiida-vasp/aiida-vasp/commits/develop) |


This is a plugin to [AiiDA] to run calculations with the ab-initio program [VASP].

Please have a look at the [AiiDA-VASP documentation] for instructions on how to install and use the plugin.

## Installing the plugin

% Keep this comment as it is used for including these steps in the install section of the docs.
% It includes everything past the next line.
% Start installation description

1. If you are already using [AiiDA], simply activate the virtual environment associated with it, here assumed to be located in `~/env/aiida-vasp`:

   ```
   $ source ~/env/aiida-vasp/bin/activate
   ```

2. Otherwise, set up a new virtual environment:

   ```
   $ python -m venv ~/env/aiida-vasp
   ```

3. And then enable the newly installed virtual environment:

   ```
   $ source ~/env/aiida-vasp/bin/activate
   ```

4. Install the [AiiDA-VASP] plugin (and [AiiDA] if that is not already installed):

   ```
   $ (aiida-vasp) pip install aiida-vasp
   ```

> If you need to install the compatibility release of [AiiDA-VASP] which works with [AiiDA] 1.6.4 you should instead install the plugin
> using `pip install aiida-vasp=2.2`, but this is not recommended and only mentioned for legacy support. For the legacy version you
> also most likely have to run `reentry scan -r aiida` after installing the plugin.

This will automatically install the [AiiDA] python package(s) as well as any other dependencies of the plugin and register all the plugin classes with [AiiDA].

Please consider that [AiiDA] have prerequisite that needs to be installed and ensured working. The steps above will not take care of this for you. Please consult [AiiDA prerequisites] and follow the instructions therein.

% End installation description

## Support

% Start support description

The development, maintenance and use of this plugin is considered a community effort. In order to facilitate for the community to contribute,
we have established a [space on Matrix] that users can use to communicate. We encourage users to help each other. In addition,
the development team is present in the space and users are free to ask.
First consult the documentation of both [AiiDA-VASP documentation] and [AiiDA documentation] and also consider that the developers are
not paid for this work. Please respect potential lead times in getting answers and be polite.

% End support description

[aiida]: https://www.aiida.net
[aiida documentation]: http://aiida-core.readthedocs.io/en/latest/
[aiida prerequisites]: https://aiida-core.readthedocs.io/en/latest/install/prerequisites.html
[aiida-vasp]: https://github.com/aiida-vasp/aiida-vasp
[aiida-vasp documentation]: https://aiida-vasp-plugin.readthedocs.io/en/latest/
[space on matrix]: https://matrix.to/#/#aiida-vasp:matrix.org
[vasp]: https://www.vasp.at
