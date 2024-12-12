# Welcome to AiiDA-VASP's documentation!

AiiDA-VASP is a plug-in for the workflow management and data provenance tracking framework [AiiDA]. It provides the classes [AiiDA] needs to run simulations using [VASP] (Vienna Ab initio Simulation Package). [VASP] is a program for atomic scale materials modelling, e.g. electronic structure calculations and quantum-mechanical molecular dynamics, from first principles. For detailed documentation on using [VASP] take a look in their [VASP wiki].

______________________________________________________________________


% Start of the card definitions
::::{grid} 1 2 2 2
:gutter: 3

:::{grid-item-card} {fa}`rocket;mr-1` Get started
:text-align: center
:shadow: md

Instructions to install, configure and setup the plugin package.

+++

```{button-ref} getting_started/index
:ref-type: doc
:click-parent:
:expand:
:color: primary
:outline:

Getting started with aiida-vasp.
```
:::


:::{grid-item-card} {fa}`bookmark;mr-1` Tutorials
:text-align: center
:shadow: md

Tutorials for running calculations and workflows.

+++

```{button-ref} tutorials/index
:ref-type: doc
:click-parent:
:expand:
:color: primary
:outline:


To the Tutorials.
```
:::


:::{grid-item-card} {fa}`question-circle;mr-1` Concepts
:text-align: center
:shadow: md

Learn the key concepts in the plugin.

+++

```{button-ref} concepts/index
:ref-type: doc
:click-parent:
:expand:
:color: primary
:outline:


Learn the concepts.
```
:::

:::{grid-item-card} {fa}`cogs;mr-1` How-to guides
:text-align: center
:shadow: md

Guides for achieving specific goals.
+++

```{button-ref} howto/index
:ref-type: doc
:click-parent:
:expand:
:color: primary
:outline:


How-to guides.
```
:::


:::{grid-item-card} {fa}`cogs;mr-1` Workflows
:text-align: center
:shadow: md

Automate your VASP calculations with pre-defined workflows.

+++



```{button-ref} workflows/index
:ref-type: doc
:click-parent:
:expand:
:color: primary
:outline:


Learn to use workflows.
```
:::

:::{grid-item-card} {fa}`cogs;mr-1` Development
:text-align: center
:shadow: md

Information on how to develop the plugin.

+++

```{button-ref} developments/index
:ref-type: doc
:click-parent:
:expand:
:color: primary
:outline:


Learn to contribute.
```
:::


::::

% End of the card definitions

AiiDA-VASP is under active development, check out the [changelog].


:::{note}
We are currently looking for additional developers. If you are interested, please open an issue on our repository on Github.
:::

Please consider to [open an issue] instead of sending an email to the AiiDA mailing list if the issue is related to this plugin.
Also, please consider that [AiiDA-VASP] is no substitute for not knowing how to use VASP. If in doubt at any point of the [VASP] parts of the
training material for [AiiDA-VASP] or when using [AiiDA-VASP], consult for instance [VASP lectures], [VASP tutorials], [VASP howtos],
[VASP tutorials using notebooks] or [VASP videos] or ask experienced [VASP] users.

Please accept that the development of this plugin is a community effort and any error, bug or missing functionality that might appear is the responsibility of the individual user. If you detect something that needs to be improved we highly encourage to [open an issue] or even better, [submit a pull request], where you try to fix the issue yourself. You can do this on our Github repository.



```{toctree}
:hidden: true
:maxdepth: 2

getting_started/index
tutorials/index
workflows/index
concepts/index
howto/index
developments/index
apidocs/index
./changelog.md
./faq.md
```


# Getting help

The development, maintenance and use of this plugin is considered a community effort. In order to facilitate for the community to contribute,
we have established a [space on Matrix] that users can use to communicate. We encourage users to help each other. In addition,
the development team is present in the space and users are free to ask.
First consult the documentation of both [AiiDA-VASP documentation] and [AiiDA documentation] and also consider that the developers are
not paid for this work. Please respect potential lead times in getting answers and be polite.




[aiida]: https://www.aiida.net
[aiida-vasp]: https://github.com/aiida-vasp/aiida-vasp
[changelog]: https://github.com/aiida-vasp/aiida-vasp/blob/develop/CHANGELOG.md
[conda]: https://docs.conda.io/en/latest/
[open an issue]: https://github.com/aiida-vasp/aiida-vasp/issues
[submit a pull request]: https://github.com/aiida-vasp/aiida-vasp/pull
[vasp]: https://www.vasp.at
[vasp howtos]: https://www.vasp.at/wiki/index.php/Category:Howto
[vasp lectures]: https://www.vasp.at/wiki/index.php/Lectures_and_presentations
[vasp tutorials]: https://www.vasp.at/wiki/index.php/Category:Tutorials
[vasp tutorials using notebooks]: https://www.vasp.at/tutorials/latest/
[vasp videos]: https://www.youtube.com/channel/UCBATkNZ7pkAXU9tx7GVhlaw
[vasp wiki]: https://cms.mpi.univie.ac.at/wiki/index.php
[aiida documentation]: http://aiida-core.readthedocs.io/en/latest/
[aiida-vasp documentation]: https://aiida-vasp.readthedocs.io/en/latest/
[space on matrix]: https://matrix.to/#/#aiida-vasp:matrix.org
