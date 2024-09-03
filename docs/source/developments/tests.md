(running-tests)=

# Running tests

[AiiDA-VASP] utilizes [pytest] which makes it easy to run and construct new tests. We also utilize the fixture
managers that are present, see the [AiiDA documentation].

In order to run all tests execute, in the root folder:

```
$ pytest
```

This is a useful test to see if the plugin works and that everything is installed correctly.

:::{note}
For contributors, one principle rule is in place: always try to construct tests when submitting a PR.
:::

:::{note}
We also utilize `tox` in order to make the test environment fully isolated. In order to run tests
in `tox`, please issue `tox` in the root folder of [AiiDA-VASP].
:::

[aiida documentation]: https://aiida.readthedocs.io/projects/aiida-core/en/latest/index.html
[aiida-vasp]: https://github.com/aiida-vasp/aiida-vasp
[pytest]: https://docs.pytest.org/en/latest/
[vasp]: https://www.vasp.at/
