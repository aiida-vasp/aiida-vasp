(running-tests)=

# Tests

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

## Workchain tests using mock-vasp

The plugin is shipped with a *mock* VASP script to simulate running VASP without actually running it. This allows us to test the workchains without having to run VASP.

The `mock-vasp` analyse the input and compute a hash, then searches a registry folder containing
pre-computed calculations results for a match. If a match is found, it copies the output files to the folder of the calculation.
Otherwise, if `MOCK_VASP_VASP_CMD` is set, it runs VASP with the given command and upload
the output to the registry.

To generate the test data for `pytest` when developing the plugin, set the following environment variables:

- `MOCK_VASP_POTCAR_PATH`: path to the directory containing the POTCAR files
- `MOCK_VASP_VASP_CMD`: command to run VASP
- `MOCK_VASP_REG_BASE`: Base of the registry for the mock VASP executable, defaults to the `<root>/tests/test_data` folder.
- `MOCK_VASP_UPLOAD_PREFIX`: A prefix to use when creating folders of containing the registered VASP calculations.

:::{note}
`mock-vasp` does not analyse nor copy the POTCAR files as they are in most cases licensed.
The POTCAR files are simply ignored when computing the hash of calculation inputs.
:::

To test a new workchain and record the calculation data, set `MOCK_VASP_VASP_CMD` to your real vasp executable (optionally, with the MPI launcher and arguments).
Set `MOCK_VASP_POTCAR_PATH` to the path containing the POTCAR files, this will upload the POTCAR and create a corresponding `PBE.54` potcar family.
The `MOCK_VASP_UPLOAD_PREFIX` should be set when constructing the calculation input, as the
custom scheduler comments.
This allows one to known which calculation data is for which test case.
Finally, run the tests as usual as `pytest`, and rerun without setting the environmental variables above to make sure they pass without using the actual VASP executable.


[aiida documentation]: https://aiida.readthedocs.io/projects/aiida-core/en/latest/index.html
[aiida-vasp]: https://github.com/aiida-vasp/aiida-vasp
[pytest]: https://docs.pytest.org/en/latest/
[vasp]: https://www.vasp.at/
