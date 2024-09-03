(writing_tutorials)=
# Writing tutorials

The tutorials are written as jupyter notebook in the markdown format.
The VASP calculations mocked using `mock-vasp` which is shipped with this package.

To write a new tutorial, one can first draft it in a jupyter notebook and then convert it to markdown format using the `jupytext` package:

```bash
jupytext my_tutorial.ipynb --to md
```

This will create a new file `my_tutorial.md` which can be edited in a text editor.

To include the mocked calculation in the tutorial (to run on machines without access to VASP)
set the `MOCK_VASP_VASP_CMD` to the command to run the VASP executable during the execution of the notebook (e.g. `make html` for building the documentation).

## Timing for the notebook execution

```{nb-exec-table}
```
