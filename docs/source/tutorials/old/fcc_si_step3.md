(tutorial-fcc-si-step3)=

# 3. Your first workchain

Here we will continue where the {ref}`previous tutorial<tutorial_fcc_si_step2>` left of.
In particular we will now create a dedicated workchain that calls the {ref}`vasp_workchain`
and thus eventually executes VASP for each volume. The workchain makes sure the total energies
are extracted and stored. In the end we will even use some Python
functions to locate the minimum and store that as well. We hope that you after this tutorial will
finally start to appreciate the convenience of [AiiDA-VASP] and maybe also start to get a feel
for what is possible. Even better, when performing calculations like this, especially when the
workchain stack becomes much more complex than showed in this example, the workflow becomes
standardized and versioned, so it is easier to share and collaborate. Also, reproducibility is finally
possible in practice. Let us start to develop the workchain. As before we assume you have
completed the steps of the {ref}`previous tutorial<tutorial_fcc_si_step2>`.

In order to make it as easy as possible, we have created a cookie cutter recipe for creating
[AiiDA-VASP] workchains. This is located in `cookiecutters/workchain` on the root of the
repository. If you have not cloned it, you can download only the cookiecutter folder using:

```
$ svn checkout https://github.com/aiida-vasp/aiida-vasp/trunk/cookiecutters
```

What we want to do is to create a workchain that calculates some equation of state (EOS), in this
case over volume. For simplicity we will call the workchain `EosWorkChain`. It is always recommended to
dissect the operations we will perform into parts before writing the workchain. In this way, it
is easier to reuse larger parts of the workchain at a later point. Also, it is easier
to understand what is going on from just inspecting the definitions of the inputs, steps and outputs
of the workchain. Please, have a look at {ref}`workchains` and references therein before continuing.
Also, consider to look into the [tutorial on writing workflows] in the [AiiDA documentation] before
continuing.

01. First, let us Make sure you have `cookiecutter` installed:

    ```
    $ pip install cookiecutter
    ```

02. Now we simply execute the `cookiecutter` program using the `workchains` folder as the template:

    ```
    $ cookiecutter cookiecutters/workchains
    workchain_name [testchain]: eos
    workchain_docstring [This workchain is a testchain]: This workchain will accept a dictionary of structures and extract the total energies for each structure. The data is saved and the energy minimum is calculated and stored.
    next_workchain_to_call [vasp.vasp]:
    ```

03. A new folder `eos` was created which contains the `eos.py` workchain file. Inspect it.
    Maybe the most important part of a workchain is the `define` function, which, as its name
    implies, defines the workchain. For this particular case it looks like:

    ```{literalinclude} ../../../tutorials/eos_base.py
    :end-before: def initialize
    :start-after: classmethod
    ```

    We have defined some inputs:

    - All inputs defined in the `_next_workchain`, which in this case is the {ref}`vasp_workchain`
      are also available for this workchain by using `expose_inputs`. You can of course expose
      inputs from other workchains as long as [AiiDA] can find them. In this way, we do not have
      to repeat the definition of more basic workchains that can easily be repeated, like the {ref}`vasp_workchain`.
    - A namespace called `structures`, which is defined as a dictionary. Since we here do not know how many structures will be supplied,
      we have defined this input as a dynamic namespace. However, we know the dictionary should only
      contain {py:class}`StructureData<aiida.orm.nodes.data.structure.StructureData>` values.

    Then we have defined some outputs:

    - We attach all outputs from `_next_workchain`, i.e. the {ref}`vasp_workchain` to this
      workchain (basically a link is created to avoid storing things twice). Similar to
      the inputs, `expose_outputs` can be leveraged.

    And some generic exit codes.

    Finally, there is the `outline` section which tells the workchain which functions to run. In
    this case we have:

    - A function that initializes what we need, `initialize`.
    - A while loop that iterates over all the structures.
    - A `init_next_workchain` function that initializes what we need for the next workchain
      at the current iteration (here we make sure the specific structure is set).
    - A `run_next_workchain` which basically executes {ref}`vasp_workchain`.
    - A `verify_next_workchain` which verifies that there is a valid {ref}`vasp_workchain`
      and inherits any exit codes present from {ref}`vasp_workchain`.
    - A `extract_results` which gets what we need from the outputs.
    - A `finalize` which stores the results.

    Notice that the `cookiecutter` gave all this automatically, which is rather useful as
    a starting point to develop new workchains.

04. The workchain is however not yet ready. We need to extract the total energies
    and specify the general workchain further. Please make the following changes
    to the generated workchain:

    ```{literalinclude} ../../../tutorials/eos.py
    :diff: ../../../tutorials/eos_base.py
    ```

    and save it as `eos.py`. Or you could also download it:

    ```
    $ wget https://github.com/aiida-vasp/aiida-vasp/raw/master/tutorials/eos.py
    ```

    The majority of changes were related to being more specific, except two things:

    - The necessity of decorating the function that generates the output array containing the
      volume and total energies in a `calcfunction`. This is to preserve data provenance,
      otherwise we would not have known how the data was collected from each of the underlying
      workchains.
    - The inclusion of a `calcfunction` that interpolates the calculated data to find
      a better estimate of the volume at the energy minima. The example uses a
      cubic fit, which is certainly not very physical and should not be used in production.
      It is only to show how simply it is to leverage the power of Python, NumPy, SciPy and [AiiDA].
      This was decorated with a `calcfunction` in order to make sure [AiiDA] can
      honor data provenance.

05. Next, download the launch script that is tailored to launch the workchain we have now developed:

    ```
    $ wget https://github.com/aiida-vasp/aiida-vasp/raw/master/tutorials/run_fcc_si_workchain.py
    ```

06. Change the `options` and `code_string` as you did in {ref}`previously<tutorial_fcc_si_step1>`.

07. Now we need to make sure the daemon can pick up the workchain. We can do this by
    making sure the daemon sees the directory where `eos.py` and `run_fcc_si_workchain.py` is
    located. The simplest approach is to add the following, to your virtual environment `activate`
    script (assuming you do not use Conda):

    ```
    $ echo "export PYTHONPATH=$PYTHONPATH:<yourdirectory>" >> ~/env/aiida-vasp/bin/activate
    ```

    assuming `<yourdirectory>` is the directory containing the `eos.py` and
    `run_fcc_si_workchain.py` files. The location of the `activate` is assumed from the
    previous steps in the tutorial. If you use Conda, please do:

    ```
    $ echo "export PYTHONPATH=$PYTHONPATH:<yourdirectory>" >> $CONDA_PREFIX/etc/conda/activate.d/env_vars.sh
    ```

    :::{warning}
    Make sure you have (re)activated your [AiiDA] virtual environment and
    (re)started the [AiiDA] daemon before continuing. Otherwise, the updated `PYTHONPATH` variable
    will not be picked up by the daemon worker. You may want to double check by running `echo $PYTHONPATH`
    in the shell before starting the daemon with `verdi daemon start`.
    :::

08. Submit the workchain by running the call script:

    ```
    $ python run_fcc_si_workchain.py
    ```

09. After a while we check the status:

    ```
    $ verdi process list -a
      PK  Created    Process label         Process State      Process status
    ----  ---------  --------------------  -----------------  -----------------------------------------------------------
    1721  22m ago    EosWorkChain          ⏹ Finished [0]
    1722  22m ago    VaspWorkChain         ⏹ Finished [0]
    1724  22m ago    VaspCalculation       ⏹ Finished [0]
    1728  19m ago    VaspWorkChain         ⏹ Finished [0]
    1730  19m ago    VaspCalculation       ⏹ Finished [0]
    1734  17m ago    VaspWorkChain         ⏹ Finished [0]
    1736  17m ago    VaspCalculation       ⏹ Finished [0]
    1740  15m ago    VaspWorkChain         ⏹ Finished [0]
    1742  15m ago    VaspCalculation       ⏹ Finished [0]
    1746  13m ago    VaspWorkChain         ⏹ Finished [0]
    1748  13m ago    VaspCalculation       ⏹ Finished [0]
    1752  10m ago    VaspWorkChain         ⏹ Finished [0]
    1754  10m ago    VaspCalculation       ⏹ Finished [0]
    1758  8m ago     VaspWorkChain         ⏹ Finished [0]
    1760  8m ago     VaspCalculation       ⏹ Finished [0]
    1764  6m ago     VaspWorkChain         ⏹ Finished [0]
    1766  6m ago     VaspCalculation       ⏹ Finished [0]
    1770  4m ago     VaspWorkChain         ⏹ Finished [0]
    1772  3m ago     VaspCalculation       ⏹ Finished [0]
    1777  1m ago     store_total_energies  ⏹ Finished [0]
    1779  1m ago     locate_minimum        ⏹ Finished [0]

    Total results: 165

    Report: last time an entry changed state: 1m ago (at 17:18:23 on 2022-12-21)
    Report: Using 2% of the available daemon worker slots.
    ```

    As you can see, seven {ref}`vasp_workchain` and {ref}`vasp_calculation` were executed,
    one for each supplied volume.
    Also, there is a separate entry for the storage of the total energies, which also performs
    a sort. The location of the minima is also listed as a separate process as we decorated that
    function with a `calcfunction` decorator.

10. Let us have a look at the output of `EosWorkChain`:

    ```
    $ verdi process show 1721
    Property     Value
    -----------  ------------------------------------
    type         EosWorkChain
    state        Finished [0]
    pk           1721
    uuid         69e9a920-c783-4f61-ac73-349b6e19059d
    label
    description
    ctime        2022-12-21 16:58:07.145071+01:00
    mtime        2022-12-21 17:18:23.379756+01:00

    Inputs              PK    Type
    ------------------  ----  -------------
    structures
        silicon_at_3_5  1703  StructureData
        silicon_at_3_6  1704  StructureData
        silicon_at_3_7  1705  StructureData
        silicon_at_3_8  1706  StructureData
        silicon_at_3_9  1707  StructureData
        silicon_at_4_0  1708  StructureData
        silicon_at_4_1  1709  StructureData
        silicon_at_4_2  1710  StructureData
        silicon_at_4_3  1711  StructureData
    clean_workdir       1720  Bool
    code                818   InstalledCode
    kpoints             1712  KpointsData
    max_iterations      1719  Int
    options             1716  Dict
    parameters          1713  Dict
    potential_family    1714  Str
    potential_mapping   1715  Dict
    settings            1717  Dict
    verbose             1718  Bool

    Outputs        PK  Type
    -----------  ----  ---------
    eos          1778  ArrayData
    eos_minimum  1780  Dict

    Called      PK  Type
    --------  ----  --------------------
    CALL      1722  VaspWorkChain
    CALL      1728  VaspWorkChain
    CALL      1734  VaspWorkChain
    CALL      1740  VaspWorkChain
    CALL      1746  VaspWorkChain
    CALL      1752  VaspWorkChain
    CALL      1758  VaspWorkChain
    CALL      1764  VaspWorkChain
    CALL      1770  VaspWorkChain
    CALL      1777  store_total_energies
    CALL      1779  locate_minimum

    Log messages
    ---------------------------------------------
    There are 9 log messages for this calculation
    Run 'verdi process report 1721' to see them
    ```

11. Inspect the total energies versus volume:

    ```
     $ verdi data core.array show 1777
     {
         "eos": [
         [
            10.71875,
            -4.42341939
         ],
         [
            11.664,
            -4.66006381
         ],
         [
            12.66325,
            -4.79595554
         ],
         [
            13.718,
            -4.86303429
         ],
         [
            14.82975,
            -4.87588353
         ],
         [
            16.0,
            -4.8481407
         ],
         [
            17.23025,
            -4.78451926
         ],
         [
            18.522,
            -4.69228837
         ],
         [
            19.87675,
            -4.58122058
         ]
         ]
    }
    ```

12. And the located minimum:

    ```
    $ verdi data core.dict show 1780
    {
        "energy": -4.8769540208841,
        "volume": 14.559367617229
    }
    ```

That concludes this tutorial. We hope at this point you have now realized
that [AiiDA-VASP] seems somewhat useful and that you would like to continue to
learn more, maybe even start to write your own {ref}`workflows` or {ref}`workchains`.
You might have noticed when running this workflow that the each volume was running sequentially
and was a bit concerned about that being not so efficient as there is no data sharing between
the different volume runs. And indeed you are right. The next tutorial will show how this can be
addressed.

[aiida]: https://www.aiida.net
[aiida documentation]: https://aiida.readthedocs.io/projects/aiida-core/en/latest/index.html
[aiida-vasp]: https://github.com/aiida-vasp/aiida-vasp
[fcc si]: https://cms.mpi.univie.ac.at/wiki/index.php/Fcc_Si
[gnuplot]: http://gnuplot.info/
[tutorial on writing workflows]: https://aiida.readthedocs.io/projects/aiida-core/en/latest/intro/tutorial.html#workflows
[vasp]: https://www.vasp.at
