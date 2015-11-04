# Open Questions about VASP Workflow

## Inputfiles

### incar
Which tags to automatically determine? From which other tags? From what other information?
Which tags influence what other params?

Example:
* nbands -> num_total_mpiprocs
* choice of potential -> other stuff
* type of material -> ismear
* type of calculation -> ismear, icharg
...

### potcar
How to determine which potpaw symbol to use (for example InAs -> In_d, not In)?
What info to retrieve and store separately?

### poscar / contcar
What causes the differences between pymatgen::structure and aiida::structure read from the same CIF?
Are the differences significant for our purposes?
Are there any benefits to use / be compatible with aiida::structure?
What to retrieve from contcar?

### kpoints
Why different kpoint labels for publication / pymatgen / hand-made kpoints file?
Are the differences significant?
Which way to support?
retrieve all, right?

### output files
What to retrieve from where???
