## Hints and tips {#hints}

### Hints on input files.

- The main CDE input file is CASE SENSITIVE. *Everything should be lower-case, except for chemical element symbols*.
- If you are running a calculation which does not use some of the input parameters, they will be read in but then ignored. However, their formatting
still needs to be correct.

### Hints on bonding definitions.

- In CDE, bonds between atoms are based on a cutoff distance given by *R = (CovRad(i) + CovRad(j)) * bondingsf*. The factors CovRad(:) are the covalent radii of each species, and can be found in *constants.f90*. The factor **bondingsf** is also found in constants.f90 - it is a scale-factor, usually with a value of 1.1.
- The CovRad(:) and bondingsf are defined as Fortran constants - if you want to change them, you need to recompile the code.
- The Covrad(:) and bondingsf values do not need to be extremely accurate - they are used to define what is bonded to what, so as long as they adequately capture typical bonds, they should be fine.

### Hints on external executables

- The input parameters **pesexecutable** and **pesoptexecutable** define the external executables which are run when evaluating energy or performing geometry optimizations.
- The CDE code uses the Fortran implicit *EXECUTE_COMMAND_LINE* to run these executables.
- It is often easiest to define an alias for executables like ORCA and DFTB+. For example, if the ORCA executable is in a directory */user/test/code/orca/bin*, then setting an alias so that the command *orca* actually executes the command */user/test/code/orca/bin/orca* is a useful way of simplifying the input file.
- If the executables can't be found when run by *EXECUTE_COMMAND_LINE*, the code will just give up and crash!
- If you want to use **pesopttype UFF** (that is, the Universal Force-Field), you need to install *OpenBabel*.
- When using **pesopttype UFF**, the executables *babel* and *obminimize* are called directly. As a result, both of these executables need to be in your PATH variable before running CDE.


### Common problems

- In several types of calculations (e.g. double-ended reaction-path finding, NEB), the first step that the CDE code performs is to evaluate the connectivity matrix of the input molecular structures. Depending on the calculation type, this connectivity matrix is then used to define the graph-restraining potential (GRP) for structure optimization and checking. However, if the connectivity matrix calculated using the input structure does not correspond to the connectivity matrix that you wanted (for example, due to bond-lengths being too long in the input sutrcture), then you might end up with some odd results. **In short, TAKE CARE with your input structures!**
