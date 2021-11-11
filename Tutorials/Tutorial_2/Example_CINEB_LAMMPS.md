### Tutorial 2: CINEB using LAMMPS {#Tutorial2}

This tutorial gives a set of example input files and instructions for a CINEB calculation which uses the ReaxFF forcefield through LAMMPS.

In this case, we are going to run a CINEB optimization of a reaction path representing association of water H2O onto a Pt atom which also has a pre-bound CO molecule. We will use the ReaxFF reactive force-field to model the potential energy surface, accessed through the LAMMPS molecular dynamics code.

**The actual files to run this example can be found in the *~/cde/examples/Tutorial_2 directory.**

In this directory, you will find several input files:

      input
      path.xyz
      CHOPtNiX.ff
      lmp.head
      lmp.min
      lmp_control

These files are all required in order to the run the CINEB/LAMMPS calculation defined in the input file.

In general, the *input* is similar to that given in *Tutorial 1* (@ref Tutorial1). The key difference is that the setup of the PES evaluations is different - in *Tutorial 1* we used *ORCA*, but here we use ReaxFF via LAMMPS. As a result, the *PES* input block of the input file looks like the following in this case:

    pestype  lammps
    pesfile   lmp.head
    pesopttype  lammps
    pesoptfile lmp.min
    pesexecutable '/Users/scott/code/lammps-22Aug18/src/lmp_serial'
    pesoptexecutable '/Users/scott/code/lammps-22Aug18/src/lmp_serial'

There are a few important things to note about running calculations with ReaxFF/LAMMPS:

- Here, the *pestype* and *pestoptype* have both been set to *lammps*, and the paths to the relevant LAMMPS executable have also been specified.

- **To get this example running, you MUST modify the executables to point to the LAMMPS executable on your own system.**

- Another difference with *Tutorial 1* is that the PES template files are different. In this case, the template files *lmp.head* and *lmp.min* are LAMMPS template files. More information can be found in the *PES templates* (@ref Templates) section of the documentation.

- Finally, as also noted in the *PES templates* (@ref Templates) section, the LAMMPS code requires two additional files to be present in the run directory. These are:
(1) A *.ff file, containing the ReaxFF parameters which should be used. In our case, this is called *CHOPtNiX.ff*.
(2) A *lmp_control* file, which contains additional control parameters used by LAMMPS!
