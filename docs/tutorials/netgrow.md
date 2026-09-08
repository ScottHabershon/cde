## Tutorial 5: Single-ended reaction generation (EXPERIMENTAL) {#netgrow}

This section gives example input files and instructions for a typical single-ended mechanism-finding calculation run through CDE.

This calculation will apply a number of randomised graph moves to propagate out a reaction from a given start point. Moves get applied at random to random steps of the reaction (this is not necessarily desired functionality, hence why this method is experimental). Note that while move probabilities should be set within moves.in, these have no effect in this routine and all moves have equal likelihood of occurring.

**Important: Note that the target reactant graph is calculated directly from the input reactant structure, so make sure that the bond lengths in your reactants are correct so that the calculated connectivity matrices are also correct!**

The allowed moves (i.e. chemical reactions) which are used in the search for a mechanism are provided in the move file.

In this example, we will use xTB to perform our electronic structure calculations.

**The actual files to run this example can be found in the *~/cde/Tutorials/Tutorial_5/* directory.**

In this directory, you will find several input files:

	input
	start.xyz
	moves.in
    xtb.head

All of these input filenames are arbitrary; we could have called them alice, bob, clive and ralph and the program would still read them in correctly. Note that the extensions have no effect for input files read in by CDE, but they can be used to remind you what the different input files are.

Let's have a look at each input file individually:

### input

The *input* file for this GDS run looks like this:

    calctype netgrow 
    optaftermove .false.
    pesfull .false.         ! Required for radicals and charges...
    nimage 8                ! Number of images in each generated reaction-step 
    startfile start.xyz
    endfile end.xyz
    idppguess .true.

    pathoptmethod cineb
    nebmethod quickmin 
    nebiter 500
    cithresh 5d-4
    nebspring 0.05
    nebstep  10.0 
    neboutfreq 5
    nebrestrend .false.

    dofconstraints 0
    atomconstraints 0

    pestype xtb
    pesfile xtb.head
    pesopttype xtb 
    pesoptfile xtb.head
    pesexecutable xtb --iterations 1000 --grad
    pesoptexecutable xtb --iterations 1000 --opt normal --grad

    movefile moves.in
    gdsthresh 0.5 
    gdsspring 0.05
    gdsrestspring 0.05
    nbstrength 0.04
    nbrange 2.5
    kradius 0.05
    ngdsrelax 10000
    gdsdtrelax 0.1
    gdsoutfreq 10

    valencerange{
    C 1 4
    H 0 1
    }

    reactiveatomtypes{  
    C 
    H 
    }

    reactiveatoms{  
    all
    }

    reactivevalence{
    }

    fixedbonds{
    }

    allowedbonds{
    } 

    nrxn 10
    nmcrxn 5000
    graphfunctype 4
    ranseed 44

The input file in this case only contains (mostly) those options which are relevant to the netgrow simulation; other input directives have been removed and will be ignored internally in CDE.

The parameter input blocks in the input file have already been covered in other tutorials for GDS and for CINEB (see, for example, @ref Annotated). In the case of a netgrow calculation, the following are the new parameters which have an impact on the algorithm:

- **calctype**: Note that the *calctype* parameter must be set to *netgrow*.

- **nrxn**: This is the total number of reactions used in the reaction mechanisms being searched and generated. Note that a "null" reaction is automatically included as a possible graph-move (reaction) during the search, so *nrxn* is the maximum number of active chemical reactions used in each proposed reaction mechanism.

- **nmcrxn**: This is the total number of mechanism updates to perform. Mechanism updates are performed on a random mechanism step each time.

### start.xyz

The *start.xyz* file is a standard XYZ format file (Cartesian coordinates in Angstroms).

In a netgrow simulation, *start.xyz* contains the coordinates of the **reactant** structure.

### moves.in

The moves.in file describes the allowed chemical reaction classes which are allowed to take place; these allowed moves are applied to configurations during a netgrow calculation in order to generate new reaction end-points.

The format of the move-file is discussed in @ref moves.

### xtb.head

The xtb.head file is interntionally blank as xtb does not require a configuration file. Instead, it is set up to automatically run calculations on an xyz file which CDE generates (xtbin.xyz), and any calculation-dependent arguments can be passed directly to the pesexecutable and pesoptexecutable input arguments. However, it must be present as a dummy file since CDE expects such a file to exist for each calculation.

### Running the calculation

With these input files, we are now able to run the pathfinding calculation. To do so, go into the *~/cde/Tutorials/Tutorial_5/* directory and type:

	cde.x input

As in the other tutorials, the above assumes that you have already made sure that CDE can be run by simply typing *cde.x*. See the [setup section] (@ref setup) for more details.

As the calculation runs, you will find lots of output files generated in the run directory.

### Optimizing structures

If one sets the input parameter

        optaftermove .true.

then this means that each of the intemediate structures generated in the *final* reaction mechanism will be optimized, using the *pesopttype* and defined in the input file. In addition, the energy of each intermediate, as calculated according to *pestype*, will also be output to the *log* file.

Using *optaftermove .true.* is not necessarily recommended within this routine, as it has no way of telling if a given mechanism step's graph has been invalidated by the optimisation and as such can lead to internally inconsitent mechanisms.