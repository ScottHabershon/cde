## Tutorial 3: Reaction path-finding {#pathfind}

This section gives example input files and instructions for a typical double-ended reaction path-finding calculation run through CDE.

This calculation will perform a simulated annealing (SA) optimization of the mechanism error-function, seeking out a reaction mechanism which connects the reactant graph (as calculated from the input reactant xyz file) to the product graph (as determined from the input product xyz file).

**Important: Note that the target reactant and product graphs are calculated directly from the input reactant and product structures. So, make sure that the bond-lengths in your reactants and products are correct so that the calculated connectivity matrices are also correct!**

The allowed moves (i.e. chemical reactions) which are used in the search for a mechanism connecting the reactant and product are provided in the move file.

In this example, we're going to look at the oxidation of carbon monoxide to carbon dioxide, occuring on a platinum cluster.

**The actual files to run this example can be found in the *~/cde/Tutorials/Tutorial_3/* directory.**

In this directory, you will find several input files:

	input
	start.xyz
	end.xyz
	moves.in


All of these input filenames are arbitrary; we could have called them alice, bob, clive and ralph and the program would still read them in correctly. Note that the extensions have no effect for input files read in by CDE, but they can be used to remind you what the different input files are.

Let's have a look at each input file individually:

### input

The *input* file for this GDS run looks like this:


    # Test input file

    optaftermove .true.
    calctype pathfind
    nimage 10
    startfile start.xyz
    endfile end.xyz
    ranseed 1

    # path optimization
    projforcetype 3
    nebmethod quickmin
    nebiter 1000
    cithresh 1d-4
    nebspring 0.1
    nebstep  10.0
    neboutfreq 5
    nebconv 1d-4

    # constraints
    dofconstraints 0
    atomconstraints 7
    1 2 3 4 5 6 7

    # PES input
    pestype  lammps
    pesfile   lmp.head
    pesopttype  lammps
    pesoptfile lmp.min
    pesexecutable '/Users/scott/code/lammps-22Aug18/src/lmp_serial'
    pesoptexecutable '/Users/scott/code/lammps-22Aug18/src/lmp_serial'

    # Graph-driven sampling (GDS) control.
    movefile moves.in
    gdsthresh 0.5
    gdsspring 0.02
    gdsrestspring 0.05
    nbstrength 0.03
    nbrange 2.2
    kradius 0.05
    ngdsiter 500
    ngdsrelax 7000
    gdsdtrelax 0.15

    # Chemical constraints.
    valencerange{
    Pt 2 12
    O 1 2
    C 1 4
    }

    reactiveatomtypes{
    Pt
    C
    O
    }

    reactiveatoms{
    all
    }

    reactivevalence{
    }

    fixedbonds{
    Pt Pt
    }

    allowedbonds{
    }

    # Pathfinder calculation setup parameters.
    nrxn 10
    nmcrxn 1000000
    mcrxntemp 1000000.0
    nmechmove 1
    graphfunctype 2


The input file in this case only contains (mostly) those options which are relevant to the pathfinder simulation; other input directives have been removed and will be ignored internally in CDE.

The parameter input blocks in the input file have already been covered in other tutorials for GDS and for CINEB (see, for example, @ref Annotated). In the case of a pathfinding calculation, the following are the new parameters which have an impact on the algorithm:

- **calctype**: Note that the *calctype* parameter must be set to *pathfind*.

- **nrxn**: This is the total number of reactions used in the reaction mechanisms being searched and generated. Note that a "null" reaction is automatically included as a possible graph-move (reaction) during the simulated-annealing search, so *nrxn* is the maximum number of active chemical reactions used in each proposed reaction mechanism.

- **nmcrxn**: This is the total number of simulated annealing iterations to search for; if a reaction is not found within this number of steps, the simulation simply stops without producing final structures.

- **mcrxntemp**: This is the *initial* temperature for the simulated-annealing search; note that this temperature is linearly-scaled to zero over the course of the *nmcrxn* iterations to drive the system to local minima with lower error functions. The actual initial temperature must be chosen to suit the graphfunctype (see below).

- **nmechmove**: This integer identifies how many of the reactions the simulated-annealing algorithm should attempt to change in *each* iteration. We haven't yet systematically investigated the impact of this number; 1 or 2 seems like sensible values.

- **graphfunctype**: This integer (0, 1, 2 or 3) defines the function-type which is optimized during the SA run. All of the graph-functions are zero when one has found a mechanisms which connects the reactants and products in *nrxn* steps or less. The difference arises due to the treatment of permutational invariance, as follows:
1. *graphfunctype 0* does not consider permutational invariance. The reactant and product structures must match in their atomic ordering, and one must also make a choice of which reactant atom ends up in which position in the prodoct structure. This is obviously not that useful for automatic reaction discovery.
2. *graphfunctype 1* accounts for permutational invariance by calculating a graph-error function based on eigenvalues of a matrix defined as \f$M_{ij} = (m_{i} m_{j})/d_{ij}\f$, where \f$m_{i}\f$ is the atomic mass of atom i and \f$d_{ij}\f$ is the shortest distance between any two atoms calculated in terms of bond connections.
3. *graphfunctype 2* accounts for permutational invariance by calculating a valence histogram for each pair of element types.
4. *graphfunctype 3* accounts for permutational invariance by calculating a graph-error function based on eigenvalues of a matrix which is identical to the connectivity matrix except with the atomic masses on the diagonal elements.


In the calculation setup for this example, we're allowing a maximum of 10 active reactions in the proposed reaction mechanisms, we're running for a maximum of 1000000 iterations, and we're starting at a temperature of 1,000,000 K. Note that this is not a "real" temperature, but is instead related to the graph error function type.

### start.xyz and end.xyz

The *start.xyz* and *end.xyz* files are standard XYZ format files (Cartesian coordinates in Angstroms).

In a reaction pathfinding simulation, *start.xyz* contains the coordinates of the **reactant** structure and *end.xyz* contains the coordinates of the **product** structure. These initial stuctures are used to calculate the reactant connectivity matrix and the target (product) connectivity matrix.


### moves.in

The moves.in file describes the allowed chemical reaction classes which are allowed to take place; these allowed moves are applied to configurations during a pathfinding calculation in order to generate new reaction end-points.

The format of the move-file is discussed in @ref moves.


### Running the calculation

With these input files, we are now able to run the pathfinding calculation. To do so, go into the *~/cde/examples/pathfind/* directory and type:

	cde.x input

As in the other tutorials, the above assumes that you have already made sure that CDE can be run by simply typing *cde.x*. See the [setup section] (@ref setup) for more details.

As the calculation runs, you will find lots of output files generated in the run directory.

The interesting output files are as follows:

- **mcopt.dat**: This file contains a running value of the total graph-error function as a function of the number of simulated annealing iterations.

- **input.log**: Once the simulated annealing calculation is complete, the *log* file will contain lots of useful information about the final reaction mechanism (if successful).

- **final_path.xyz**: Contains an *xyz* file with molecular snapshots of the intermediates generated along the final reaction-path.

- **final_path_rx_ZZZ.xyz**: Contains an approximation to the reaction-path for reaction-step *ZZZ* in the final mechanism.

The *xyz* files can obviously be plotted in VMD in order to visualize the reaction.

### Optimizing structures

If one sets the input parameter

        optaftermove .true.

then this means that each of the intemediate structures generated in the *final* reaction mechanism will be optimized, using the *pesopttype* and defined in the input file. In addition, the energy of each intermediate, as calculated according to *pestype*, will also be output to the *log* file.

Using *optaftermove .true.* is advised if one is intending to perform further calculations on the final reaction mechanism, such as NEB calculations for each reaction-path file (*final_path_rx_ZZZ.xyz*). In this case, the initial reaction-paths created for each reaction are generated such that they connect the optimized geomtries.

### Final output and next steps.

- Each of the reaction-path files ((*final_path_rx_ZZZ.xyz*) generated by a pathfinding could, in principle, be used as the starting-point of a CINEB refinement calculation.

- By comparing CINEB calculations for a number of different pathfinding simulation outputs, it should also be possible to identify the "most likely" reaction mechanism (for example, with most favourable thermodynamic and kinetic properties).
