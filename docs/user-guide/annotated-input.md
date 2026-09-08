## Annotated input file description ## {#Annotated}

This page describes the allowed input parameters, and a description of what they do. Please note the following:

- Note that input-file keywords are case-sensitive - everything should be lower-case!
- In the input file, lines beginning with '#' are comments.
- Blank lines are ignored.
- Keywords can appear in any order.

## Main input parameters

These parameters are used to control input, number of images, and other general aspects of the calculation.

 *nimage 10*           ::  Number of images.

 *calctype  optpath*   :: Calculation type (*optpath* for path optimization, and *pathfind* for double-ended reaction-path finding calculations.)

 *startfile start.xyz* :: Start-point coordinate file (used when startfrompath = .FALSE.)

 *endfile end.xyz*     :: End-point coordinate file (used when startfrompath = .FALSE.)

 *pathfile path.xyz*   :: Initial path file (used when startfrompath = .TRUE.)

 *ranseed  13*     :: Random number seed (positive integer)

 *startfrompath   .FALSE.*   :: Flag indicating start from end-point files only (.FALSE.) or full path (.TRUE.)

 *pathinit  linear*   :: If startfrompath = .FALSE., pathinit indicates initial interpolation ( must be *linear* for now ...)

## Path optimization control

The following parameters are used to control the calculation when performing a path-optimization simulation (
**calctype optpath**), for example using climbing-image nudged elastic band.

 *pathoptmethod cineb*  :: Path-optimization method (currently cineb only)

 *nebmethod quickmin*   :: (CI)NEB optimization method (*fire* for fire algorithm (Bitzek *et al* Phys. Rev. Lett. 97, 170201 (2006)), *quickmin* for velocity verlet method (competitive), *steepest* for steepest-descents (SLOW))

 *nebiter 500*        :: Number of (CI)NEB optimization iterations.

 *cithresh 5d-4*      :: RMS force threshold at which climbing-image is activated.

 *nebspring 0.05*     :: (CI)NEB spring constant in au (Eh/bohr**2)

 *nebstep  10.1*    :: Step-size for steepest descent or quickmin (au)

 *neboutfreq 10*     :: Output frequency during (CI)NEB optimization.

 *optendsbefore .FALSE.* :: Flag indicating whether to optimize end-points BEFORE NEB refinement.

 *optendsduring .TRUE.* :: Flag indicating whether to optimize end-point SIMULTANEOUSLY with path.

 *nebrestrend .FALSE.* :: Flag indicating whether to apply graph-restrain potential during NEB (also used in GDS).

 *vsthresh 1d-3* :: RMS force threshold at which the variable spring strength is switched on (see Henkelman 00).

 *reconnect .true.* :: For 2->N-1 images, linearly interpolates the coordinates of beads 1 and N, AFTER a minimization has been performed

 *idppguess .true.* :: performs NEB under a Image-Dependent Pair Potential (IDPP) starting from a linear interpolation to generate a better initial guess for NEB. See Article Smidstrup *et al* J.
 Chem. Phys. 140, 214106 (2016)

 *projforcetype 2* :: Selects which type of NEB projected force to use (1: original, 2: Henkelman *et al*, J. Chem. Phys. 113, 9901 (2000), 3: Kolsbjerg *et al* J. Chem. Phys. 145, 94107 (2016) )

 *stripinactive .TRUE.* :: This is a logical flag which controls whether or not "inactive" molecules are removed from the initial path of a NEB calculation before NEB refinement.

 *nebconv 1d-3*    ::

 *nebmaxconv 5d-3*     ::


## Constraints control

This block controls the constrained atoms and degrees-of-freedom in GDS calculations and geometry-optimization calculations.

The *dofconstraints* input gives the number of fixed DOFs (in this case, 6), with following line in the input file listing the integer-values of the fixed DOFS. Here, counting starts at 1 for the *x*-coordinate of atom 1. The *z*-coordinate of atom 2 has the integer value of 6, and so on...

     dofconstraints 6
     1 2 3 5 6 9

The *atomconstraints* has the same format as *dofconstraints*, but the integers now refer to the number of fixed atoms, and the indices of the fixed atoms. In the following, 1 atom is fixed - atom number 13.

     atomconstraints 1
     13

Finally, the *alignedatoms* input allows one to define three atoms which will be used to position and orient any input molecular structure in space. This is very useful when performing NEB calculations as it helps remove overall translations and orientations of the molecules in the optimized string. Note that the NEB routines will automatically select 3 atoms to position and orient if you don't input anything as *alignedatoms*. In the following example, atoms 1, 2 and 3 are chosen to define the relative position and orientation.

     alignedatoms
     1 2 3

## PES control

This block controls potential energy evaluations during GDS and NEB, and also controls geometry-optimization.

 *pesfull .TRUE.*   :: If .TRUE., we perform PES evaluations and geometry optimizations for the full system in one go. If .FALSE., we separately calculate the energy for each independent molecule.

 *pestype orca*        :: Defines single-point PES evaluation code to use (*orca*, *dftb*, *null* or *lammps*)

 *pesfile  orca.head*  :: Defines the inputfile containing the header for pes single-point calculation.

 *pesopttype orca* :: Defines geometry optimization PES code being used (*orca*, *dftb*, *null*, *uff*)

 *pesoptfile orca.min* :: Defines the inputfile containing the header for geometry-optimization calculation.

 *pesexecutable orca* :: Defines the executable to use to for single-point PES evaluation (full path or alias)

 *pesoptexecutable orca* :: Defines the executable to use to for geometry optimization (full path or alias)


## Graph-restraint potential control

Controls GDS simulations.

**movefile moves.in** :: Identifies the movefile used to input allowed graph moves.

**gdsspring 0.025** :: Inter-bead spring constant for GRP (Eh/Bohr**2)

**gdsrestspring 0.1** :: Spring constant for harmonic graph enforcement terms (Eh/Bohr**2)

**nbstrength 0.02**  :: Repulsion strength for non-bonded atoms [V = nbstrength * exp(-r_{ij}**2 / (2 * nbrange**2))]

**nbrange 2.0** :: Range parameter for non-bonded atoms (see above).

**kradius 0.1**  :: Spring constant for non-interacting molecules in GDS simulations.

**ngdsrelax 2000** :: Number of steepest-descent optimization steps to minimize graph restraint potential

**gdsdtrelax 0.05** :: Step-size (au) for steepest-descent relaxation above.

**gdsoutfreq 10** :: Information output frequency during GDS simulation.

**nebrestrend .FALSE.** :: Flag indicating weather to apply to it duing the optimization of molecules during GDS, if optaftermove is .true. (also used in CINEB).

**optaftermove .TRUE.** :: Optimizes every molecule involved in the path formation - if the resulting optimised molecule does not conform to the graph,  the proposed path is rejected. If nebrestrend .TRUE., the graph constraining potential are also included during the optimization, and then switched off  during a second, in the hope to find a minima in agreement with the graphs (would recommend to set nebrestrend .TRUE. for this).

**valencerange** :: Defines the allowed valence ranges of each element as follows:

	valencerange{
	C 1 4
	O 1 2
	H 0 1
	Pt 0 3 fz
	}
The fz keyword applies only to frozen atoms, either by being defined in the fixedbonds{} input (see below), or by atomcontraints keyword. In the example,
only Pt-X bonds (X= whatever) which can be formed or broken during the GDS calculation are counted into the 'valence'. If a cluster of Pt atoms are not
moving during a GDS calculation, and they are bonded to each other (to varying degrees), the above keyword will allow for up to 3 more bonds, not counting the
Pt-Pt bonds, to be formed on the Pt atoms.

**reactiveatomtypes** :: Defines which types of elements are allowed to react, as follows:

	reactiveatomtypes{
	C
	O
	H
	}

**reactiveatoms** :: Defines which atom numbers are allowed to react. These can be given as an inclusive range (*range XX YY*), or by individual id number (*id XX*).

	reactiveatoms{
	range 1 4
	id 7
	}

**reactivevalence** :: Defines the valence range one element can have with another

	reactivevalence{
	Fe Fe 6 8
        H C 0 1
        H O 0 1
	}

In the above, *Fe* must be bonded between 6 to 8 other Fe atoms to be an acceptable molecule. Similarly, the second line indicates that Hydrogen can only
be bonded up to 1 carbon atom.

**fixedbonds** :: Defines fixed bonds, which are not allowed to change during a chemical reaction (although they can vibrate, translate, etc.)

	fixedbonds{
	C O
	3 4
	}

In the above example, *all* C-O bond orders are fixed at whatever they are in the starting structure. So, if a C-O bond is present in the starting structure, it cannot change during the GDS simulation. Similarly, the bond between atoms *3* and *4* is also fixed at the starting structure value.


**essentialmoveatoms** :: Defines a list of atoms of which AT LEAST ONE OF THEM MUST be included in any possible graph moves. The formatting options are the same as the *reactiveatoms* block described above, for example:

	essentialmoveatoms{
	range 1 2
	}

In the above example, the atoms in the range 1-2 must be included amogst the moves chosen.

**essentialatoms** :: Defines atoms of which AT LEAST ONE OF THEM MUST be included in molecules that are involved in the reactions. The formatting options are the same as the *reactiveatoms* block described above, for example:

	essentialatoms{
	id 5
	id 7
	}

In the above example, the atoms 5 and 7 must be in any of the molecules that are involved in the reaction (from chemical species image 1 changing to nimage)

**allowedbonds** :: Allows definition of bond-number constraints in the structures generated by graph moves. The format is as given in the following example:

	allowedbonds{
	O O 1
	}

Here, the constraint indicates that an oxygen atom cannot be bonded to more than 1 other oxygen atom.


**forbidgraphs** :: This is a logical flag which indicates whether or not to use the forbidden graph pattern file. If .TRUE., then a GDS simulation will not allow any graph-patterns to form which are input into the forbidfile defined below. By adding entries to the forbid file, users can force GDS to stop generating user-defined bonding patterns.


**forbidfile forbid.in** :: Identifies the file (here, forbid.in) which contains the library of forbidden bonding patterns which are used if forbidgraphs = .TRUE.

## Path-finding controls

The following parameters control the double-ended reaction-path finding algorithm.

**nrxn** :: Total number of elementary steps allowed in mechanism.

**nmcrxn** :: Maximum number of Monte Carlo moves to try during simulated annealing optimization.

**nmechmove** :: Maximum number of mechanism updates in each MC step.

**mcrxntemp** :: Initial temperature (in K) for simulated annealing optimization (decreases linearly during run).

**graphfunctype** :: Determines the optimization function for the path-finding calculation. This should be an integer: (1) Standard element-wise comparison, (2) Eigenvalue comparison, (3) Histogram comparison, (4) Eigenvalue comparison for single molecule. *You should preferably use either (2) or (4)!*

**minmolcharge** :: Minimum allowed molecular charge (if electron transfer moves are allowed in the movefile).

**maxmolcharge** :: As above, but maximum value.

**nchargemol** :: Maximum number of molecules which can be charged.

**maxstepcharge** :: Maximum number of reation steps which can involve charge changes.

**maxtotalcharge** :: Maximum total charge in a given reaction step.
 
