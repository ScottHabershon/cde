##Tutorial 4: Another reaction path-finding example {#pathfind2}

This section gives a second example input file set for a typical double-ended reaction path-finding calculation run through CDE.

The calculation set-up is similar to Tutorial 3, but uses a different optimization function (type 4) instead. So, in this example, the target product is a single molecule, whereas the reactants contain a set of different molecular species.

This example can be run in the same way as tutorial 3, and similar input files are output.

The interesting output files are as follows:

- **mcopt.dat**: This file contains a running value of the total graph-error function as a function of the number of simulated annealing iterations.

- **input.log**: Once the simulated annealing calculation is complete, the *log* file will contain lots of useful information about the final reaction mechanism (if successful).

- **final_path.xyz**: Contains an *xyz* file with molecular snapshots of the intermediates generated along the final reaction-path.

- **final_path_rx_ZZZ.xyz**: Contains an approximation to the reaction-path for reaction-step *ZZZ* in the final mechanism.

- **adjusted_path_ZZZ.xyz**: These adjusted xyz files contain only those molecules which are relevant to the formation of the product of interest.

The *xyz* files can obviously be plotted in VMD in order to visualize the reaction.


