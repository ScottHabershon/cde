## Forbid-file example {#forbid}

The forbid-file is used in GDS to define which graph-patterns cannot be generated as new reaction-path end-points during each GDS graph-change step.

A typical forbid-file (usually called *forbid.in*, although can be called anything) looks like this:

	forbid
	natom 3
	-
	0 1 0
	1 0 1
	0 1 0
	-
	labels C O O
	
	forbid
	natom 3
	-
	0 1 0
	1 0 1
	0 1 0
	-
	labels C O C


Note that the format is very similar to that of the moves files. 

Each forbidden graph-patter is defined in a block, which looks like this.

	forbid	
	natom 3
	-
	0 1 0
	1 0 1
        0 1 0
	-
	labels C O O

In the example above, the bonding pattern specified by the graph (e.g. bond between atoms 1 and 2, bond between atoms 2 and 3) is NOT allowed to be generated for the atom-labels C,O,O. In othert words, this pattern would stop C-O-O, and any molecules containing this pattern, from forming. 

Note that, in contrast to the move-file, the atom labels **must** be defined in the forbid-file.
