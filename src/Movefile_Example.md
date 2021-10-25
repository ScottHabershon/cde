## Movefile example {#moves}

The move-file is used in GDS to define which graph-moves can be used to generate new reaction-path end-points during each GDS graph-change step.

A typical move-file (usually called *moves.in*, although can be called anything) looks like this:


    move
    natom 2
    -
    0 1
    1 0
    -
    0 0
    0 0
    -
    labels * *
    prob 0.4

    move
    natom 2
    -
    0 0
    0 0
    -
    0 1
    1 0
    -
    labels * *
    prob 0.4

    move
    natom 3
    -
    0 1 0
    1 0 1
    0 1 0
    -
    0 0 1
    0 0 1
    1 1 0
    -
    labels * * *
    prob 0.1

    move
    natom 3
    -
    0 1 0
    1 0 1
    0 1 0
    -
    0 0 1
    0 0 0
    1 0 0
    -
    labels * * *
    prob 0.1


Each move is defined in a block, which looks like this.

	move
	natom 2
	-
	0 1
	1 0
	-
	0 0
	0 0
	-
	labels * *
	prob 0.4

In the example above, the proposed graphmove involves two atoms; these are randomly selected in the code. The first 2x2 matrix is the required bonding pattern of the two atoms; in this case, the moves requires that the two atoms are bonded (*i.e.* there is a 1, indicating bonding, on the off-diagonal element). The second 2x2 matrix indicates the target final graph; in thise case, there is a *0* on the off-diagonal, indicating that there is no bond between the two atoms in the final structure. Overall, therefore, this is a simple bond-breaking reaction.

The *labels* line indicates the required atom labels for this move; in this case, we have indicated " * ", which means that **any** atom can participate. An alternative would be

	labels C O

which would imply that only C and O can be involved for this graph move.

The *prob* value gives the relative probability of this move taking place during a given GDS run; not that the probabilities of all moves are all normalised in the code (*i.e.* they are re-scaled so that they sum to 1), so they do not have to add to 1 in the movefile.

A second example, this time involving 3 atoms, is as follows:


	move
	natom 3
	-
	0 1 0
	1 0 1
	0 1 0
	-
	0 0 1
	0 0 0
	1 0 0
	-
	labels * * *
	prob 0.1

The above reaction is A-B-C --> A-C + B.
