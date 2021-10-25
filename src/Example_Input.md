## Example input file for CDE ## {#Example}

This is an example input file for CDE.


     # Example input file - note that all possible options are included here, but are NOT required in the input file if not being used.

     nimage 8
     temperature 10
     calctype pathfind 
     startfile start.xyz
     endfile end.xyz
     pathfile path.xyz
     ranseed  1
     startfrompath   .FALSE.
     pathinit  linear

     # path optimization

     pathoptmethod cineb
     nebmethod quickmin
     nebiter 500
     cithresh 5d-4
     nebspring 0.05
     nebstep  10.1
     neboutfreq 10
     stripinactive .TRUE.
     optendsduring .FALSE.
     optendsbefore .TRUE.
     reconnect .TRUE.
     idppguess .TRUE.
     projforcetype 3
     vsthresh 1d-3
     nebconv 1d-3
     nebmaxconv 5d-3

     #constraints

     dofconstraints 6
     1 2 3 5 6 9

     atomconstraints 0

     alignedatoms
     1 2 3

     #PES input

     pesfull .TRUE.
     pestype  null
     pesfile   orca.head
     pesopttype  null
     pesoptfile orca.min
     pesexecutable orca
     pesoptexecutable orca

     # Graph-driven sampling (GDS) control - by default, everything inactive!

     movefile moves.in
     gdsrestspring 0.1
     nbstrength 0.02
     nbrange 2.0
     kradius 0.1
     ngdsiter 500
     ngdsrelax 2000
     gdsdtrelax 0.05
     gdsoutfreq 10

     valencerange{
     C 1 4
     O 1 2
     H 0 1
     }

     reactiveatomtypes{
     C
     O
     H
     }

     reactiveatoms{
     range 1 4
     }

     reactivevalence{
     Fe Fe 6 8
     }

     fixedbonds{
     C O
     3 4
     }

     essentialatoms{
     range 1 5
     }

     essentialmoveatoms{
     id 1
     }

     reactiveatomtypes{
     C
     H
     }

     forbidgraphs .TRUE.
     forbidfile forbid.in

     # Pathfinder calculation
     nrxn 15
     nmcrxn 2000000
     mcrxntemp 1000000.0
     graphfunctype 4
     nmechmove 1
     minmolcharge 0
     maxmolcharge 0
     nchargemol 0
     maxstepcharge 0
     maxtotalcharge 0
