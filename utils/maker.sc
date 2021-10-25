#!/bin/bash -f
######################################
# Example script to create 100 independent runs with different random number seeds, and then
# copy all the relevant files for a CDE run into each directory (assuming they all exist
# in a root directory from which this script is run.
#
# nfiles controls the number of directories to create.
#####################################

nfiles=100
for i in $(seq 1 $nfiles) ; do
  mkdir run_$i
  echo "ranseed $i" > out
  cat input.root out > ./run_$i/input
  cp lmp.min ./run_$i/.
  cp lmp.head ./run_$i/.
  cp lmp_control ./run_$i/.
  cp CHOPtNiX.ff ./run_$i/.
  cp moves.in ./run_$i/.
  cp start.xyz ./run_$i/.
  cp end.xyz ./run_$i/.
  cp script.sc ./run_$i/.
done

