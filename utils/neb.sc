#!/bin/bash -f
###################################################################
# This example script creates directories for all nimage=10 NEB runs following a previous
# double-ended GDS calculation. The initial NEB paths for each nimage reaction are in
# final_path_rx_$i.xyz, and directories path_$i are created for each run.
# script_neb contains the job submission script.
###################################################################
for i in {1..10}
do
mkdir path_$i
cp input_neb ./path_$i/.
cp dftb.min ./path_$i/.
cp dftb.head ./path_$i/.
cp moves_new.in ./path_$i/.
cp final_path_rx_$i.xyz ./path_$i/path.xyz
cp start.xyz ./path_$i/.
cp end_new.xyz ./path_$i/.
cp script_neb ./path_$i/.
cd path_$i
sbatch script_neb
cd ../.
done


