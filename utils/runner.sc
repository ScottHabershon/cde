#!/bin/bash -f
#############################
# After using the 'maker.sc' script, one can then submit all the jobs to the
# COW (in this case) using the script below.
# We assume that there is a jobscript 'script.sc' which controls job submission.
#############################

nfiles=100
for i in $(seq 1 $nfiles) ; do
  cd run_$i
  sbatch script.sc
  cd ../.
done

