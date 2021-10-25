#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=1096mb
#SBATCH --time=06:00:00

/storage/chem/msslap/Permute/src/cde.x input

