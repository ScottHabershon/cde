#/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=2012mb
#SBATCH --time=24:00:00

/storage/chem/msslap/ISM/src/cde.x input > output 
