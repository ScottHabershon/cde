###Useful scripts {#scripts}

The directory *~/utils/* in the main CDE installation directory contains a series of useful *bash* and *python* scripts as follows:

- **CompareProfiles.py**: Uses python and matplotlib to compare two energy profiles resulting from NEB calculations.
- **WholeEnergyProfile.py**: Plots a sequence of NEB results as a single continuous energy profile using python and matplotlib.
- **maker.sc**: A bash script for creating multiple directories with identicial copies of a set of input files, except with different random number seeds.
- **runner.sc**: A bash script to submit multiple jobs to the HPC queue. This can be used after using *maker.sc* to generate the input directories and files.
- **neb.sc**: A bash script to perform multiple NEB runs. This can be used after a double-ended GDS calculation to calculate the energy along the minimum energy path for a whole sequence of reactions. The final result of using *neb.sc* can then be plotted using *WholeEnergyProfile.py*.
- **job_script.sc**: An example queue script for submission of a job to the Warwick COW.
