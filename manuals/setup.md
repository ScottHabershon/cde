## Setting up CDE ## {#setup}

## Organization of CDE

The main CDE directory (referred to hereafter as ~/cde) has a series of subdirectories containing the following:
- ./src: contains the main source files
- ./docs: contains the documentation (generated using doxygen)
- ./bin: a directory where the cde.x executable is stored
- ./examples: a directory containing the CDE tutorials
- ./utils: a directory containining a few (possibly) useful plotting and analysis scripts, mostly written in
  python and using matplotlib.

## Compiling CDE

CDE is generally written in pretty standard Fortran90 throughout - it should be relatively easy to compile, and does not require any fancy libraries beyond LAPACK.

To compile CDE, just do the following:

- Go to the *src* directory.
- Edit the *Makefile* so that your Fortran compiler and LAPACK libraries are picked up. You can also edit the *Makefile* (if you want) to change the executable name and the location of the compiled executable (the default is *~/cde/bin/cde.x*).
- Type *make clean*
- Type *make*
- After compilation, you should now have an executable labelled *cde.x* in the specified location (probably *~/cde/bin/).
- At this point, it's a good idea (although not essential) to make sure that either *~/bin/* is in your environment's *PATH* variable, or you add an alias to run cde.x. From now on, the rest of the documentation will assume that typing *cde.x* will run the compiled CDE executable.

## Running CDE

To run the code, you'll need some input files. First, you should go and read the sections about the input files that CDE requires. 

Once your input files are ready, you simply type:

        cde.x input

where *input* is the name of your input file.

See the sections @ref Annotated and @ref Example for descriptions of input files.

## Adding an alias to your environment

To make things easier, you can add an alias to your environment's configuration file in your home directory. To do so (assuming a linux/unix system using BASH shell):

- Type *cd* to change to your home directory;
- Type *vi .bashrc*
- Add the line


	alias cde='YYY/cde.x'


where *YYY* is the full directory path to the *cde.x* executable.

- Save the *.bashrc* file, then type


	source .bashrc


- If successful, you should then be able to type *cde.x* to run the CDE code.

