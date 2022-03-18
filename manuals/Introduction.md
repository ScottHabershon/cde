\mainpage CDE: Chemical discovery engine

CDE is a set of fortran routines which implement several 
different calculation types associated with chemical reaction-path analysis.

*This is a legacy version of our reaction discovery code - an updated python version is in progress.*

CDE can perform several types of calculations, including:
- Single-ended Graph-driven sampling for chemical reactions (CURRENTLY UNAVAILABLE),
- Double-ended mechanism searching,
- Generation of initial approximate MEPs connecting user-defined reactant and
  product structures using a variety of methods, including linear interpolation and 
  image-dependent pair potential (IDPP) method.
- Nudged elastic band calculations for refinement of chemical reaction paths.

CDE interfaces with several external programs to perform energy evaluations and geometry optimization, including:
- ORCA
- Psi-4
- LAMMPS
- DFTB+
- Molpro

See the **references** section for further references to methods employed.
