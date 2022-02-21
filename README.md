**CDE: Chemistry discovery engine**

CDE is a set of fortran routines which implements several
different calculation types associated with chemical reaction-path analysis.

CDE can perform several types of calculations, including:
- Single-ended Graph-driven sampling for chemical reactions,
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

The graph-driven sampling of reaction-paths implemented in CDE is initially described in the following articles:
 - Automated prediction of catalytic mechanism and rate law using graph-based reaction path sampling, S. Habershon, Journal of chemical theory and computation 12 (4), 1786-1798
 - Sampling reactive pathways with random walks in chemical space: Applications to molecular dissociation and catalysis, S Habershon, Journal of chemical physics 143 (9), 094106
