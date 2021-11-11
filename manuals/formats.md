### I/O structure formats {#formats}

Some important things to note about structure formats used in CDE:

- Throughout CDE, we generally use simple *.xyz files as input/output formats.

- You can use *openbabel* (http://openbabel.org/wiki/Main_Page) to convert from a wide range of files to/from xyz files.

- A convenient program for creating initial structures in xyz format is *Avogadro* (https://avogadro.cc).

- The coordinates in an XYZ format file in units of Angstroms.

- If performing NEB calculations, the number of snapshots (or configurations) in the xyz file for the full path should be the same as that defined in *nimage* in the input file.

- Each xyz frame or file looks like this:

		2

		C 0.000  0.000 0.0000
		O 1.50   1.670 5.6708


- The first line contains the number of atoms, the second line is blank (or a comment line) and then there is a list of atoms with their atomic labels and (x,y,z) coordinates in Angstroms.

- Different frames in an xyz file follow each other immediately, like this:

		2

		C 0.000  0.000 0.0000
		O 1.50   1.670 5.6708
		2

		C 0.1500  0.2166 0.007
		O 1.56  1.679  5.7098
