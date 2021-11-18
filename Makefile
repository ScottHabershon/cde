# Makefile for CDE

.SUFFIXES:
EXE = cde.x
FC = ifort
LIBS = -lmkl_intel_lp64 -lmkl_lapack -lmkl_sequential -lmkl_core -lpthread -static
FFLAGS = -O3 -fopenmp

MF = Makefile
SRC = \
	src/constants.f90 \
	src/globaldata.f90 \
	src/functions.f90 \
	src/io.f90 \
	src/structure.f90 \
	src/pes.f90 \
	src/rpath.f90 \
	src/pathopt.f90 \
	src/pathfinder.f90 \
	src/main.f90

OBJ = $(SRC:%.f90=%.o)

$(EXE):	$(OBJ)
	$(FC) $(FFLAGS) -o $@ $(OBJ) $(LIBS)

%.o: %.f90
	$(FC) $(FFLAGS) -c $< -o $@ -module src

tar:
	tar cvf $(EXE).tar $(MF) $(SRC)

.PHONY: clean
clean:
	rm src/*.o src/*.mod
