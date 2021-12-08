# Makefile for CDE

.SUFFIXES:
EXE = ./bin/cde.x
FC = ifort
LIBS =  -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_ilp64.a ${MKLROOT}/lib/intel64/libmkl_intel_thread.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -liomp5 -lpthread -lm -ldl
FFLAGS = -O3 -i8 -I"${MKLROOT}/include"

ifeq "$(FC)" "gfortran"
MODCMD = -J
else ifeq "$(FC)" "ifort"
MODCMD = -module
endif

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

prog: $(EXE) clean

$(EXE):	$(OBJ)
	mkdir -p bin
	$(FC) $(FFLAGS) -o $@ $(OBJ) $(LIBS)

%.o: %.f90
	$(FC) $(FFLAGS) -c $< -o $@ $(MODCMD) src

tar:
	tar cvf $(EXE).tar $(MF) $(SRC)

.PHONY: clean
clean:
	rm src/*.o src/*.mod
