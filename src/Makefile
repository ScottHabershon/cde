#
# Makefile for CDE
#

#
# COW compilation
#
# FC = ifort
# LIBS = -L/warwick/intel/mkl/10.0.5.025/lib/em64t/ -lmkl_intel_lp64 -lmkl_lapack -lmkl_sequential -lmkl_core -lpthread -lguide -static /warwick/mathlib/intel/x86_64/lib/libfftw3.a

#
# MAC compilation
#
mac: FC = gfortran-mp-9
mac: LIBS =  -framework Accelerate
mac: FFLAGS= -funroll-loops -Ofast   
mac: LFLAGS=	$(FFLAGS)

######

EXE = ../bin/cde.x
#FC = gfortran
#FFLAGS= -llapack
#FFLAGS= -O0 -g -fbacktrace -ffpe-trap=invalid,zero,overflow  -fbounds-check -llapack
FFLAGS= -g -O3 -funroll-loops -llapack
LFLAGS=	$(FFLAGS)

#EXE = cde.x

SRC= \
	constants.f90 \
	globaldata.f90 \
	functions.f90 \
	io.f90 \
	structure.f90 \
	pes.f90 \
	rpath.f90 \
	pathopt.f90 \
	pathfinder.f90 \
	main.f90


#############################################################
# CHANGE THINGS BELOW HERE AT YOUR PERIL
#############################################################

MF=	Makefile

.SUFFIXES:
.SUFFIXES: .f90 .o

OBJ=	$(SRC:.f90=.o)

.f90.o:
	$(FC) $(FFLAGS) -c $<

all:  $(EXE)
mac:	$(EXE)

$(EXE):	$(OBJ)
	$(FC) $(LFLAGS) -o $@ $(OBJ) $(LIBS)

$(OBJ):	$(MF)

tar:
	tar cvf $(EXE).tar $(MF) $(SRC)

clean:
	rm -f $(OBJ) $(EXE) core
