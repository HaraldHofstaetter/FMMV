include ../Makefile.inc

LIBDIR = ../lib

all: test_huang

HUANG_OBJS    = \
	fmmadaplap.o \
	treeadap.o \
	laptable.o \
	lapoperators.o \
	lapadap.o \
	prini.o \
	second.o \
	functions.o  

# 5. SUFFIX RULE
.f.o:
	$(FC) $(FFLAGS) -c $*.f

fmmv3d.o fmmv3d_module.mod: ../fortran/fmmv3d.f90
	$(FC) $(FFLAGS) -o fmmv3d.o -c ../fortran/fmmv3d.f90

test_utilities.o test_utilities.mod: ../fortran/test_utilities.F90
	$(FC) $(FFLAGS) -o test_utilities.o -c ../fortran/test_utilities.F90

test_huang.o: test_huang.f90 test_utilities.mod fmmv3d_module.mod
	$(FC) $(FFLAGS) -o test_huang.o -c test_huang.f90


test_huang: fmmv3d.o test_huang.o test_utilities.o $(HUANG_OBJS)
	$(FC) $(TEST_LN_FLAGS) -o test_huang test_huang.o \
	fmmv3d.o test_utilities.o $(HUANG_OBJS) \
        -L$(LIBDIR) -lfmmv3d  $(LIBS)

clean:
	rm -f *.o
	rm -f *.mod
	rm -f test_huang

