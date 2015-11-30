include ../Makefile.inc

all: test_fmmv3d test_fmmv3df  \
     test_fmmv2d test_fmmv2df 


test_utilities.o test_utilities.mod: test_utilities.F90
	$(FC) $(FFLAGS) -o test_utilities.o -c test_utilities.F90

test_utilities_f.o test_utilities_f.mod: test_utilities.F90
	$(FC) $(FFLAGS) -DUSE_SINGLE_PRECISION \
	-o test_utilities_f.o -c test_utilities.F90


fmmv3d_module.mod: fmmv3d.f90
	$(FC) $(FFLAGS) -c fmmv3d.f90

test_fmmv3d.o: test_fmmv.F90 fmmv3d_module.mod test_utilities.mod
	$(FC) $(FFLAGS) -DFMM_DIM=3 \
        -o test_fmmv3d.o -c test_fmmv.F90

test_fmmv3df.o: test_fmmv.F90 fmmv3d_module.mod test_utilities_f.mod
	$(FC) $(FFLAGS) -DFMM_DIM=3 -DUSE_SINGLE_PRECISION \
	-o test_fmmv3df.o -c test_fmmv.F90

test_fmmv3d: test_fmmv3d.o test_utilities.o
	$(FC) $(TEST_LN_FLAGS) -o test_fmmv3d test_fmmv3d.o \
	test_utilities.o \
        -L$(LIB_INSTALLDIR) -lfmmv3d  $(LIBS)

test_fmmv3df: fmmv3d.o test_fmmv3df.o test_utilities_f.o
	$(FC) $(TEST_LN_FLAGS) -o test_fmmv3df test_fmmv3df.o \
	fmmv3d.o test_utilities_f.o \
        -L$(LIB_INSTALLDIR) -lfmmv3df  $(LIBS)



fmmv2d_module.mod: fmmv2d.f90
	$(FC) $(FFLAGS) -c fmmv2d.f90

test_fmmv2d.o: test_fmmv.F90 fmmv2d_module.mod test_utilities.mod
	$(FC) $(FFLAGS) -DFMM_DIM=2 \
        -o test_fmmv2d.o -c test_fmmv.F90

test_fmmv2df.o: test_fmmv.F90 fmmv2d_module.mod test_utilities_f.mod
	$(FC) $(FFLAGS) -DFMM_DIM=2 -DUSE_SINGLE_PRECISION \
	-o test_fmmv2df.o -c test_fmmv.F90

test_fmmv2d:  test_fmmv2d.o test_utilities.o
	$(FC) $(TEST_LN_FLAGS) -o test_fmmv2d test_fmmv2d.o \
	test_utilities.o \
        -L$(LIB_INSTALLDIR) -lfmmv2d  $(LIBS)

test_fmmv2df: fmmv2d.o test_fmmv2df.o test_utilities_f.o
	$(FC) $(TEST_LN_FLAGS) -o test_fmmv2df test_fmmv2df.o \
	fmmv2d.o test_utilities_f.o \
        -L$(LIB_INSTALLDIR) -lfmmv2df  $(LIBS)


clean:
	rm -f *.o
	rm -f *.mod
	rm -f test_fmmv3d test_fmmv3df
	rm -f test_fmmv2d test_fmmv2df

