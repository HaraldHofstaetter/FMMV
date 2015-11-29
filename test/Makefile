include ../Makefile.inc

INCDIR = ../include
LIBDIR = ../lib

all: test_fmmv3d test_fmmv3df \
     test_fmmv2d test_fmmv2df \
     test_M_L_3d test_M_L_3df \
     test_M_L_2d test_M_L_2df \
     test_X_3d test_X_3df \
     test_X_2d test_X_2df \


test_utilities_f.o: test_utilities.c
	$(CC) $(CFLAGS) -I$(INCDIR)  -o test_utilities_f.o -c test_utilities.c -DUSE_SINGLE_PRECISION 


test_utilities.o: test_utilities.c
	$(CC) $(CFLAGS) -I$(INCDIR)  -o test_utilities.o -c test_utilities.c 

test_fmmv3df.o: test_fmmv.c
	$(CC) $(CFLAGS) -I$(INCDIR) -DFMM_DIM=3 -DUSE_SINGLE_PRECISION \
        -o test_fmmv3df.o -c test_fmmv.c 

test_fmmv3d.o: test_fmmv.c
	$(CC) $(CFLAGS) -I$(INCDIR) -DFMM_DIM=3  \
        -o test_fmmv3d.o -c test_fmmv.c 


test_fmmv3d: test_fmmv3d.o test_utilities.o
	$(LN) $(TEST_LN_FLAGS) -o test_fmmv3d \
	test_fmmv3d.o test_utilities.o \
        -L$(LIBDIR) -lfmmv3d  $(LIBS)

test_fmmv3df: test_fmmv3df.o test_utilities_f.o
	$(LN) $(TEST_LN_FLAGS) -o test_fmmv3df \
	test_fmmv3df.o test_utilities_f.o \
        -L$(LIBDIR) -lfmmv3df  $(LIBS)


test_fmmv2df.o: test_fmmv.c
	$(CC) $(CFLAGS) -I$(INCDIR) -DFMM_DIM=2 -DUSE_SINGLE_PRECISION \
        -o test_fmmv2df.o -c test_fmmv.c 

test_fmmv2d.o: test_fmmv.c
	$(CC) $(CFLAGS) -I$(INCDIR) -DFMM_DIM=2  \
        -o test_fmmv2d.o -c test_fmmv.c 


test_fmmv2d: test_fmmv2d.o test_utilities.o
	$(LN) $(TEST_LN_FLAGS) -o test_fmmv2d \
	test_fmmv2d.o test_utilities.o \
        -L$(LIBDIR) -lfmmv2d  $(LIBS)

test_fmmv2df: test_fmmv2df.o test_utilities_f.o
	$(LN) $(TEST_LN_FLAGS) -o test_fmmv2df \
	test_fmmv2df.o test_utilities_f.o \
        -L$(LIBDIR) -lfmmv2df  $(LIBS)







test_M_L_3d.o: test_M_L_3d.c
	$(CC) $(CFLAGS) -I$(INCDIR) -I.. -I../3d -o test_M_L_3d.o -c test_M_L_3d.c 

test_M_L_3df.o: test_M_L_3d.c
	$(CC) $(CFLAGS) -I$(INCDIR) -I.. -I../3d -o test_M_L_3df.o -c test_M_L_3d.c -DUSE_SINGLE_PRECISION

test_M_L_3d: test_M_L_3d.o test_utilities.o
	$(LN) $(TEST_LN_FLAGS) -o test_M_L_3d \
	test_M_L_3d.o test_utilities.o \
        -L$(LIBDIR) -lfmmv3d  $(LIBS)

test_M_L_3df: test_M_L_3df.o test_utilities_f.o
	$(LN) $(TEST_LN_FLAGS) -o test_M_L_3df \
	test_M_L_3df.o test_utilities_f.o \
        -L$(LIBDIR) -lfmmv3df  $(LIBS)

test_M_L_2d.o: test_M_L_2d.c
	$(CC) $(CFLAGS) -I$(INCDIR) -I.. -I../2d -o test_M_L_2d.o -c test_M_L_2d.c 

test_M_L_2df.o: test_M_L_2d.c
	$(CC) $(CFLAGS) -I$(INCDIR) -I.. -I../2d -o test_M_L_2df.o -c test_M_L_2d.c -DUSE_SINGLE_PRECISION

test_M_L_2d: test_M_L_2d.o test_utilities.o
	$(LN) $(TEST_LN_FLAGS) -o test_M_L_2d \
	test_M_L_2d.o test_utilities.o \
        -L$(LIBDIR) -lfmmv2d  $(LIBS)

test_M_L_2df: test_M_L_2df.o test_utilities_f.o
	$(LN) $(TEST_LN_FLAGS) -o test_M_L_2df \
	test_M_L_2df.o test_utilities_f.o \
        -L$(LIBDIR) -lfmmv2df  $(LIBS)



test_X_3d.o: test_X_3d.c
	$(CC) $(CFLAGS) -I$(INCDIR) -I.. -I../3d -o test_X_3d.o -c test_X_3d.c 

test_X_3df.o: test_X_3d.c
	$(CC) $(CFLAGS) -I$(INCDIR) -I.. -I../3d -o test_X_3df.o -c test_X_3d.c -DUSE_SINGLE_PRECISION

test_X_3d: test_X_3d.o test_utilities.o
	$(LN) $(TEST_LN_FLAGS) -o test_X_3d \
	test_X_3d.o test_utilities.o \
        -L$(LIBDIR) -lfmmv3d  $(LIBS)

test_X_3df: test_X_3df.o test_utilities_f.o
	$(LN) $(TEST_LN_FLAGS) -o test_X_3df \
	test_X_3df.o test_utilities_f.o \
        -L$(LIBDIR) -lfmmv3df  $(LIBS)


test_X_2d.o: test_X_2d.c
	$(CC) $(CFLAGS) -I$(INCDIR) -I.. -I../2d -o test_X_2d.o -c test_X_2d.c 

test_X_2df.o: test_X_2d.c
	$(CC) $(CFLAGS) -I$(INCDIR) -I.. -I../2d -o test_X_2df.o -c test_X_2d.c -DUSE_SINGLE_PRECISION

test_X_2d: test_X_2d.o test_utilities.o
	$(LN) $(TEST_LN_FLAGS) -o test_X_2d \
	test_X_2d.o test_utilities.o \
        -L$(LIBDIR) -lfmmv2d  $(LIBS)

test_X_2df: test_X_2df.o test_utilities_f.o
	$(LN) $(TEST_LN_FLAGS) -o test_X_2df \
	test_X_2df.o test_utilities_f.o \
        -L$(LIBDIR) -lfmmv2df  $(LIBS)



clean:
	rm -f *.o
	rm -f *.mod
	rm -f test_fmmv3d test_fmmv3df
	rm -f test_fmmv2d test_fmmv2df
	rm -f test_M_L_3d test_M_L_3df 
	rm -f test_M_L_2d test_M_L_2df 
	rm -f test_X_3d test_X_3df
	rm -f test_X_2d test_X_2df


