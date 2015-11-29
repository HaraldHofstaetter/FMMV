CC = gcc
CFLAGS =  -g -fPIC 

MATLABPATH = /usr/local/MATLAB/R2013a/
MEX = $(MATLABPATH)/bin/mex  -O -largeArrayDims
MATLABINC = $(MATLABPATH)/extern/include

#MEXSFX = mexglx
MEXSFX = mexa64 
FMMVINC = ../include
FMMV2DLIB = -L ../lib -lfmmv2d
FMMV3DLIB = -L ../lib -lfmmv3d

all: fmmv2d.$(MEXSFX) fmmv2d_initialize.$(MEXSFX) fmmv2d_evaluate.$(MEXSFX) fmmv2d_finalize.$(MEXSFX) \
     fmmv3d.$(MEXSFX) fmmv3d_initialize.$(MEXSFX) fmmv3d_evaluate.$(MEXSFX) fmmv3d_finalize.$(MEXSFX)

fmmv2d.o : fmmv.c
	$(CC) $(CFLAGS) -DDIM=2 -I$(FMMVINC) -I$(MATLABINC) -c -o fmmv2d.o fmmv.c

fmmv2d_initialize.o : fmmv_initialize.c
	$(CC) $(CFLAGS) -DDIM=2 -I$(FMMVINC) -I$(MATLABINC) -c -o fmmv2d_initialize.o fmmv_initialize.c

fmmv2d_evaluate.o : fmmv_evaluate.c
	$(CC) $(CFLAGS) -DDIM=2 -I$(FMMVINC) -I$(MATLABINC) -c -o fmmv2d_evaluate.o fmmv_evaluate.c

fmmv2d_finalize.o : fmmv_finalize.c
	$(CC) $(CFLAGS) -DDIM=2 -I$(FMMVINC) -I$(MATLABINC) -c -o fmmv2d_finalize.o fmmv_finalize.c

statistics2d.o : statistics.c
	$(CC) $(CFLAGS) -DDIM=2 -I$(FMMVINC) -I$(MATLABINC) -c -o statistics2d.o statistics.c
	
fmmv2d.$(MEXSFX): fmmv2d.o statistics2d.o
	$(MEX) -output fmmv2d fmmv2d.o statistics2d.o $(FMMV2DLIB)

fmmv2d_initialize.$(MEXSFX): fmmv2d_initialize.o statistics2d.o
	$(MEX) -output fmmv2d_initialize fmmv2d_initialize.o statistics2d.o $(FMMV2DLIB)

fmmv2d_evaluate.$(MEXSFX): fmmv2d_evaluate.o statistics2d.o
	$(MEX) -output fmmv2d_evaluate fmmv2d_evaluate.o statistics2d.o $(FMMV2DLIB)

fmmv2d_finalize.$(MEXSFX): fmmv2d_finalize.o statistics2d.o
	$(MEX) -output fmmv2d_finalize fmmv2d_finalize.o statistics2d.o $(FMMV2DLIB)



fmmv3d.o : fmmv.c
	$(CC) $(CFLAGS) -DDIM=3 -I$(FMMVINC) -I$(MATLABINC) -c -o fmmv3d.o fmmv.c

fmmv3d_initialize.o : fmmv_initialize.c
	$(CC) $(CFLAGS) -DDIM=3 -I$(FMMVINC) -I$(MATLABINC) -c -o fmmv3d_initialize.o fmmv_initialize.c

fmmv3d_evaluate.o : fmmv_evaluate.c
	$(CC) $(CFLAGS) -DDIM=3 -I$(FMMVINC) -I$(MATLABINC) -c -o fmmv3d_evaluate.o fmmv_evaluate.c

fmmv3d_finalize.o : fmmv_finalize.c
	$(CC) $(CFLAGS) -DDIM=3 -I$(FMMVINC) -I$(MATLABINC) -c -o fmmv3d_finalize.o fmmv_finalize.c

statistics3d.o : statistics.c
	$(CC) $(CFLAGS) -DDIM=3 -I$(FMMVINC) -I$(MATLABINC) -c -o statistics3d.o statistics.c
	
fmmv3d.$(MEXSFX): fmmv3d.o statistics3d.o
	$(MEX) -output fmmv3d fmmv3d.o statistics3d.o $(FMMV3DLIB)

fmmv3d_initialize.$(MEXSFX): fmmv3d_initialize.o statistics3d.o
	$(MEX) -output fmmv3d_initialize fmmv3d_initialize.o statistics3d.o $(FMMV3DLIB)

fmmv3d_evaluate.$(MEXSFX): fmmv3d_evaluate.o statistics3d.o
	$(MEX) -output fmmv3d_evaluate fmmv3d_evaluate.o statistics3d.o $(FMMV3DLIB)

fmmv3d_finalize.$(MEXSFX): fmmv3d_finalize.o statistics3d.o
	$(MEX) -output fmmv3d_finalize fmmv3d_finalize.o statistics3d.o $(FMMV3DLIB)

clean:
	rm -f *.o *.$(MEXSFX)

