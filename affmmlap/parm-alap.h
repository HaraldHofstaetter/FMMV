cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
c
c  header file parm-alap.h
c
c    adjustable parameters for adaptive fast multipole method.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
c
c-----the boundary condition, currently only 0 allowed
c
      INTEGER *4 IFLAG
      PARAMETER (IFLAG=0)
c
c-----a large working array for generating the tree structure, the following
c       parameter gives the length of the structure.
c
      INTEGER *4 LW
      PARAMETER (LW=300000000)
c
c-----max allowed number of particles per box.
c
      INTEGER *4 NBOX
      PARAMETER (NBOX=80)
c
c-----fast multiple expansion parameters.
c
      INTEGER *4 NTERMS,NLAMBS
c
c-----for 3 digits accuracy.
c
      parameter (nterms=9,nlambs=9)
c
c-----for 6 digits accuracy.
c
c     PARAMETER (NTERMS=18,NLAMBS=18)
c
c-----the conversion ratio from integer to real.
c
      INTEGER *4 NINIRE
      PARAMETER (NINIRE=2)
c
c-----check if two particles are too close to each other
c
      REAL *8 EPSCLOSE
      PARAMETER (EPSCLOSE=1D-12)
c
c----------------------------------------------------------------------d
c
c     the following parameters should not be changed.
c     fixed memory requirements.
c
c----------------------------------------------------------------------d
c
c-----tree structure
c
      INTEGER *4 LADDR(2,200)
c
c-----local memory for charge, can be eliminated later.
c
      REAL *8 ZAT2(3,NBOX),CHARG2(NBOX),DN2(3,NBOX)
c
c-----precomputed fmm stuff
c
      INTEGER *4 NUMPHYS(NLAMBS),NUMFOUR(NLAMBS)
c
      REAL *8 WHTS(NLAMBS),RLAMS(NLAMBS)
      REAL *8 CARRAY((4*NTERMS+1)*(4*NTERMS+1))
      REAL *8 DC((2*NTERMS+1)*(2*NTERMS+1)*(2*NTERMS+1))
      REAL *8 RDPLUS((NTERMS+1)*(NTERMS+1)*(2*NTERMS+1))
      REAL *8 RDMINUS((NTERMS+1)*(NTERMS+1)*(2*NTERMS+1))
      REAL *8 RDSQ3((NTERMS+1)*(NTERMS+1)*(2*NTERMS+1))
      REAL *8 RDMSQ3((NTERMS+1)*(NTERMS+1)*(2*NTERMS+1))
      REAL *8 C(0:60,0:60),CS(0:60,0:60),CSINV(0:60,0:60)
      REAL *8 RLSC(NLAMBS*(1+NTERMS)*(1+NTERMS))
      REAL *8 YTOP(0:60,0:60)
c
c-----expansion stuff.
c
      COMPLEX *16 MW1((NTERMS+1)*(NTERMS+1))
      COMPLEX *16 MW2((NTERMS+1)*(NTERMS+1))
      COMPLEX *16 MW3((NTERMS+1)*(NTERMS+1))
