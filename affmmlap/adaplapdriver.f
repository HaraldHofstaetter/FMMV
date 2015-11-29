cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
c
c  this is the driver of the adaptive fast multipole code.
c
c  the following variables should be defined in the driver.
c
c    natoms: total number of atoms.
c    zat(3,natoms): particle locations.
c    charge(natoms) : charge each particle carries.
c    pot(natoms): potential at each particle location.
c    field(3,natoms): field at each particle location.
c    nlev: number of levels in the octree structure.
c    ier: error message.
c
c  a subroutine is provided to compare fmm results with direct results,
c  so you also provide the following vectors for such comparisons.
c
c    dpot(natoms): potential using direct summation.
c    dfield(3,natoms): field using direct summation.
c
c  for advanced users:
c    if you want to change the precision of this code, you can check
c    parm-alap.h.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
c
      IMPLICIT NONE
c
      INTEGER *4 NATOMS,IER
      PARAMETER (NATOMS=50000)
      REAL *8 ZAT(3,NATOMS),CHARGE(NATOMS)
      REAL *8 POT(NATOMS),DPOT(NATOMS)
      REAL *8 FIELD(3,NATOMS),DFIELD(3,NATOMS)
c
c-----functions called.
c
      REAL *8 SECOND
c
c-----more local variables.
c
      INTEGER *4 I,IMAX
      REAL *8 CTOT,TIME0,TIME1,SALG,STOT,SALG2,STOT2
      REAL *8 TT,EALG,ERRMAX
c
c-----call initialization routines
c     1. initialize the output units.
c        6: screen.
c        13: fort.13.
c
      CALL PRINI(6,13)
c
c-----output parameters.
c
      CALL PRINF('natoms = *',NATOMS,1)
c
c-----2. generate charges and their locations.
c
      CALL DUMMY(NATOMS,ZAT,CHARGE)
c
c-----3. main fmm call, call fmm to calculate the potential
c
      TIME0 = SECOND()
      CALL FMMLAP_A(NATOMS,ZAT,CHARGE,POT,FIELD,IER)
      TIME1=SECOND()
c
      CALL PRIN2(' time for adaptive fmm is *',TIME1-TIME0,1)
      WRITE(11,554)TIME1-TIME0
554   FORMAT(' time for expansion work is ',F8.2)
c
c-----end of the fast multipole algorithm.
c
c-----4. finally, comparison with direct method.
c       write out first imax data points
c
      IMAX = NATOMS
      IF (NATOMS .GT. 400) IMAX = 400
c
c-----write out first imax data points
c
c      call prin2( 'from adapfmm, pot = *',pot,imax)
c      call prin2( 'from adapfmm, field = *',field,3*imax)
c
c------direct calculation.
c
      TIME0 = SECOND()
      DO 1100 I=1,IMAX
        CALL DIRECI(NATOMS,ZAT,CHARGE,I,DFIELD(1,I),DPOT(I))
        SALG = SALG + (DPOT(I) - POT(I))**2
        STOT = STOT + DPOT(I)*DPOT(I)
        IF (DABS(DPOT(I)-POT(I)).GT.ERRMAX) THEN
          ERRMAX=DABS(DPOT(I)-POT(I))
        ENDIF
1100  CONTINUE
      TIME1 = SECOND()
c
c-----error analysis
c
      SALG = 0.0D0
      SALG2= 0.0D0
      STOT = 0.0D0
      STOT2= 0.0D0
      ERRMAX=0.0D0
      DO 1200 I=1,IMAX
        SALG = SALG + (DPOT(I) - POT(I))**2
        STOT = STOT + DPOT(I)*DPOT(I)
        SALG2= SALG2+ (DFIELD(1,I)-FIELD(1,I))**2+
     1    (DFIELD(2,I)-FIELD(2,I))**2+(DFIELD(3,I)-FIELD(3,I))**2
        STOT2= STOT2+ DFIELD(1,I)*DFIELD(1,I)+
     1    DFIELD(2,I)*DFIELD(2,I)+DFIELD(3,I)*DFIELD(3,I)
        IF (DABS(DPOT(I)-POT(I)).GT.ERRMAX) THEN
          ERRMAX=DABS(DPOT(I)-POT(I))
        ENDIF
1200  CONTINUE
c
c-----output results.
c
      CALL PRIN2(' time for imax points directly is *',TIME1-TIME0,1)
      TT = (TIME1-TIME0)*NATOMS/IMAX
      CALL PRIN2(' time for all points directly is *',TT,1)
      WRITE(11,556)TT
556   FORMAT(' time for direct calc is ',F8.2)
c
c      call prin2( 'direct pot = *',dpot,imax)
c      call prin2( 'direct field = *',dfield,3*imax)
c
      EALG = DSQRT(SALG/STOT)
      CALL PRIN2( ' l2 error of potential = *',EALG,1)
      CALL PRIN2( ' max error of potential = *',ERRMAX,1)
      EALG = DSQRT(SALG2/STOT2)
      CALL PRIN2( ' l2 error of field = *',EALG,1)
c
      STOP
      END
