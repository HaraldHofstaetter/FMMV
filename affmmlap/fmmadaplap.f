cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
      SUBROUTINE FMMLAP_A(NATOMS,ZAT,CHARGE,POT,FIELD,IER)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
c
c  purpose:
c    this is the main calling subroutine. the purpose of this
c    subroutine is to provide an interface between the fmm main
c    subroutines and how users prefer to use the code.
c
c  on input:
c    natoms: total number of atoms.
c    zat(3,natoms): particle locations.
c    charge(natoms) : charge each particle carries.
c
c  on output:
c    pot(natoms): potential at each particle location.
c    field(3,natoms): the force field at each particle location.
c    ier: error message.
c
c  other parameters:
c    a few parameters are specified in the file parm-alap.h,
c    including:
c      iflag: for boundary condition, currently only free space
c             boundary condition is implemented. periodic will be
c             added later.
c      lw: the length of the work space for generating the adaptive tree
c          structure.
c      nbox: the max number of particles in the childless box.
c      nterms: the number of terms in the multipole and local expansion.
c      nlambs: the number of terms in the exponential expansion.
c
c    the default is set to nterms=nlambs=9 for 3 digits accuracy.
c
c  memory allocation.
c    the memory is divided to two parts, fixed size part and allocated
c    part.
c
c    for fixed size part, check parm-alap.h.
c    memory allocation is done by calculating the total memory required
c      for integer, real and complex, and three big vectors are
c      allocated.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
c
      IMPLICIT NONE
c
c-----include parm-alap.h for additional parameters and fixed lengh
c       stuff.
c
      INCLUDE "parm-alap.h"
c
      INTEGER *4 NATOMS,IER
      REAL *8 ZAT(3,NATOMS),CHARGE(NATOMS)
      REAL *8 POT(NATOMS),FIELD(3,NATOMS)
c
c-----local variables important for memory allocation.
c
      INTEGER *4 NEXPMAX,NEXPTOT,NEXPTOTP,NTHMAX,NBOXES
      INTEGER *4 IIEXP1,IIEXP2,IIZ,IIWORK,IZS
      INTEGER *4 IXS,IYS,IMEXPF1,IMEXPF2,IMEXPP1,IMEXPP2,IMEXPPALL,
     1           IFEXPE,IFEXPO,IFEXPBACK,IMPOLE,ILOCAL,ILEXP1,ILEXP2
      INTEGER *4 IPTR,LASTI,LASTR,LASTC
c
c-----fmm working array.
c
      INTEGER *4, ALLOCATABLE :: IWORK(:)
      INTEGER *4, ALLOCATABLE :: IXFMM(:)
      REAL *8, ALLOCATABLE :: XRFMM(:)
      COMPLEX *16, ALLOCATABLE :: XCFMM(:)
c
c-----parameters for the adaptive tree structure.
c
      INTEGER *4 NLEV,LUSED
      REAL *8 CENTER0(3),SIZE
c
c-----functions called.
c
      REAL *8 SECOND
c
c-----more local variables.
c
      INTEGER *4 I,JJ
      REAL *8 TIME0,TIME1,TOTMEM,RMEM
c
c-----0. calculate the memory used. it contains the
c        following parts.
c        (a)input/output. (b) fixed fmm memory. (c) allocated fmm memory.
c        the first two parts will be calculated first. and (c) will be
c        added later.
c
      TOTMEM=DBLE(64*NATOMS)/1024D0/1024D0
      CALL PRIN2('fmm input/output memory (mb) *',TOTMEM,1)
      RMEM=DBLE(4*(170+NLAMBS*2)+8*((4*NTERMS+1)**2+
     1  3*(2*NTERMS+1)*(2*NTERMS+1)*(2*NTERMS+1)+3600*3+
     2  NLAMBS*(1+NTERMS)*(1+NTERMS))+16*((NTERMS+1)**2)
     3  *3)/1024D0/1024D0
      TOTMEM=TOTMEM+RMEM
      CALL PRIN2('fmm fixed memory (mb) *',RMEM,1)
c
c-----1. allocate a huge memory for generating tree structures, based on input
c       source and targets.
c
      TIME0=SECOND()
      ALLOCATE(IWORK(LW),STAT=IER)
      RMEM=DBLE(LW*4)/1024D0/1024D0
      CALL PRIN2('tree memory allocated (mb) *',RMEM,1)
      TOTMEM=TOTMEM+RMEM
      CALL PRIN2('total memory (mb) *',TOTMEM,1)
c
      IF (IER .NE. 0) THEN
        PRINT *, 'ier is', IER
        STOP 'allocation error when generating adaptive tree'
      ENDIF
c
c-----2. generate adaptive tree structure.
c
      CALL D3MSTRCR(IER,ZAT,NATOMS,NBOX,NBOXES,IWORK(1),LADDR,NLEV,
     1  CENTER0,SIZE,IWORK(NATOMS+5),LW-NATOMS-5,LUSED,EPSCLOSE,NINIRE)
      TIME1=SECOND()
      CALL PRINF('adaptive tree generated, ier= *',IER,1)
      CALL PRINF('number of levels = *',NLEV, 1)
      CALL PRINF('number of boxes = *',NBOXES, 1)
      CALL PRIN2('time to generate tree =*',TIME1-TIME0,1)
      IF (IER .NE. 0) THEN
        CALL PRINF('fatal: not enough memory in adaptive tree',
     2    IER,1)
        STOP
      ENDIF
c
c-----3. generate precomputed matrices.
c        this part should be optimized using quadrature code.
c
      CALL FRMINI(C,CS,CSINV)
      CALL ROTGEN(NTERMS,CARRAY,RDPLUS,RDMINUS,RDSQ3,RDMSQ3,DC)
      CALL VWTS(RLAMS,WHTS,NLAMBS)
      CALL NUMTHETAHALF(NUMFOUR,NLAMBS)
      CALL NUMTHETAFOUR(NUMPHYS,NLAMBS)
      CALL RLSCINI(RLSC,NLAMBS,RLAMS,NTERMS)
c
      NEXPTOT = 0
      NTHMAX = 0
      NEXPTOTP = 0
      DO I = 1,NLAMBS
        NEXPTOT = NEXPTOT + NUMFOUR(I)
        IF (NUMFOUR(I).GT.NTHMAX) NTHMAX = NUMFOUR(I)
        NEXPTOTP = NEXPTOTP + NUMPHYS(I)
      ENDDO
      NEXPTOTP = NEXPTOTP/2
c
      NEXPMAX=MAX(NEXPTOT,NEXPTOTP)+1
c
c-----4. now make pointers to fmm memory blocks.
c
      IPTR=1
      CALL ADJ_MEM_PTR(IPTR,IIEXP1,NBOXES)
      CALL ADJ_MEM_PTR(IPTR,IIEXP2,NBOXES)
      CALL ADJ_MEM_PTR(IPTR,IIZ,NATOMS)
      CALL ADJ_MEM_PTR(IPTR,IIWORK,LUSED)
      LASTI=IPTR
c
      IPTR=1
      CALL ADJ_MEM_PTR(IPTR,IZS,3*NEXPMAX)
      LASTR=IPTR
c
      IPTR=1
      CALL ADJ_MEM_PTR(IPTR,IXS,3*NEXPMAX)
      CALL ADJ_MEM_PTR(IPTR,IYS,3*NEXPMAX)
      CALL ADJ_MEM_PTR(IPTR,IMEXPF1,NEXPMAX)
      CALL ADJ_MEM_PTR(IPTR,IMEXPF2,NEXPMAX)
      CALL ADJ_MEM_PTR(IPTR,IMEXPP1,NEXPMAX)
      CALL ADJ_MEM_PTR(IPTR,IMEXPP2,NEXPMAX)
      CALL ADJ_MEM_PTR(IPTR,IMEXPPALL,16*NEXPMAX)
      CALL ADJ_MEM_PTR(IPTR,IFEXPE,15000)
      CALL ADJ_MEM_PTR(IPTR,IFEXPO,15000)
      CALL ADJ_MEM_PTR(IPTR,IFEXPBACK,15000)
      CALL ADJ_MEM_PTR(IPTR,IMPOLE,(NTERMS+1)*(NTERMS+1)*NBOXES)
      CALL ADJ_MEM_PTR(IPTR,ILOCAL,(NTERMS+1)*(NTERMS+1)*NBOXES)
      CALL ADJ_MEM_PTR(IPTR,ILEXP1,NEXPMAX*NBOXES)
      CALL ADJ_MEM_PTR(IPTR,ILEXP2,NEXPMAX*NBOXES)
      LASTC=IPTR
c
c-----5. allocate integer memory and deallocate the tree structure memory.
c
      ALLOCATE(IXFMM(LASTI),STAT=IER)
      IF (IER .NE. 0) THEN
        STOP 'integer array allocation error'
      ENDIF
c
      CALL ICOPY(NATOMS,IWORK(1),IXFMM(IIZ))
      CALL ICOPY(LUSED,IWORK(NATOMS+5),IXFMM(IIWORK))
c
      DEALLOCATE(IWORK)
      RMEM=DBLE(LW*4D0)/1024D0/1024D0
      TOTMEM=TOTMEM-RMEM
      CALL PRIN2('tree memory deallocated (mb) *',RMEM,1)
      CALL PRIN2('total memory (mb) *',TOTMEM,1)
c
      ALLOCATE(XRFMM(LASTR),XCFMM(LASTC),STAT=IER)
      RMEM=DBLE(LASTR*8D0+LASTC*16D0)/1024D0/1024D0
      TOTMEM=TOTMEM+RMEM
      CALL PRIN2('fmm real and complex memory (mb) allocated*',RMEM,1)
      CALL PRIN2('total memory (mb) *',TOTMEM,1)
c
      IF (IER .NE. 0) THEN
        CALL PRINF('allocation error in solvpb, ier=*',IER,1)
        STOP
      ENDIF
c
c-----6. more precomputing stuff.
c
      CALL MKEXPS(RLAMS,NLAMBS,NUMPHYS,NEXPTOTP,XCFMM(IXS),
     1  XCFMM(IYS),XRFMM(IZS))
      CALL MKFEXP(NLAMBS,NUMFOUR,NUMPHYS,XCFMM(IFEXPE),
     1  XCFMM(IFEXPO),XCFMM(IFEXPBACK))
c
c-----7. now the fast multipole calculation.
c
      TIME0 = SECOND()
      CALL LADAPFMM(IFLAG,ZAT,NATOMS,CHARGE,NTERMS,
     1  ZAT2,CHARG2,XCFMM(IMPOLE),XCFMM(ILOCAL),IXFMM(IIEXP1),
     2  XCFMM(ILEXP1),IXFMM(IIEXP2),XCFMM(ILEXP2),POT,FIELD,DC,RDPLUS,
     3  RDMINUS,RDMSQ3,RDSQ3,WHTS,RLAMS,NLAMBS,NUMFOUR,NUMPHYS,NEXPTOT,
     4  NEXPTOTP,NTHMAX,XCFMM(IMEXPF1),XCFMM(IMEXPF2),XCFMM(IMEXPP1),
     5  XCFMM(IMEXPP2),XCFMM(IMEXPPALL),XCFMM(IXS),XCFMM(IYS),
     6  XRFMM(IZS),MW1,MW2,MW3,XCFMM(IFEXPE),XCFMM(IFEXPO),
     7  XCFMM(IFEXPBACK),RLSC,C,CS,CSINV,NBOX,NBOXES,IXFMM(IIZ),LADDR,
     8  NLEV,IXFMM(IIWORK),LUSED,NINIRE,SIZE)
c
      TIME1 = SECOND()
      WRITE(11,554)TIME1-TIME0
554   FORMAT(' time for fmm work is ',F8.2)
c
c-----memory deallocation.
c
      deallocate(ixfmm,xrfmm,xcfmm)
c
      RETURN
      END
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
      SUBROUTINE LADAPFMM(IFLAG,ZAT,NATOMS,CHARGE,NTERMS,
     1  ZAT2,CHARG2,MPOLE,LOCAL,IEXP1,LEXP1,IEXP2,LEXP2,
     2  POT,FIELD,DC,RDPLUS,RDMINUS,RDMSQ3,RDSQ3,WHTS,RLAMS,
     3  NLAMBS,NUMFOUR,NUMPHYS,NEXPTOT,NEXPTOTP,NTHMAX,MEXPF1,
     4  MEXPF2,MEXPP1,MEXPP2,MEXPPALL,XS,YS,ZS,MW1,MW2,MW3,
     5  FEXPE,FEXPO,FEXPBACK,RLSC,C,CS,CSINV,
     6  NBOX,NBOXES,IZ,LADDR,NLEV,IWORK,LW,NINIRE,SIZE0)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
c
c  purpose:
c
c    the main subroutine of the new fmm, based on multipole and
c      exponential expansions. two passes are executed. in the
c      first pass, multipole expansions for all boxes at all
c      levels are computed. in the second pass, interactions are
c      computed at successively finer levels.
c
c    this is the uniform code.
c
c  on input:
c
c    iflag: free space or periodic boundary condition.
c           iflag = 0 for free-space problems
c           iflag = 1 for periodic problems
c    nlev: the number of levels in the uniform tree structure.
c    xat,yat,zat(natoms): the location of the particles,
c    natoms: total number of atoms.
c    charge(natoms): the charge each particle carries.
c    nterms: the number of terms in the multipole/local expansion.
c    nlambs: the number of terms in the exponential expansion.
c    nexptot: the total number of exponential expansion terms.
c    nexptotp: the number of exponential expansions.
c    nthmax: the max number of fourier terms in the exponential
c      expansion.
c
c  adaptive tree structure:
c    nbox: the maximum number of points allowed in a childless box
c    nboxes: the total number of boxes created
c    iz(natoms): the integer array addressing the particles in all
c      boxes.
c    laddr(2,200): an integer array (2,nlev), describing the
c      numbers of boxes on various levels of subdivision, so that
c      the first box on level (i-1) has sequence number laddr(1,i),
c      and there are laddr(2,i) boxes on level i-1
c    nlev: the number of levels in the adaptive tree structure.
c      the maximim number possible is 200.
c    iwork(lw) - the array containing all tables describing boxes,
c      lists, etc. it is a link-list (for the most part), and can
c      only be accessed via the entries d3mgetb, d3mgetl, d3mlinfo,
c      of this subroutine.
c    lw - the amount of memory in the array w (in integer*4 elements)
c    ninire: the ratio for storing real*8 and integer *4. this is
c      required as only iwork is used for storage. note: for ifort
c      compiler: real *8 = 8 bytes. integer *4 = 4 bytes. so ninire=2.
c
c  precomputed tables:
c    dc(0:nterms,0:nterms,0:nterms): precomputed array containing
c      coefficients for the local translation along the z axis.
c      this is precomputed by the subroutine lcshftcoef().
c    rdplus(0:nterms,0:nterms,-nterms:nterms): rotation matrix,
c      y to z, see subroutine rotgen<- fstrtn().
c    rdminus(0:nterms,0:nterms,-nterms:nterms): similar to rdplus,
c      z to x.
c    rdsq3(0:nterms,0:nterms,-nterms:nterms): similar to rdplus,
c      shifts the multipole and local expansions, +z direction.
c    rdmsq3(0:nterms,0:nterms,-nterms:nterms): similar to rdsq3,
c      -z direction.
c    whts(nlambs): the weights for the plane wave expansion.
c    rlams(nlambs): the nodes for the plane wave expansion.
c    numfour(nlambs): number of fourier modes in the expansion.
c    numphys(nlambs): number of modes in the plane wave expansion.
c    c(0:60,0:60),cs(0:60,0:60),csinv(0:60,0:60): precomputed
c       vectors for factorials.
c
c  on output:
c    pot(natoms): the potential at different particle locations.
c    field(3,natoms) : the force field.
c
c  variables
c    xat2,yat2,zat2,charg2(nbox): for generating multipole and
c      local expansions.
c    mpole(0:nterms+1,0:nterms+1,nboxes): multipole expansions for
c      all boxes.
c    local(0:nterms+1,0:nterms+1,nboxes): local expansions for
c      all boxes.
c    lexp1(mnexptotp,nboxes): exponential expansions in the
c      first direction.
c    lexp2(mnexptotp,nboxes): similar to lexp1, for the second direction.
c    icnt,icnt2: workspace array for assigning particles.
c    mexpf1(nexptot): used for exponential expansion, for mexpup.
c    mexpf2(nexptot): similar to mexpf1, for mexpdown.
c    mexpp1(mnexptotp): similar to mexpf1, for mexpuphys.
c    mexpp2(mnexptotp): similar to mexpf1, for mexpdnphys.
c    mexppall(mnexptotp,16): used for exponential expansions. the
c      expansions are first merged and then translated to different
c      locations.
c    xs(3,mnexptotp): stores the diagonal translation operators when
c      shifting the exponential expansion. this is complex *16.
c    ys(3,mnexptotp): similar to xs.
c    zs(3,mnexptotp): similar to xs, but this is real *8.
c    mw1(0:nterms,0:nterms): temporary working space for storing multipole
c      and local expansions.
c    mw2(0:nterms,0:nterms): similar to mw1.
c    mw3(0:nterms,0:nterms): similar to mw1, for "point and shoot" technique.
c    fexpe(?),fexpo(?),fexpback(?): used for merging exponential expansions.
c      note that the size of these vectors changes depending on beta and
c      accuracy requirements. therefore the memory must be allocated
c      correctly. the following is an estimate of the size.
c      size of fexpe, fexpo, fexpback = 40000 for nlambs = 39
c      size of fexpe, fexpo, fexpback = 15000 for nlambs = 30
c      size of fexpe, fexpo, fexpback =  4000 for nlambs = 20
c      size of fexpe, fexpo, fexpback =   400 for nlambs = 10
c    rlsc(0:nterms,0:nterms,nlambs): stores p_n^m for different lambda_k.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
c
      IMPLICIT NONE
c
      INTEGER *4 IFLAG,NTERMS,NATOMS,NLAMBS
      INTEGER *4 NEXPTOT,NEXPTOTP,NTHMAX
      INTEGER *4 NBOX,NBOXES,IZ(NATOMS),LADDR(2,200),NLEV
      INTEGER *4 NUMPHYS(1), NUMFOUR(1)
      INTEGER *4 LW,IWORK(LW),NINIRE
      INTEGER *4 IEXP1(NBOXES),IEXP2(NBOXES)
c
      REAL *8 ZAT(3,1),CHARGE(1)
      REAL *8 ZAT2(3,1),CHARG2(1)
      REAL *8 POT(1),FIELD(3,1)
      REAL *8 RLAMS(1)
      REAL *8 C(0:60,0:60),CS(0:60,0:60),CSINV(0:60,0:60)
      REAL *8 ZS(3,NEXPTOTP)
      REAL *8 RLSC(NLAMBS,0:NTERMS,0:NTERMS)
      REAL *8 RDPLUS(0:NTERMS,0:NTERMS,-NTERMS:NTERMS)
      REAL *8 RDMINUS(0:NTERMS,0:NTERMS,-NTERMS:NTERMS)
      REAL *8 RDSQ3(0:NTERMS,0:NTERMS,-NTERMS:NTERMS)
      REAL *8 RDMSQ3(0:NTERMS,0:NTERMS,-NTERMS:NTERMS)
      REAL *8 DC(0:2*NTERMS,0:2*NTERMS)
      REAL *8 WHTS(NLAMBS)
      REAL *8 SIZE0
c
      COMPLEX *16 MPOLE(0:NTERMS,0:NTERMS,*)
      COMPLEX *16 LOCAL(0:NTERMS,0:NTERMS,*)
      COMPLEX *16 LEXP1(NEXPTOTP,*)
      COMPLEX *16 LEXP2(NEXPTOTP,*)
      COMPLEX *16 MEXPF1(NEXPTOT)
      COMPLEX *16 MEXPF2(NEXPTOT)
      COMPLEX *16 MEXPP1(NEXPTOTP)
      COMPLEX *16 MEXPP2(NEXPTOTP)
      COMPLEX *16 MEXPPALL(NEXPTOTP,16)
      COMPLEX *16 XS(3,NEXPTOTP)
      COMPLEX *16 YS(3,NEXPTOTP)
      COMPLEX *16 FEXPE(1)
      COMPLEX *16 FEXPO(1)
      COMPLEX *16 FEXPBACK(1)
      COMPLEX *16 MW1(0:NTERMS,0:NTERMS)
      COMPLEX *16 MW2(0:NTERMS,0:NTERMS)
      COMPLEX *16 MW3(0:NTERMS,0:NTERMS)
c
c-----local variables
c
      INTEGER *4 I,J,K,KK,IBOX,JBOX,IAT,JAT,IFL
      INTEGER *4 BOX(15),BOX2(15),NKIDS,NKIDS2
      INTEGER *4 NINBOX,KAT,ILEV,MYCHILD
      INTEGER *4 NLIST,LIST(2000),MYLIST,LUSED
      INTEGER *4 JJ
      INTEGER *4 IEND,ISTART2,IEND2,NTERMS2
      INTEGER *4 IER
      INTEGER *4 IEXP(16)
c
      REAL *8 CENTER(3),CENTER2(3),SIZE,SIZE2
c
      INTEGER *4 I1(36),N1,IX1(36),IY1(36)
      INTEGER *4 I2(16),N2,IX2(16),IY2(16)
      INTEGER *4 I3(4),N3,IX3(4),IY3(4)
      INTEGER *4 I4(4),N4,IX4(4),IY4(4)
      INTEGER *4 I5(4),N5,IX5(4),IY5(4)
      INTEGER *4 I6(4),N6,IX6(4),IY6(4)
      INTEGER *4 I7(4),N7,IX7(4),IY7(4)
      INTEGER *4 I8(4),N8,IX8(4),IY8(4)
      INTEGER *4 ISTART,IOFF
c
      REAL *8 P(1000)
      REAL *8 ZERO
      REAL *8 SCALE(20),SHIFT(3)
c
      COMPLEX *16 IMAG
c
c-----for timing purposes.
c
      REAL *8 TIMEA,TIMEB
      REAL *8 TIME_INIT,TIME_UPPASS,TIME_S2,TIME_S3
      REAL *8 TIME_S4,TIME_S5,TIME_S6,TIME_S7UD
      REAL *8 TIME_S7NS,TIME_S7EW,TIME_S7ALL
c
c-----functions called.
c
      REAL *8 SECOND
c
      DATA IMAG/(0.0D0,1.0D0)/
      DATA ZERO/0.0D0/
c
c-----initialize time counters:
c
      TIME_INIT = 0
      TIME_UPPASS = 0
      TIME_S2 = 0
      TIME_S3 = 0
      TIME_S4 = 0
      TIME_S5 = 0
      TIME_S6 = 0
      TIME_S7UD = 0
      TIME_S7NS = 0
      TIME_S7EW = 0
      TIME_S7ALL = 0
c
c-----initialize multipole and local expansions to zero.
c
      TIMEB = SECOND()
      DO I = 1,NBOXES
        DO J = 0,NTERMS
          DO K = 0,NTERMS
            MPOLE(J,K,I) = ZERO
            LOCAL(J,K,I) = ZERO
          ENDDO
        ENDDO
      ENDDO
c
c-----initialize the potential and field at all locations.
c
      DO I=1,NATOMS
        POT(I)=0.0D0
        FIELD(1,I)=0.0D0
        FIELD(2,I)=0.0D0
        FIELD(3,I)=0.0D0
      ENDDO
c
c-----set scale factors for all levels.
c
      DO 100 I = 1,NLEV+1
        SCALE(I) = 2.0D0**(I-1)/SIZE0
100   CONTINUE
c
      NTERMS2=NTERMS*NTERMS
c
      TIMEA=SECOND()
      TIME_INIT=TIME_INIT+TIMEA-TIMEB
      CALL PRIN2(' init time is *',TIME_INIT,1)
c
c=============================================================
c-----begin upward pass
c======================================================================
c     step 1: form multipole expansions at all levels
c       if childless, then form multipole expansion.
c       if has children, then merge children's expansions.
c       the step will loop over all boxes.
c
c       note: we combined steps 1 and 2 as described in cheng's
c         3d adaptive fmm paper.
c=============================================================
c
      CALL PRINF(' beginning upward pass *',NLEV,0)
c
      ILEV=NLEV+1
c
      TIMEB=SECOND()
      DO 300 IBOX=NBOXES,1,-1
c
c-------get the box information.
c
        CALL D3MGETB(IER,IBOX,BOX,NKIDS,CENTER,SIZE,IWORK)
c
        IF (NKIDS.EQ.0) THEN
c
c======================================================================
c---------case 1: construct the multipole expansion for childless boxes.
c======================================================================
c
          NINBOX=BOX(15)
          DO 200 K = 1,BOX(15)
            KAT=IZ(BOX(14)+K-1)
            ZAT2(1,K) = ZAT(1,KAT)
            ZAT2(2,K) = ZAT(2,KAT)
            ZAT2(3,K) = ZAT(3,KAT)
            CHARG2(K) = CHARGE(KAT)
200       CONTINUE
          CALL FORMMP(CENTER,ZAT2,CHARG2,
     1      NINBOX,MPOLE(0,0,IBOX),NTERMS,SCALE(BOX(1)+1),P,C)
        ELSE
c
c=======================================================================
c---------case 2: the box has kids, so merge children's multipole expansions
c           to current box.
c======================================================================
c
c---------now loop over -z children (if exist).
c
          DO K=6,9
            MYCHILD=BOX(K)
            IF (MYCHILD.GT.0) THEN
c
c-------------the child's multipole expansion is already created.
c
              IF (K.EQ.6) THEN
                IFL=3
              ELSEIF (K.EQ.7) THEN
                IFL=4
              ELSEIF (K.EQ.8) THEN
                IFL=2
              ELSEIF (K.EQ.9) THEN
                IFL=1
              ENDIF
c
c-------------shift the child's multipole expansion.
c
              CALL MPSHIFT(IFL,MPOLE(0,0,MYCHILD),MW1,MW2,
     1          NTERMS,DC,RDSQ3,SCALE(BOX(1)+2),SCALE(BOX(1)+1))
              CALL ADDEXP(MW1,MPOLE(0,0,IBOX),NTERMS)
            ENDIF
          ENDDO
c
c---------now loop over +z children (if exist).
c
          DO K=10,13
            MYCHILD=BOX(K)
            IF (MYCHILD.GT.0) THEN
c
c-------------the child's multipole expansion is already created.
c
              IF (K.EQ.10) THEN
                IFL=3
              ELSEIF (K.EQ.11) THEN
                IFL=4
              ELSEIF (K.EQ.12) THEN
                IFL=2
              ELSEIF (K.EQ.13) THEN
                IFL=1
              ENDIF
              CALL MPSHIFT(IFL,MPOLE(0,0,MYCHILD),MW1,MW2,
     1          NTERMS,DC,RDMSQ3,SCALE(BOX(1)+2),SCALE(BOX(1)+1))
              CALL ADDEXP(MW1,MPOLE(0,0,IBOX),NTERMS)
            ENDIF
          ENDDO
        ENDIF
c
300   CONTINUE
c
      TIMEA=SECOND()
      TIME_UPPASS=TIME_UPPASS+TIMEA-TIMEB
c
c=======================================================================
c
c     upward pass complete: all multipole expansions are available.
c
c=======================================================================
c
c     begin downward pass
c     at each level, do the following steps:
c
c     step 2: for all boxes, compute interactions due to list 4.
c     step 3: shift parent's local to its children, if has kids.
c     step 4: compute interactions due to direct interactions, if childless.
c     step 5: compute interactions due to list 3, if childless.
c
c     step 6: evaluate the local expansion.
c
c     step 7: compute interactions due to interaction list for children's
c             level.
c
c
c     note that steps 2, 3, 4 and 5 are combined.
c     reference: cheng's 3d adaptive paper.
c
c======================================================================
c
      DO 2000 ILEV=1,NLEV+1
c
        CALL PRINF(' downward pass, ilev = *',ILEV,1)
c
c-------loop over all boxes in this level.
c
        ISTART=LADDR(1,ILEV)
        IEND=LADDR(1,ILEV)+LADDR(2,ILEV)-1
c
        DO 3000 IBOX=ISTART,IEND
c
c---------get the box information.
c
          CALL D3MGETB(IER,IBOX,BOX,NKIDS,CENTER,SIZE,IWORK)
c
c======================================================================
c--step 2: for all boxes, process list #4.
c======================================================================
c
          TIMEB=SECOND()
          CALL D3MGETLIST(IER,IBOX,4,LIST,NLIST,IWORK)
c
          IF (BOX(15) .GT. NTERMS2) THEN
c
c-----------box # ibox has many particles. the most efficient way
c               is to form the local expansion for current box using
c               big box's particle information.
c

            DO JJ=1,NLIST
              JBOX=LIST(JJ)
c
              CALL D3MGETB(IER,JBOX,BOX2,NKIDS2,CENTER2,SIZE,IWORK)
c
              NINBOX=BOX2(15)
              DO 210 K = 1,BOX2(15)
                KAT=IZ(BOX2(14)+K-1)
                ZAT2(1,K) = ZAT(1,KAT)
                ZAT2(2,K) = ZAT(2,KAT)
                ZAT2(3,K) = ZAT(3,KAT)
                CHARG2(K)  = CHARGE(KAT)
210           CONTINUE
              CALL FORMLC(CENTER,ZAT2,CHARG2,NINBOX,
     1          LOCAL(0,0,IBOX),NTERMS,SCALE(BOX(1)+1),P,C)
            ENDDO
          ELSE
c
c-----------direct interaction is more efficient.
c
            DO JJ=1,NLIST
              JBOX=LIST(JJ)
c
              CALL D3MGETB(IER,JBOX,BOX2,NKIDS2,CENTER2,SIZE2,IWORK)
c
              CALL BRNBRF(NATOMS,ZAT,CHARGE,POT,FIELD,BOX(14),
     1          BOX(15),BOX2(14),BOX2(15),IZ)
            ENDDO
          ENDIF
c
          TIMEA=SECOND()
          TIME_S2=TIME_S2+TIMEA-TIMEB
c
          IF (NKIDS.NE.0) THEN
c
c
c======================================================================
c--step 3: box with children, send out parent's information.
c======================================================================
c             shift its local expansion to its children. first loop
c             over -z children (if exist), and shift the local expansion.
c
            TIMEB=SECOND()
c
            DO K=6,9
              MYCHILD=BOX(K)
              IF (MYCHILD.GT.0) THEN
                IF (K.EQ.6) THEN
                  IFL=1
                ELSEIF (K.EQ.7) THEN
                  IFL=2
                ELSEIF (K.EQ.8) THEN
                  IFL=4
                ELSEIF (K.EQ.9) THEN
                  IFL=3
                ENDIF
c
c---------------shift the local expansion to child.
c
                CALL LCSHIFT(IFL,LOCAL(0,0,IBOX),MW1,MW2,
     1            NTERMS,DC,RDSQ3,SCALE(BOX(1)+1),SCALE(BOX(1)+2))
                CALL ADDEXP(MW1,LOCAL(0,0,MYCHILD),NTERMS)
              ENDIF
            ENDDO
c
c-----------now loop over +z children (if exist) and shift the local to child.
c
            DO K=10,13
              MYCHILD=BOX(K)
              IF (MYCHILD.GT.0) THEN
                IF (K.EQ.10) THEN
                  IFL=1
                ELSEIF (K.EQ.11) THEN
                  IFL=2
                ELSEIF (K.EQ.12) THEN
                  IFL=4
                ELSEIF (K.EQ.13) THEN
                  IFL=3
                ENDIF
                CALL LCSHIFT(IFL,LOCAL(0,0,IBOX),MW1,MW2,
     1            NTERMS,DC,RDMSQ3,SCALE(BOX(1)+1),SCALE(BOX(1)+2))
                CALL ADDEXP(MW1,LOCAL(0,0,MYCHILD),NTERMS)
              ENDIF
            ENDDO
c
            TIMEA=SECOND()
            TIME_S3=TIME_S3+TIMEA-TIMEB
          ELSE
c
c======================================================================
c--step 4: childless box, collect information from direct
c             interaction list. get list 1 for current box.
c======================================================================
c
            TIMEB=SECOND()
c
            CALL D3MGETLIST(IER,IBOX,1,LIST,NLIST,IWORK)
            DO JJ=1,NLIST
              JBOX=LIST(JJ)
              CALL D3MGETB(IER,JBOX,BOX2,NKIDS2,CENTER2,SIZE2,IWORK)
              CALL BRNBRF(NATOMS,ZAT,CHARGE,POT,FIELD,BOX(14),
     1          BOX(15),BOX2(14),BOX2(15),IZ)
            ENDDO
c
c-----------interaction with self.
c
            CALL BRNBRF(NATOMS,ZAT,CHARGE,POT,FIELD,BOX(14),
     1        BOX(15),BOX(14),BOX(15),IZ)
c
            TIMEA=SECOND()
            TIME_S4=TIME_S4+TIMEA-TIMEB
c
c======================================================
c======================================================================
c--step 5: childless box, collect information from
c             smaller well-separated box in list #3.
c======================================================================
c
c           note: both the particle information and the
c             multipole expansion are available, so choose the
c             most efficient one.
c
c           first get the list.
c
            TIMEB=SECOND()
c
            CALL D3MGETLIST(IER,IBOX,3,LIST,NLIST,IWORK)
            DO JJ=1,NLIST
              JBOX=LIST(JJ)
              CALL D3MGETB(IER,JBOX,BOX2,NKIDS2,CENTER2,SIZE2,IWORK)
              IF (BOX2(15).GT. NTERMS2) THEN
c
c---------------too many particles, evaluating multipole expansion
c                 is more efficient.
c
                DO KK=1,BOX(15)
                  IAT=IZ(BOX(14)+KK-1)
                  CALL BRMPEV(MPOLE(0,0,JBOX),CENTER2,ZAT(1,IAT),
     1              NTERMS,POT(IAT),FIELD(1,IAT),SCALE(BOX2(1)+1),P,
     2              C,CS,CSINV)
                ENDDO
              ELSE
c
c---------------direct particle particle interaction.
c
                CALL BRNBRF(NATOMS,ZAT,CHARGE,POT,FIELD,BOX(14),
     1            BOX(15),BOX2(14),BOX2(15),IZ)
              ENDIF
            ENDDO
c
            TIMEA=SECOND()
            TIME_S5=TIME_S5+TIMEA-TIMEB
c
          ENDIF
3000    ENDDO
c
c======================================================================
c--step 6: evaluate the local expansion for childless boxes.
c======================================================================
c
        TIMEB=SECOND()
c
        DO 6000 IBOX=ISTART,IEND
c
c---------get the box information.
c
          CALL D3MGETB(IER,IBOX,BOX,NKIDS,CENTER,SIZE,IWORK)
c
          IF (NKIDS.EQ.0) THEN
            DO 6100 J = 1,BOX(15)
              IAT = IZ(BOX(14)+J-1)
              CALL BRTAEV(LOCAL(0,0,IBOX),CENTER,ZAT(1,IAT),
     1          NTERMS,POT(IAT),FIELD(1,IAT),SCALE(BOX(1)+1),P,
     2          C,CS,CSINV)
6100        CONTINUE
          ENDIF
6000    ENDDO
c
        TIMEA=SECOND()
        TIME_S6=TIME_S6+TIMEA-TIMEB
c
c======================================================================
c--step 7: for current box's children, compute interaction in
c       the interaction list, and send resulting local expansions
c       to children.
c======================================================================
c
c       not sure here if the merging then shifting idea is good or
c       not, for nonadaptive, this is definitely good, but for
c       adaptive...
c
c       i am implementing this merging version first.
c
c-------initialize the exponential expansions.
c
        ISTART2=LADDR(1,ILEV+1)
        IEND2=LADDR(1,ILEV+1)+LADDR(2,ILEV+1)-1
c
c-------process the up-down list.
c
        TIMEB=SECOND()
c
        DO IBOX=ISTART2,IEND2
          IEXP1(IBOX)=0
          IEXP2(IBOX)=0
          DO JJ = 1,NEXPTOTP
            LEXP1(JJ,IBOX) = 0.0D0
            LEXP2(JJ,IBOX) = 0.0D0
          ENDDO
        ENDDO
c
        DO 7200 IBOX = ISTART,IEND
c
c---------get the box information.
c
          CALL D3MGETB(IER,IBOX,BOX,NKIDS,CENTER,SIZE,IWORK)
c
c---------for childless box, do nothing.
c
          IF (NKIDS.EQ.0) GOTO 7200
c
c---------process the interaction list for current box's children.
c
          CALL MKUDEXP(IBOX,BOX,NTERMS,MPOLE,RLAMS,NLAMBS,
     1      NUMFOUR,NUMPHYS,NTHMAX,NEXPTOT,NEXPTOTP,MEXPF1,
     2      MEXPF2,MEXPP1,MEXPP2,MEXPPALL(1,1),MEXPPALL(1,2),
     3      MEXPPALL(1,3),MEXPPALL(1,4),IEXP,XS,YS,ZS,
     4      FEXPE,FEXPO,RLSC)
c
          CALL MKUPLIST(IBOX,BOX,CENTER,SIZE,IWORK,
     1      I1,N1,IX1,IY1,I2,N2,IX2,IY2)
          CALL PROCESSUP(SCALE(BOX(1)+2),LEXP1,I1,N1,IX1,IY1,I2,N2,
     1      IX2,IY2,MEXPPALL(1,1),MEXPPALL(1,2),XS,YS,ZS,NEXPTOTP,
     2      IEXP,IEXP1)
c
          CALL MKDNLIST(IBOX,BOX,CENTER,SIZE,IWORK,
     1      I1,N1,IX1,IY1,I2,N2,IX2,IY2)
          CALL PROCESSDN(SCALE(BOX(1)+2),LEXP2,I1,N1,IX1,IY1,I2,N2,
     1      IX2,IY2,MEXPPALL(1,3),MEXPPALL(1,4),XS,YS,ZS,NEXPTOTP,
     2      IEXP,IEXP2)
7200    CONTINUE
c
c-------convert lexp1 and lexp2 to local expansion.
c
        DO 7300 IBOX = ISTART2,IEND2
          IF (IEXP1(IBOX).GT.0) THEN
            CALL PHYSTOF(MEXPF1,NLAMBS,RLAMS,NUMFOUR,NUMPHYS,
     1        NTHMAX,LEXP1(1,IBOX),FEXPBACK)
          ENDIF
          IF (IEXP2(IBOX).GT.0) THEN
            CALL PHYSTOF(MEXPF2,NLAMBS,RLAMS,NUMFOUR,NUMPHYS,
     1        NTHMAX,LEXP2(1,IBOX),FEXPBACK)
          ENDIF
          IF (IEXP1(IBOX).GT.0 .OR. IEXP2(IBOX).GT.0 ) THEN
            CALL EXPTOLOCAL(MW1,NTERMS,RLAMS,WHTS,NLAMBS,
     1        NUMFOUR,NTHMAX,NEXPTOT,IEXP1(IBOX),MEXPF1,
     2        IEXP2(IBOX),MEXPF2,CS)
            CALL ADDEXP(MW1,LOCAL(0,0,IBOX),NTERMS)
          ENDIF
7300    CONTINUE
c
        TIMEA=SECOND()
        TIME_S7UD=TIME_S7UD+TIMEA-TIMEB
c
c-------next process north and south lists
c
        TIMEB=SECOND()
c
        DO IBOX=ISTART2,IEND2
          IEXP1(IBOX)=0
          IEXP2(IBOX)=0
          DO JJ = 1,NEXPTOTP
            LEXP1(JJ,IBOX) = 0.0D0
            LEXP2(JJ,IBOX) = 0.0D0
          ENDDO
        ENDDO
c
        DO 7400 IBOX = ISTART,IEND
c
c---------get the box information.
c
          CALL D3MGETB(IER,IBOX,BOX,NKIDS,CENTER,SIZE,IWORK)
c
c---------for childless box, do nothing.
c
          IF (NKIDS.EQ.0) GOTO 7400
c
c---------process the interaction list for current box's children.
c
          CALL MKNSEXP(IBOX,BOX,NTERMS,MPOLE,MW1,MW2,
     1      RLAMS,NLAMBS,NUMFOUR,NUMPHYS,NTHMAX,NEXPTOT,NEXPTOTP,
     2      MEXPF1,MEXPF2,MEXPP1,MEXPP2,RDMINUS,MEXPPALL(1,1),
     3      MEXPPALL(1,2),MEXPPALL(1,3),MEXPPALL(1,4),
     4      MEXPPALL(1,5),MEXPPALL(1,6),MEXPPALL(1,7),
     5      MEXPPALL(1,8),IEXP,XS,YS,ZS,FEXPE,FEXPO,RLSC)
         CALL MKNOLIST(IBOX,BOX,CENTER,SIZE,IWORK,
     1      I1,N1,IX1,IY1,I2,N2,IX2,IY2,I3,N3,IX3,IY3,I4,N4,IX4,IY4)
          CALL PROCESSNO(SCALE(BOX(1)+2),LEXP1,I1,N1,IX1,IY1,I2,N2,
     1      IX2,IY2,I3,N3,IX3,IY3,I4,N4,IX4,IY4,
     2      MEXPPALL(1,1),MEXPPALL(1,2),MEXPPALL(1,3),
     3      MEXPPALL(1,4),XS,YS,ZS,NEXPTOTP,IEXP,IEXP1)
c
          CALL MKSOLIST(IBOX,BOX,CENTER,SIZE,IWORK,
     1      I1,N1,IX1,IY1,I2,N2,IX2,IY2,I3,N3,IX3,IY3,I4,N4,IX4,IY4)
          CALL PROCESSSO(SCALE(BOX(1)+2),LEXP2,I1,N1,IX1,IY1,I2,N2,
     1      IX2,IY2,I3,N3,IX3,IY3,I4,N4,IX4,IY4,
     2      MEXPPALL(1,5),MEXPPALL(1,6),MEXPPALL(1,7),
     3      MEXPPALL(1,8),XS,YS,ZS,NEXPTOTP,IEXP,IEXP2)
7400    CONTINUE
c
c-------convert lexp1 and lexp2 to local multipole expansion.
c
        DO 7500 IBOX = ISTART2,IEND2
c
          IF (IEXP1(IBOX).GT.0) THEN
            CALL PHYSTOF(MEXPF1,NLAMBS,RLAMS,NUMFOUR,NUMPHYS,
     1        NTHMAX,LEXP1(1,IBOX),FEXPBACK)
          ENDIF
          IF (IEXP2(IBOX).GT.0) THEN
            CALL PHYSTOF(MEXPF2,NLAMBS,RLAMS,NUMFOUR,NUMPHYS,
     1        NTHMAX,LEXP2(1,IBOX),FEXPBACK)
          ENDIF
          IF (IEXP1(IBOX).GT.0 .OR. IEXP2(IBOX).GT.0) THEN
            CALL EXPTOLOCAL(MW1,NTERMS,RLAMS,WHTS,NLAMBS,
     1        NUMFOUR,NTHMAX,NEXPTOT,IEXP1(IBOX),MEXPF1,
     2        IEXP2(IBOX),MEXPF2,CS)
            CALL ROTYTOZ(NTERMS,MW1,MW3,MW2,RDPLUS)
            CALL ADDEXP(MW2,LOCAL(0,0,IBOX),NTERMS)
          ENDIF
7500    CONTINUE
c
        TIMEA=SECOND()
        TIME_S7NS=TIME_S7NS+TIMEA-TIMEB
c
c-------next process east and west lists
c
        TIMEB=SECOND()
c
        DO IBOX=ISTART2,IEND2
          IEXP1(IBOX)=0
          IEXP2(IBOX)=0
          DO JJ = 1,NEXPTOTP
            LEXP1(JJ,IBOX) = 0.0D0
            LEXP2(JJ,IBOX) = 0.0D0
          ENDDO
        ENDDO
c
        DO 7600 IBOX = ISTART,IEND
c
c---------get the box information.
c
          CALL D3MGETB(IER,IBOX,BOX,NKIDS,CENTER,SIZE,IWORK)
c
c---------for childless box, do nothing.
c
          IF (NKIDS.EQ.0) GOTO 7600
c
c---------process the interaction list for current box's children.
c
          CALL MKEWEXP(IBOX,BOX,NTERMS,MPOLE,MW1,RLAMS,
     1      NLAMBS,NUMFOUR,NUMPHYS,NTHMAX,NEXPTOT,NEXPTOTP,
     2      MEXPF1,MEXPF2,MEXPP1,MEXPP2,RDPLUS,MEXPPALL(1,1),
     3      MEXPPALL(1,2),MEXPPALL(1,3),MEXPPALL(1,4),
     4      MEXPPALL(1,5),MEXPPALL(1,6),MEXPPALL(1,7),
     5      MEXPPALL(1,8),MEXPPALL(1,9),MEXPPALL(1,10),
     6      MEXPPALL(1,11),MEXPPALL(1,12),MEXPPALL(1,13),
     7      MEXPPALL(1,14),MEXPPALL(1,15),MEXPPALL(1,16),
     8      IEXP,XS,YS,ZS,FEXPE,FEXPO,RLSC)
          CALL MKEALIST(IBOX,BOX,CENTER,SIZE,IWORK,
     1      I1,N1,IX1,IY1,I2,N2,IX2,IY2,I3,N3,IX3,IY3,I4,N4,
     2      IX4,IY4,I5,N5,IX5,IY5,I6,N6,IX6,IY6,I7,N7,IX7,IY7,
     3      I8,N8,IX8,IY8)
          CALL PROCESSEA(SCALE(BOX(1)+2),LEXP1,I1,N1,IX1,IY1,I2,N2,
     1      IX2,IY2,I3,N3,IX3,IY3,I4,N4,IX4,IY4,I5,N5,IX5,IY5,
     2      I6,N6,IX6,IY6,I7,N7,IX7,IY7,I8,N8,IX8,IY8,
     3      MEXPPALL(1,1),MEXPPALL(1,2),MEXPPALL(1,3),
     4      MEXPPALL(1,4),MEXPPALL(1,5),MEXPPALL(1,6),
     5      MEXPPALL(1,7),MEXPPALL(1,8),XS,YS,ZS,NEXPTOTP,
     6      IEXP,IEXP1)
c
          CALL MKWELIST(IBOX,BOX,CENTER,SIZE,IWORK,I1,N1,IX1,
     1      IY1,I2,N2,IX2,IY2,I3,N3,IX3,IY3,I4,N4,IX4,IY4,
     2      I5,N5,IX5,IY5,I6,N6,IX6,IY6,I7,N7,IX7,IY7,
     3      I8,N8,IX8,IY8)
          CALL PROCESSWE(SCALE(BOX(1)+2),LEXP2,I1,N1,IX1,IY1,I2,N2,
     1      IX2,IY2,I3,N3,IX3,IY3,I4,N4,IX4,IY4,I5,N5,IX5,IY5,
     2      I6,N6,IX6,IY6,I7,N7,IX7,IY7,I8,N8,IX8,IY8,
     3      MEXPPALL(1,9),MEXPPALL(1,10),MEXPPALL(1,11),
     4      MEXPPALL(1,12),MEXPPALL(1,13),MEXPPALL(1,14),
     5      MEXPPALL(1,15),MEXPPALL(1,16),XS,YS,ZS,NEXPTOTP,
     6      IEXP,IEXP2)
7600    CONTINUE
c
c-------convert lexp1 and lexp2 to local multipole expansion.
c
        DO 7700 IBOX = ISTART2,IEND2
c
          IF (IEXP1(IBOX).GT.0) THEN
            CALL PHYSTOF(MEXPF1,NLAMBS,RLAMS,NUMFOUR,NUMPHYS,
     1        NTHMAX,LEXP1(1,IBOX),FEXPBACK)
          ENDIF
          IF (IEXP2(IBOX).GT.0) THEN
            CALL PHYSTOF(MEXPF2,NLAMBS,RLAMS,NUMFOUR,NUMPHYS,
     1        NTHMAX,LEXP2(1,IBOX),FEXPBACK)
          ENDIF
          IF (IEXP1(IBOX).GT.0 .OR. IEXP2(IBOX).GT.0) THEN
            CALL EXPTOLOCAL(MW1,NTERMS,RLAMS,WHTS,NLAMBS,
     1        NUMFOUR,NTHMAX,NEXPTOT,IEXP1(IBOX),MEXPF1,
     2        IEXP2(IBOX),MEXPF2,CS)
            CALL ROTZTOX(NTERMS,MW1,MW2,RDMINUS)
            CALL ADDEXP(MW2,LOCAL(0,0,IBOX),NTERMS)
          ENDIF
7700    CONTINUE
c
        TIMEA=SECOND()
        TIME_S7EW=TIME_S7EW+TIMEA-TIMEB
c
c======================================================================
c     downward pass complete
c     this include all the interactions in the adaptive interactions.
c======================================================================
c
2000  ENDDO
c
c-----output cpu time for different steps.
c
      CALL PRIN2(' init time (assign etc.) is *',TIME_INIT,1)
      CALL PRIN2(' upward pass time is *',TIME_UPPASS,1)
      CALL PRIN2(' step 2 for list 4 is *',TIME_S2,1)
      CALL PRIN2(' step 3 for tata is *',TIME_S3,1)
      CALL PRIN2(' step 4 for list 1 is *',TIME_S4,1)
      CALL PRIN2(' step 5 for list 3 is *',TIME_S5,1)
      CALL PRIN2(' step 6 for taev is *',TIME_S6,1)
      CALL PRIN2(' step 7 for ud is *',TIME_S7UD,1)
      CALL PRIN2(' step 7 for ns is *',TIME_S7NS,1)
      CALL PRIN2(' step 7 for ew is *',TIME_S7EW,1)
      TIME_S7ALL=TIME_S7EW+TIME_S7NS+TIME_S7UD
      CALL PRIN2(' step 7 for all is *',TIME_S7ALL,1)
c
      TIMEB=TIME_UPPASS+TIME_INIT+TIME_S2+TIME_S3+TIME_S4+TIME_S5
     1  +TIME_S6+TIME_S7ALL
      CALL PRIN2(' total time is *',TIMEB,1)
c
      WRITE(11,*)'n     nlev   ntemrs  nboxes   uppass  step2   step3
     1  STEP4    STEP5    STEP6    STEP7       TOTAL'
      WRITE(11,1221)NATOMS,NLEV,NTERMS,NBOXES,TIME_UPPASS,TIME_S2,
     1  TIME_S3,TIME_S4,TIME_S5,TIME_S6,TIME_S7ALL,TIMEB
1221  FORMAT(I7,2X,I2,4X,I2,4X,I6,1X,8F8.2)
c
      RETURN
      END
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
      SUBROUTINE BRNBRF(NATOMS,ZAT,CHARGE,POT,FIELD,ISTART,
     1  IINBOX,JSTART,JINBOX,IZ)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
c
c  purpose: computes interactions directly for two cluster of particles.
c
c  on input :
c    beta: the frequency.
c    natoms: total number of particles.
c    zat: particle locations
c    charge: particle charges
c    istart: starting particle number in box ibox (receiving box).
c    iinbox: total number of particles in box ibox.
c    jstart: starting particle location in box jbox (sending box)
c    jinbox: the total number of particles in box jbox.
c
c  on output:
c    pot : the potential
c
c  note: error checking may be necessary in the future.
c
c  functions called: datan,dsqrt,dexp
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
c
      IMPLICIT NONE
      INTEGER *4 NATOMS,ISTART,IINBOX,JSTART,JINBOX
      INTEGER *4 IZ(NATOMS)
c
      REAL *8 ZAT(3,NATOMS),CHARGE(NATOMS),POT(NATOMS)
      REAL *8 FIELD(3,NATOMS)
c
c-----local variables
c
      INTEGER *4 I,IAT,K,JAT
      REAL *8 PI,RX,RY,RZ,RR,RDIS,RMUL
c
c-----function called
c
      REAL *8 DSQRT,DEXP,DATAN
c
      IF (IINBOX.EQ.0 .OR. JINBOX.EQ.0) THEN
        PRINT *, "problem in direct interaction"
        STOP
      ENDIF
c
      PI = 4.0D0*DATAN(1.0D0)
c
c-----loop through all particles.
c
      DO 3000 I=1,IINBOX
        IAT=IZ(ISTART+I-1)
        DO 2400 K = 1,JINBOX
          JAT = IZ(JSTART+K-1)
c
          IF ( IAT .EQ. JAT) GOTO 2000
          RX = ZAT(1,IAT) - ZAT(1,JAT)
          RY = ZAT(2,IAT) - ZAT(2,JAT)
          RZ = ZAT(3,IAT) - ZAT(3,JAT)
c
          RR =  RX*RX + RY*RY + RZ*RZ
          RDIS = DSQRT(RR)
          POT(IAT) = POT(IAT) + CHARGE(JAT)/RDIS
          RMUL = CHARGE(JAT)/(RDIS*RR)
          FIELD(1,IAT) = FIELD(1,IAT)+RMUL*RX
          FIELD(2,IAT) = FIELD(2,IAT)+RMUL*RY
          FIELD(3,IAT) = FIELD(3,IAT)+RMUL*RZ
2000      CONTINUE
2400    CONTINUE
3000  CONTINUE
c
      RETURN
      END
