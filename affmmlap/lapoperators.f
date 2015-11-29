cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
c
c  this file contains the translation operators for both uniform and
c   adaptive fast multipole laplace solvers.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
      SUBROUTINE FORMMP(X0Y0Z0,ZPARTS,CHARGE,NPARTS,
     1  MPOLE,NTERMS,SCALE,P,C)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
c
c  purpose:
c
c    this subroutine forms the multipole expansion cause by the nparts
c      particles in the box.
c
c  on input:
c
c    x0y0z0: center of the expansion
c    nparts: number of sources
c    zparts(3,nparts): array of coordinates of sources
c    charge(nparts): array of strengths of sources
c    nparts: the total number of particles.
c    nterms: order of desired expansion
c    scale: the scaling factor.
c    c: precomputed numbers.
c
c  on output:
c
c    mpole: coefficients of multipole expansion
c
c  working space :
c
c    p: used for storing the associate legendre polynomials.
c
c  subroutine called : dsqrt(), lgndr()
c  called from : ladapfmm()
c
c  note 1: this subroutine needs the precomputed variables c(,)
c          derived from entry frmini()
c
c       2: the multipole expansion is scaled to avoid over- and
c          under-flow.
c
c       3: only the n=0, ... ,nterms, m=0, ..., nterms coefficients
c          are calculated.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
c
      IMPLICIT NONE
c
      INTEGER *4 NTERMS,NPARTS
      REAL *8 SCALE
      REAL *8 X0Y0Z0(3),ZPARTS(3,NPARTS),CHARGE(NPARTS)
      REAL *8 P(0:NTERMS,0:NTERMS),C(0:60,0:60)
      COMPLEX *16 MPOLE(0:NTERMS,0:NTERMS)
c
c-----local variables.
c
      INTEGER *4 I,L,M,NCALC
      REAL *8 POWERS(0:60)
      REAL *8 PRECIS,D,CP,PROJ
      REAL *8 RX,RY,RZ,RR,CTHETA
      COMPLEX *16 IMAG,EPHI(1:60),CPZ
c
      DATA IMAG/(0.0D0,1.0D0)/
      DATA PRECIS/1.0D-14/
c
      DO 4900 I = 1, NPARTS
        RX = ZPARTS(1,I) - X0Y0Z0(1)
        RY = ZPARTS(2,I) - X0Y0Z0(2)
        RZ = ZPARTS(3,I) - X0Y0Z0(3)
c
c-------compute  distance rr
c         ctheta = cos(theta)
c         ephi(1)=e^(i*phi)
c
        PROJ = RX*RX+RY*RY
        RR = PROJ+RZ*RZ
        PROJ = DSQRT(PROJ)
        D = DSQRT(RR)
c
c-------note: here is a hack. when computing cos(theta) as
c             rz/d, we have to be careful about the possibility
c             of d  being 0 (in which case ctheta is not well
c             defined - we arbitrarily set it to 1.)
c
        IF ( D .LE. PRECIS ) THEN
          CTHETA = 1.0D0
        ELSE
          CTHETA = RZ/D
        ENDIF
        IF ( PROJ .LE. PRECIS*D ) THEN
          EPHI(1) = 1.0D0
        ELSE
          EPHI(1) = RX/PROJ + IMAG*RY/PROJ
        ENDIF
c
c-------create array of powers of r and powers of e^(i*phi).
c
        D = D*SCALE
        POWERS(0) = 1.0D0
c
        DO 4100 L = 1,NTERMS+1
          POWERS(L) = POWERS(L-1)*D
          EPHI(L+1) = EPHI(L)*EPHI(1)
4100    CONTINUE
        MPOLE(0,0) = MPOLE(0,0) + CHARGE(I)
c
c-------compute legendre polynomials of argument cos(theta) = ctheta
c       and add contributions from legendre polynomials
c
        CALL LGNDR(NTERMS,CTHETA,P)
        DO 4300 L = 1,NTERMS
          CP = CHARGE(I)*POWERS(L)*P(L,0)
          MPOLE(L,0) = MPOLE(L,0) + CP
4300    CONTINUE
c
c-------add contributions from associated legendre functions.
c
        DO 4500 L = 1,NTERMS
          DO 4400 M=1,L
            CP = CHARGE(I)*POWERS(L)*C(L,M)*P(L,M)
            MPOLE(L,M) = MPOLE(L,M) + CP*DCONJG(EPHI(M))
4400      CONTINUE
4500    CONTINUE
4900  CONTINUE
c
      RETURN
      END
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
      SUBROUTINE BRTAEV(LOCAL,X0Y0Z0,POINT,NTERMS,RPOT,FIELD,SCALE,
     1             P,C,CS,CSINV)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
c
c  purpose:
c
c    evaluates local expansion at arbitrary point.
c
c  on input:
c
c    scale: a scaling factor to scale the current box to
c      size 1.
c    local: coefficients of local expansion(scaled)
c    x0y0z0: the center of the expansion
c    point: point of evaluation
c    nterms: order of expansion
c    c,cs,csinv: precomputed coefficients.
c
c  on output:
c    rpot: computed potential
c    field: the computed field.
c
c  working array:
c
c    p: work arrays to hold legendre polynomials
c       and associated legendre functions, respectively.
c
c  subroutine called :
c
c  called from : ladapfmm()
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
c
      IMPLICIT NONE
c
      INTEGER *4 NTERMS
      REAL *8 SCALE
      REAL *8 X0Y0Z0(3),P(0:NTERMS,0:NTERMS),POINT(3),RPOT
      REAL *8 FIELD(3),RLOC,RPOTZ
      REAL *8 C(0:60,0:60),CS(0:60,0:60),CSINV(0:60,0:60)
      COMPLEX *16 LOCAL(0:NTERMS,0:NTERMS)
c
c-----local variables.
c
      INTEGER *4 I,L,M,NTP1,NCALC
      REAL *8 RX,RY,RZ,PROJ,RR,CTHETA
      REAL *8 PRECIS,D,DD,CP,POWERS(0:60)
      REAL *8 FIELDTEMP(3)
      COMPLEX *16 ZS1,ZS2,ZS3
      COMPLEX *16 IMAG,EPHI(1:60),CPZ
c
      DATA IMAG/(0.0D0,1.0D0)/
      DATA PRECIS/1.0D-14/
c
c-----functions called.
c
      REAL *8 DREAL,DSQRT,DBLE
c
      RPOTZ = 0.0D0
      ZS1 = 0.0D0
      ZS2 = 0.0D0
      ZS3 = 0.0D0
c
c-----compute relevant functions of spherical coordinates
c     d = distance, ctheta = cos(theta), ephi = exp(i*phi)
c
      RX = POINT(1) - X0Y0Z0(1)
      RY = POINT(2) - X0Y0Z0(2)
      RZ = POINT(3) - X0Y0Z0(3)
      PROJ = RX*RX+RY*RY
      RR = PROJ+RZ*RZ
      PROJ = DSQRT(PROJ)
      D = DSQRT(RR)
      IF (D .LE. PRECIS) THEN
        CTHETA = 0.0D0
      ELSE
        CTHETA = RZ/D
      ENDIF
      IF ( PROJ .LE. PRECIS*D ) THEN
        EPHI(1) = 1.0D0
      ELSE
        EPHI(1) = RX/PROJ + IMAG*RY/PROJ
      ENDIF
c
c-----create array of powers of r and e^(i*m*phi).
c
      D = D*SCALE
      DD = D
      POWERS(0) = 1.0D0
      DO 5600 L = 1,NTERMS+2
        POWERS(L) = DD
        DD = DD*D
        EPHI(L+1) = EPHI(L)*EPHI(1)
5600  CONTINUE
c
c-----compute values of legendre functions
c
      CALL LGNDR(NTERMS,CTHETA,P)
c
c-----add contribution from constant term
c
      RPOT = RPOT + DREAL(LOCAL(0,0))
c
c-----add contributions from legendre polynomials
c
      FIELDTEMP(3)=0.0D0
      DO 5700 L=1,NTERMS
        RLOC=DREAL(LOCAL(L,0))
        CP=RLOC*POWERS(L)*P(L,0)
        RPOT=RPOT + CP
        CP=POWERS(L-1)*P(L-1,0)*CS(L-1,0)
        CPZ=LOCAL(L,1)*(CP*CSINV(L,1))
        ZS2=ZS2+CPZ
        CP=RLOC*CP*CSINV(L,0)
        FIELDTEMP(3) = FIELDTEMP(3)+CP
5700  CONTINUE
c
c-----add contributions from associated legendre functions.
c
      DO 5900 L = 1,NTERMS
c
c-------first compute potential
c
        DO 5800 M=1,L
          CPZ = LOCAL(L,M)*EPHI(M)
          RPOTZ=RPOTZ+DREAL(CPZ)*POWERS(L)*C(L,M)*P(L,M)
5800    CONTINUE
c
c-------then term zs3 for z derivative
c
        DO 5820 M=1,L-1
          ZS3 = ZS3 + LOCAL(L,M)*EPHI(M)*
     1      (POWERS(L-1)*C(L-1,M)*P(L-1,M)*CS(L-1,M)*CSINV(L,M))
5820    CONTINUE
c
c-------then term zs2 for  x,y derivatives.
c
        DO 5840 M=2,L
          ZS2 = ZS2 + LOCAL(L,M)*EPHI(M-1)*(CS(L-1,M-1)*CSINV(L,M)*
     1      POWERS(L-1)*C(L-1,M-1)*P(L-1,M-1))
5840    CONTINUE
c
c-------then term zs1 for  x,y derivatives.
c
        DO 5860 M=0,L-2
          ZS1 = ZS1 + LOCAL(L,M)*EPHI(M+1)*(CS(L-1,M+1)*
     1      CSINV(L,M)*POWERS(L-1)*C(L-1,M+1)*P(L-1,M+1))
5860    CONTINUE
5900  CONTINUE
c
      RPOT = RPOT + 2D0*RPOTZ
      FIELDTEMP(1) = DREAL( ZS2 - ZS1 )
      FIELDTEMP(2) = -DIMAG( ZS2 + ZS1 )
      FIELDTEMP(3) = FIELDTEMP(3)+2D0*DREAL(ZS3)
      FIELD(1) = FIELD(1)+FIELDTEMP(1)*SCALE
      FIELD(2) = FIELD(2)+FIELDTEMP(2)*SCALE
      FIELD(3) = FIELD(3)-FIELDTEMP(3)*SCALE
c
      RETURN
      END
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
      SUBROUTINE BRTAEVPOT(LOCAL,X0Y0Z0,POINT,NTERMS,RPOT,SCALE,
     1             P,C,CS,CSINV)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
c
c  purpose:
c
c    evaluates local expansion at arbitrary point
c    this version only evaluates the potential, not the force.
c
c  on input:
c
c    scale: a scaling factor to scale the current box to
c      size 1.
c    local: coefficients of local expansion(scaled)
c    x0y0z0: the center of the expansion
c    point: point of evaluation
c    nterms: order of expansion
c    c,cs,csinv: precomputed coefficients.
c
c  on output:
c    rpot: computed potential
c
c  working array:
c
c    p: work arrays to hold legendre polynomials
c       and associated legendre functions, respectively.
c
c  subroutine called :
c
c  called from : ladapfmm()
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
c
      IMPLICIT NONE
c
      INTEGER *4 NTERMS
      REAL *8 SCALE
      REAL *8 X0Y0Z0(3),P(0:NTERMS,0:NTERMS),POINT(3),RPOT
      REAL *8 RLOC,RPOTZ
      REAL *8 C(0:60,0:60),CS(0:60,0:60),CSINV(0:60,0:60)
      COMPLEX *16 LOCAL(0:NTERMS,0:NTERMS)
c
c-----local variables.
c
      INTEGER *4 I,L,M,NTP1,NCALC
      REAL *8 RX,RY,RZ,PROJ,RR,CTHETA
      REAL *8 PRECIS,D,DD,CP,POWERS(0:60)
      COMPLEX *16 IMAG,EPHI(1:60),CPZ
c
      DATA IMAG/(0.0D0,1.0D0)/
      DATA PRECIS/1.0D-14/
c
c-----functions called.
c
      REAL *8 DREAL,DSQRT,DBLE
c
      RPOTZ = 0.0D0
c
c-----compute relevant functions of spherical coordinates
c     d = distance, ctheta = cos(theta), ephi = exp(i*phi)
c
      RX = POINT(1) - X0Y0Z0(1)
      RY = POINT(2) - X0Y0Z0(2)
      RZ = POINT(3) - X0Y0Z0(3)
      PROJ = RX*RX+RY*RY
      RR = PROJ+RZ*RZ
      PROJ = DSQRT(PROJ)
      D = DSQRT(RR)
      IF (D .LE. PRECIS) THEN
        CTHETA = 0.0D0
      ELSE
        CTHETA = RZ/D
      ENDIF
      IF ( PROJ .LE. PRECIS*D ) THEN
        EPHI(1) = 1.0D0
      ELSE
        EPHI(1) = RX/PROJ + IMAG*RY/PROJ
      ENDIF
c
c-----create array of powers of r and e^(i*m*phi).
c
      D = D*SCALE
      DD = D
      POWERS(0) = 1.0D0
      DO 5600 L = 1,NTERMS+2
        POWERS(L) = DD
        DD = DD*D
        EPHI(L+1) = EPHI(L)*EPHI(1)
5600  CONTINUE
c
c-----compute values of legendre functions
c
      CALL LGNDR(NTERMS,CTHETA,P)
c
c-----add contribution from constant term
c
      RPOT = RPOT + DREAL(LOCAL(0,0))
c
c-----add contributions from legendre polynomials
c
      DO 5700 L=1,NTERMS
        RLOC=DREAL(LOCAL(L,0))
        CP=RLOC*POWERS(L)*P(L,0)
        RPOT=RPOT + CP
5700  CONTINUE
c
c-----add contributions from associated legendre functions.
c
      DO 5900 L = 1,NTERMS
        DO 5800 M=1,L
          CPZ = LOCAL(L,M)*EPHI(M)
          RPOTZ=RPOTZ+DREAL(CPZ)*POWERS(L)*C(L,M)*P(L,M)
5800    CONTINUE
5900  CONTINUE
c
      RPOT = RPOT + 2D0*RPOTZ
      RETURN
      END
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
      SUBROUTINE MPOLETOEXP(MPOLE,NTERMS,NLAMBS,NUMTETS,
     1                      NEXPTOT,MEXPUP,MEXPDOWN,RLSC)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
c
c     this subroutine converts a multipole expansion mpole into the
c     corresponding exponential moment function mexp for the
c     both the +z direction and the -z direction.
c
c     u(x,y,z) = \sum_{n=0}^{nterms} \sum_{m=-n,n}
c                mpole(n,m) y_n^m(cos theta) e^{i m \phi}/r^{n+1}
c
c              = (1/2pi) \int_0^\infty e^{-\lambda z}
c                \int_0^{2\pi} e^{i\lambda(xcos(alpha)+ysin(alpha))}
c                mexpup(lambda,alpha) dalpha dlambda
c
c     for +z direction and
c
c              = (1/2pi) \int_0^\infty e^{\lambda z}
c                \int_0^{2\pi} e^{-i\lambda(xcos(alpha)+ysin(alpha))}
c                mexpdown(lambda,alpha) dalpha dlambda
c
c     for -z direction.
c
c     note: the expression for the -z direction corresponds to the
c     mapping (x,y,z) -> (-x,-y,-z), i.e. reflection through the origin.
c     one could also use rotation about the y axis, for which
c     (x,y,z) -> (-x,y,-z) but we stick to the reflected convention.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
c     note: the multipole expansion is assumed to have been rescaled
c           so that the box containing sources has unit dimension.
c
c     note: we only store mpole(n,m) for n,m >= 0, since mpole(n,-m)=
c           dconjg(mpole(n,m)). since we store the exponential
c           moment function in the fourier domain (w.r.t. the alpha
c           variable), we compute
c
c       m_lambda(m) = (i)**m \sum_{n=m}^n c(n,m) mpole(n,m) lambda^n
c
c           for m >= 0 only, where c(n,m) = 1/sqrt((n+m)!(n-m)!).
c
c       for possible future reference, it should be noted that
c       it is not true that m_lamb(-m) = dconjg(m_lamb(m)).
c       inspection of the integral formula for y_n^{-m} shows that
c       m_lamb(-m) = dconjg(m_lamb(m)) * (-1)**m.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
c
c  on input:
c
c     mpole(0:nterms,0:nterms): the multipole expansion
c     rlams(nlambs):  discretization points in lambda integral
c     nlambs:         number of discretization pts.
c     numtets(nlambs): number of fourier modes needed in expansion
c                    of alpha variable for each lambda value.
c                    note : the numtets is given by numthehalf().
c     nexptot =      sum_j numtets(j)
c     rlsc() : p_n^m for different lambda_k
c
c  on output:
c
c     mexpf(nexptot): fourier coefficients of the function
c                     mexp(lambda,alpha) for successive discrete
c                     lambda values. they are ordered as follows:
c                 mexpf(1,...,numtets(1)) = fourier modes
c                             for lambda_1
c                 mexpf(numtets(1)+1,...,numtets(2)) = fourier modes
c                             for lambda_2
c                 etc.
c
c     note by huangjf : in return, we will output
c       in mexpup, sum_{n=m}^{nterms} m_n^m*p_n^m. (all are scaled)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
c
c     note by jingfang :
c     1.this subroutine will compute the inner sum.
c       instead of compute all the nterms modes, only
c       the necessary modes are calculated. the number of mode
c       needed are provided by the subroutine numthehalf().
c       note that the number is always less than nterms needed,
c       and the worst case is nterms+1?.
c
c     2.
c       subroutine called :
c       called from : mkudexp(), mknsexp(), mkewexp()
c
c     3. the down-list will have the same fourier modes as the
c        up-list if we only change the sign of the z. so we don't need to
c        compute them separately.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
c
      IMPLICIT NONE
c
      INTEGER *4 NTERMS,NLAMBS,NUMTETS(NLAMBS),NEXPTOT
      REAL *8 RLSC(NLAMBS,0:NTERMS,0:NTERMS)
      COMPLEX *16 MPOLE(0:NTERMS,0:NTERMS)
      COMPLEX *16 MEXPUP(NEXPTOT), MEXPDOWN(NEXPTOT)
c
c-----local variables.
c
      INTEGER *4 NTOT,NL,MTH,NCURRENT,NM
      REAL *8 SGN
      COMPLEX *16 ZEYEP,ZTMP1,ZTMP2
c
c-----loop over multipole order to generate mexpup and mexpdown values.
c
      NTOT = 0
      DO NL = 1,NLAMBS
        SGN = -1.0D0
        ZEYEP = 1.0D0
        DO MTH = 0,NUMTETS(NL)-1
          NCURRENT = NTOT+MTH+1
          ZTMP1 = 0.0D0
          ZTMP2 = 0.0D0
          SGN = -SGN
          DO NM = MTH,NTERMS,2
            ZTMP1 = ZTMP1 +
     1        RLSC(NL,NM,MTH)*MPOLE(NM,MTH)
          ENDDO
          DO NM = MTH+1,NTERMS,2
            ZTMP2 = ZTMP2 +
     1        RLSC(NL,NM,MTH)*MPOLE(NM,MTH)
          ENDDO
          MEXPUP(NCURRENT) = (ZTMP1 + ZTMP2)*ZEYEP
          MEXPDOWN(NCURRENT) = SGN*(ZTMP1 - ZTMP2)*ZEYEP
          ZEYEP = ZEYEP*DCMPLX(0.0D0,1.0D0)
        ENDDO
        NTOT = NTOT+NUMTETS(NL)
      ENDDO
c
      RETURN
      END
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
      SUBROUTINE RLSCINI(RLSC,NLAMBS,RLAMS,NTERMS)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
c
c  purpose:
c
c    precomputes the coefficients for mpoletoexp for different levels.
c
c  on input :
c    scal : the scaling factor for the p_n^m. otherwise p_n^m will
c           overflow.
c    beta : the scaled frequency/beta for the current level.
c    nlambs : the total number of lamdas for the first integral.
c    rlams : the nodes/lamdas.
c    nterms : the total number of terms in the multipole expansion.
c
c  on output :
c    rlsc(n,m,nlamda) : the required information for the whole level.
c
c  subroutine called :
c
c  called from : ladapfmm()
c
c  note : the rlsc() will be used in the subroutine mpoletoexp().
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
c
      IMPLICIT NONE
c
      INTEGER *4 NLAMBS,NTERMS
      REAL *8 RLSC(NLAMBS,0:NTERMS,0:NTERMS)
      REAL *8 RLAMS(NLAMBS)
c
c-----local variables.
c
      INTEGER *4 I,J,K,NL
      REAL *8 RLAMPOW(0:100)
      REAL *8 FACTS(0:100)
      REAL *8 RMUL
c
      FACTS(0) = 1.0D0
      DO I = 1,100
        FACTS(I) = FACTS(I-1)*DSQRT(I+0.0D0)
      ENDDO
c
      DO NL = 1,NLAMBS
c
c-------compute powers of lambda_nl
c
        RLAMPOW(0) = 1.0D0
        RMUL = RLAMS(NL)
        DO J = 1,NTERMS
          RLAMPOW(J) = RLAMPOW(J-1)*RMUL
        ENDDO
        DO J = 0,NTERMS
          DO K = 0,J
            RLSC(NL,J,K) = RLAMPOW(J)/(FACTS(J-K)*FACTS(J+K))
          ENDDO
        ENDDO
      ENDDO
c
      RETURN
      END
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
      SUBROUTINE FTOPHYS(MEXPF,NLAMBS,RLAMS,NUMFOUR,NUMPHYS,
     1  NTHMAX,MEXPPHYS,FEXPE,FEXPO)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
c
c  purpose:
c
c    this subroutine evaluates the fourier expansion of the
c    exponential moment function m(\lambda,\alpha) at equispaced
c    nodes.
c
c  on input:
c
c    mexpf(*):     fourier coefficients of the function
c                  mexp(lambda,alpha) for discrete lambda values.
c                  they are ordered as follows:
c
c               mexpf(1,...,numfour(1)) = fourier modes for lambda_1
c               mexpf(numfour(1)+1,...,numfour(2)) = fourier modes
c                                              for lambda_2
c               etc.
c
c    nlambs:        number of discretization pts. in lambda integral
c    rlams(nlambs): discretization points in lambda integral.
c    numfour(j):   number of fourier modes in the expansion
c                      of the function m(\lambda_j,\alpha)
c    nthmax =      max_j numfour(j)
c    numphys : number of fourier modes in the plane wave expansion.
c    fexpe =      precomputed array of exponentials needed for
c                 fourier series evaluation. even terms.
c    fexpo =      precomputed array of exponentials needed for
c                 fourier series evaluation. odd terms.
c
c  note : we will keep these two terms because in
c         the helmholtz equation, we will need these.
c         however, in yukawa, it is not necessary to have
c         them separated.--huangjf
c
c  on output:
c
c    mexpphys(*):  discrete values of the moment function
c                  m(\lambda,\alpha), ordered as follows.
c
c        mexpphys(1),...,mexpphys(numphys(1)) = m(\lambda_1,0),...,
c             m(\lambda_1, 2*pi*(numphys(1)-1)/numphys(1)).
c        mexpphys(numphys(1)+1),...,mexpphys(numphys(2)) =
c             m(\lambda_2,0),...,
c                 m(\lambda_2, 2*pi*(numphys(2)-1)/numphys(2)).
c        etc.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
c
c      this subroutine computes the outer sum, it is possible
c      to do this using fft. but the current version will not do that.
c
c  subroutine called :
c
c  called from : mkudexp(), mknsexp(), mkewexp()
c
c  note :
c
c    the current subroutine computes sum_{m=-numfour, numfour}
c      e^{im*alpha} * i^|m| *inner(m)
c
c    the constant will left to the pw_local.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
c
      IMPLICIT NONE
c
      INTEGER *4 NLAMBS,NUMFOUR(NLAMBS),NUMPHYS(NLAMBS),NTHMAX
      REAL *8 RLAMS(NLAMBS)
c
      COMPLEX *16 MEXPF(1)
      COMPLEX *16 MEXPPHYS(1)
      COMPLEX *16 FEXPE(1)
      COMPLEX *16 FEXPO(1)
c
c-----local variables.
c
      INTEGER *4 I,IVAL,MM
      INTEGER *4 NFTOT,NPTOT,NEXTE,NEXTO
      REAL *8 SGN,PI,RT1,RT2,RTMP
      REAL *8 ALPHAS(0:200)
      COMPLEX *16 IMA
c
c-----functions
c
      REAL *8 DATAN
c
      DATA IMA/(0.0D0,1.0D0)/
c
      PI=DATAN(1.0D0)*4.0D0
c
      NFTOT = 0
      NPTOT  = 0
      NEXTE = 1
      NEXTO = 1
      DO 2000 I=1,NLAMBS
        DO 1200 IVAL=1,NUMPHYS(I)/2
          MEXPPHYS(NPTOT+IVAL) = MEXPF(NFTOT+1)
          DO MM = 2,NUMFOUR(I),2
            RT1 = DIMAG(FEXPE(NEXTE))*DREAL(MEXPF(NFTOT+MM))
            RT2 = DREAL(FEXPE(NEXTE))*DIMAG(MEXPF(NFTOT+MM))
            RTMP = 2*(RT1+RT2)
            NEXTE = NEXTE + 1
            MEXPPHYS(NPTOT+IVAL) = MEXPPHYS(NPTOT+IVAL) +
     1        DCMPLX(0.0D0,RTMP)
          ENDDO
c
          DO MM = 3,NUMFOUR(I),2
            RT1 = DREAL(FEXPO(NEXTO))*DREAL(MEXPF(NFTOT+MM))
            RT2 = DIMAG(FEXPO(NEXTO))*DIMAG(MEXPF(NFTOT+MM))
	    RTMP = 2*(RT1-RT2)
            NEXTO = NEXTO + 1
            MEXPPHYS(NPTOT+IVAL) = MEXPPHYS(NPTOT+IVAL) +
     1        DCMPLX(RTMP,0.0D0)
          ENDDO
1200    CONTINUE
c
        NFTOT = NFTOT+NUMFOUR(I)
        NPTOT = NPTOT+NUMPHYS(I)/2
2000  CONTINUE
c
      RETURN
      END
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
      SUBROUTINE PHYSTOF(MEXPF,NLAMBS,RLAMS,NUMFOUR,NUMPHYS,
     1  NTHMAX,MEXPPHYS,FEXPBACK)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
c
c  purpose:
c
c    this subroutine converts the discretized exponential moment function
c    into its fourier expansion.
c    it calculates the inner sum of the exp->local expansion.
c      (/sum_{j=1}^{m(k)} w(k,j)*e^{-im*alpha_j})/m(k) for k=1, nlambs,
c      and m=0, numfour.
c      numfour is the total number of the fourier modes.
c      or in other words, those l_n^m <>0.
c      the summation is over the numphys.
c
c  on input:
c
c    mexpphys(*):  discrete values of the moment function
c                  m(\lambda,\alpha), ordered as follows.
c
c        mexpphys(1),...,mexpphys(numphys(1)) = m(\lambda_1,0),...,
c             m(\lambda_1, 2*pi*(numphys(1)-1)/numphys(1)).
c        mexpphys(numphys(1)+1),...,mexpphys(numphys(2)) =
c             m(\lambda_2,0),...,
c                 m(\lambda_2, 2*pi*(numphys(2)-1)/numphys(2)).
c        etc.
c
c    nlambs:        number of discretization pts. in lambda integral
c    rlams(nlambs): discretization points in lambda integral.
c    numfour(j):   number of fourier modes in the expansion
c                      of the function m(\lambda_j,\alpha)
c    nthmax =      max_j numfour(j)
c    fexpback : contains the precomputed e^{-im *alpha_j}
c
c  on output:
c
c    mexpf(*):     fourier coefficients of the function
c                  mexp(lambda,m) for discrete lambda values.
c                  they are ordered as follows:
c
c               mexpf(1,...,numfour(1)) = fourier modes for lambda_1
c               mexpf(numfour(1)+1,...,numfour(2)) = fourier modes
c                                              for lambda_2
c               etc.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
c
      IMPLICIT NONE
c
      INTEGER *4 NLAMBS,NUMFOUR(NLAMBS),NUMPHYS(NLAMBS),NTHMAX
      REAL *8 RLAMS(NLAMBS)
c
      COMPLEX *16 MEXPF(1)
      COMPLEX *16 MEXPPHYS(1)
      COMPLEX *16 FEXPBACK(1)
c
c-----local variables.
c
      INTEGER *4 I,J,IVAL,MM
      INTEGER *4 NFTOT,NPTOT,NEXT,NALPHA,NALPHA2
      REAL *8 RTMP,PI
      COMPLEX *16 IMA,ZTMP
c
c-----functions
c
      REAL *8 DATAN
c
      DATA IMA/(0.0D0,1.0D0)/
c
      PI=DATAN(1.0D0)*4.0D0
      NFTOT = 0
      NPTOT  = 0
      NEXT  = 1
c
      DO 2000 I=1,NLAMBS
        NALPHA = NUMPHYS(I)
        NALPHA2 = NALPHA/2
        MEXPF(NFTOT+1) = 0.0D0
        DO IVAL=1,NALPHA2
          MEXPF(NFTOT+1) = MEXPF(NFTOT+1) +
     1      2D0*DREAL(MEXPPHYS(NPTOT+IVAL))
        ENDDO
        MEXPF(NFTOT+1) = MEXPF(NFTOT+1)/NALPHA
        DO MM = 3,NUMFOUR(I),2
          MEXPF(NFTOT+MM) = 0.0D0
          DO IVAL=1,NALPHA2
            RTMP = 2D0*DREAL(MEXPPHYS(NPTOT+IVAL))
            MEXPF(NFTOT+MM) = MEXPF(NFTOT+MM) +
     1        FEXPBACK(NEXT)*RTMP
            NEXT = NEXT+1
          ENDDO
          MEXPF(NFTOT+MM) = MEXPF(NFTOT+MM)/NALPHA
        ENDDO
        DO MM = 2,NUMFOUR(I),2
          MEXPF(NFTOT+MM) = 0.0D0
          DO IVAL=1,NALPHA2
            ZTMP = 2D0*DCMPLX(0.0D0,DIMAG(MEXPPHYS(NPTOT+IVAL)))
            MEXPF(NFTOT+MM) = MEXPF(NFTOT+MM) +
     1        FEXPBACK(NEXT)*ZTMP
            NEXT = NEXT+1
          ENDDO
          MEXPF(NFTOT+MM) = MEXPF(NFTOT+MM)/NALPHA
        ENDDO
c
        NFTOT = NFTOT+NUMFOUR(I)
        NPTOT = NPTOT+NUMPHYS(I)/2
2000  CONTINUE
c
      RETURN
      END
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
      SUBROUTINE MPSHIFT(IFL,MPOLE,MPOLEN,MARRAY,NTERMS,DC,RD,SC1,SC2)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
c
c  purpose:
c
c    this subroutine shifts the center of a child box multipole
c      expansion to the parent center, via the rotation scheme.
c      we rotate the coordinate system, shift along the z-axis,
c      and then rotate back.
c      there are eight possible child locations, defined in the
c      calling sequence of this routine by the parameters ifl
c      and rd (see below).
c
c  on input:
c
c     integer *4  ifl    = flag which describes the quadrant in which
c                          the child box lies (1,2,3 or 4).
c     complex *16 mpole  = coefficients of original multipole exp.
c     integer *4  nterms = integer indicates the terms retained in the
c                          expansion.
c     real *8     dc     = precomputed array containing
c                          the shifting coefficients
c                          along the z-axis. this is precomputed by
c                          the subroutine mpshftcoef() at the beginning
c                          of different levels.
c     real *8     rd     = precomputed array containing rotation matrix
c                          about y-axis.
c                 there are two possible y rotations, depending on
c                 whether the child box lies in +z half space or the
c                 -z half space. they are referred to in the calling
c                 program as rdp and rdm, respectively.
c                 this is precomputed in the subroutine rotgen().
c     real *8     sc1,sc2 = scaling parameters for mpole and mpolen,
c                 respectively.
c
c     complex *16 marray = work array
c
c  on output:
c
c     complex *16 mpolen = coefficients of shifted multipole exp.
c
c     note 1 : the rotation part is the same as the old subroutine
c              of the laplace equation. the shifting along the z-axis
c              is changed to the new version. the rotation matrix
c              is precomputed at the very beginning and the shifting
c              matrix is computed at different levels.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
c
      IMPLICIT NONE
c
      INTEGER *4 IFL,NTERMS
c
      REAL *8 SC1,SC2
      REAL *8 DC(0:2*NTERMS,0:2*NTERMS)
      REAL *8 RD(0:NTERMS,0:NTERMS,-NTERMS:NTERMS)
c
      COMPLEX *16 MPOLE(0:NTERMS,0:NTERMS)
      COMPLEX *16 MPOLEN(0:NTERMS,0:NTERMS)
      COMPLEX *16 MARRAY(0:NTERMS,0:NTERMS)
c
c-----local variables.
c
      INTEGER *4 L,M,JNEW,KNEW,MP
      REAL *8 ARG,DD,POWERS(0:60)
      COMPLEX *16 EPHI(0:60),IMAG
      DATA IMAG/(0.0D0,1.0D0)/
c
c-----functions called.
c
      COMPLEX *16 DCONJG
      COMPLEX *16 DCMPLX
c
      EPHI(0)=1.0D0
      ARG = DSQRT(2.0D0)/2.0D0
      IF (IFL.EQ.1) THEN
        EPHI(1) = DCMPLX(-ARG,ARG)
      ELSE IF (IFL.EQ.2) THEN
        EPHI(1) = DCMPLX(ARG,ARG)
      ELSE IF (IFL.EQ.3) THEN
        EPHI(1) = DCMPLX(ARG,-ARG)
      ELSE IF (IFL.EQ.4) THEN
        EPHI(1) = DCMPLX(-ARG,-ARG)
      ELSE
        CALL PRINF('error in lolo with ifl = *',IFL,1)
      ENDIF
c
c-----create array of powers of r and e^(i*m*phi).
c
      DD = -DSQRT(3.0D0)/2.0D0
      POWERS(0) = 1.0D0
      DO L = 1,NTERMS+1
        POWERS(L) = POWERS(L-1)*DD
        EPHI(L+1) = EPHI(L)*EPHI(1)
      ENDDO
c
c-----a rotation of phi radians about the z-axis in the
c     original coordinate system.
c
      DO L=0,NTERMS
        DO M=0,L
          MPOLEN(L,M)=DCONJG(EPHI(M))*MPOLE(L,M)
        ENDDO
      ENDDO
c
c-----a rotation about the y'-axis  in the rotated system.
c
      DO L=0,NTERMS
        DO M=0,L
          MARRAY(L,M)=MPOLEN(L,0)*RD(L,0,M)
          DO MP=1,L
            MARRAY(L,M)=MARRAY(L,M)+MPOLEN(L,MP)*RD(L,MP,M)+
     1        DCONJG(MPOLEN(L,MP))*RD(L,MP,-M)
          ENDDO
        ENDDO
      ENDDO
c
c-----shift along z-axis.
c
      DO JNEW=0,NTERMS
        DO KNEW=0,JNEW
          MPOLEN(JNEW,KNEW) = MARRAY(JNEW,KNEW)
          DO L=1,JNEW-IABS(KNEW)
            MPOLEN(JNEW,KNEW)=MPOLEN(JNEW,KNEW)+MARRAY(JNEW-L,KNEW)*
     1        POWERS(L)*DC(JNEW-KNEW,L)*DC(JNEW+KNEW,L)
          ENDDO
        ENDDO
      ENDDO
c
c-----reverse rotation about the y'-axis.
c
      DO L=0,NTERMS
        DO M=0,L,2
          MARRAY(L,M)=MPOLEN(L,0)*RD(L,0,M)
          DO MP=1,L,2
            MARRAY(L,M)=MARRAY(L,M)-(MPOLEN(L,MP)*RD(L,MP,M)+
     1        DCONJG(MPOLEN(L,MP))*RD(L,MP,-M))
          ENDDO
          DO MP=2,L,2
            MARRAY(L,M)=MARRAY(L,M)+(MPOLEN(L,MP)*RD(L,MP,M)+
     1        DCONJG(MPOLEN(L,MP))*RD(L,MP,-M))
          ENDDO
        ENDDO
        DO M=1,L,2
          MARRAY(L,M)=-MPOLEN(L,0)*RD(L,0,M)
          DO MP=1,L,2
            MARRAY(L,M)=MARRAY(L,M)+(MPOLEN(L,MP)*RD(L,MP,M)+
     1        DCONJG(MPOLEN(L,MP))*RD(L,MP,-M))
          ENDDO
          DO MP=2,L,2
            MARRAY(L,M)=MARRAY(L,M)-(MPOLEN(L,MP)*RD(L,MP,M)+
     1        DCONJG(MPOLEN(L,MP))*RD(L,MP,-M))
          ENDDO
        ENDDO
      ENDDO
c
c-----rotate back phi radians about the z-axis in the above system.
c
      DO L=0,NTERMS
        DO M=0,L
          MPOLEN(L,M)=EPHI(M)*MARRAY(L,M)
        ENDDO
      ENDDO
c
c-----rescale expansion from sc1 to sc2
c
      POWERS(0) = 1.0D0
      DD = SC2/SC1
      DO L = 1,NTERMS+1
        POWERS(L) = POWERS(L-1)*DD
      ENDDO
      DO L=0,NTERMS
        DO M=0,L
          MPOLEN(L,M)=MPOLEN(L,M)*POWERS(L)
        ENDDO
      ENDDO
c
      RETURN
      END
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
      SUBROUTINE LCSHIFT(IFL,LOCAL,LOCALN,MARRAY,NTERMS,DC,RD,SC1,SC2)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
c
c  purpose:
c
c     this subroutine shifts the local expansion of a parent cell
c     to the center of one of its children, via the rotation scheme.
c     that is, we rotate the coordinate system, shift along the z-axis,
c     and then rotate back.
c     there are eight possible child locations, defined in the
c     calling sequence of this routine by the parameters ifl
c     and rd (see below).
c
c on input:
c
c     integer *4  ifl    = flag which describes the quadrant in which
c                          the child box lies (1,2,3 or 4).
c     complex *16 local  = coefficients of original multipole exp.
c     integer *4  nterms = integer indicates the terms retained in the
c                          expansion.
c     real *8     dc     = precomputed array containing coefficients
c                          for the local translation along the z axis.
c                          this is precomputed by the subroutine
c                          lcshftcoef()
c     real *8     rd     = precomputed array containing rotation matrix
c                          about y-axis.
c                 there are two possible y rotations, depending on
c                 whether the child box lies in +z half space or the
c                 -z half space. they are referred to in the calling
c                 program as rdp and rdm, respectively.
c                          this is precomputed by the subroutine
c                          rotgen().
c     real *8     sc1,sc2 = scaling parameters for local and localn,
c                 respectively.
c     complex *16 marray = work array
c
c on output:
c
c     complex *16 localn = coefficients of shifted multipole exp.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
c
      IMPLICIT NONE
c
      INTEGER *4 IFL,NTERMS
      REAL *8 DC(0:2*NTERMS,0:2*NTERMS)
      REAL *8 RD(0:NTERMS,0:NTERMS,-NTERMS:NTERMS)
c
      COMPLEX *16 LOCAL(0:NTERMS,0:NTERMS)
      COMPLEX *16 LOCALN(0:NTERMS,0:NTERMS)
      COMPLEX *16 MARRAY(0:NTERMS,0:NTERMS)
c
c-----local variables.
c
      INTEGER *4 L,M,MP,JNEW,KNEW,LL
      REAL *8 ARG
      REAL *8 DD,POWERS(0:55),SC1,SC2
      COMPLEX *16 EPHI(0:60),IMAG
c
c-----functions.
c
      REAL *8 DSQRT
      COMPLEX *16 DCOMPLX,DCONJG
c
      DATA IMAG/(0.0D0,1.0D0)/
c
      EPHI(0)=1.0D0
      ARG = DSQRT(2.0D0)/2.0D0
      IF (IFL.EQ.1) THEN
        EPHI(1) = DCMPLX(ARG,-ARG)
      ELSE IF (IFL.EQ.2) THEN
        EPHI(1) = DCMPLX(-ARG,-ARG)
      ELSE IF (IFL.EQ.3) THEN
        EPHI(1) = DCMPLX(-ARG,ARG)
      ELSE IF (IFL.EQ.4) THEN
        EPHI(1) = DCMPLX(ARG,ARG)
      ELSE
        CALL PRINF('error in lolo with ifl = *',IFL,1)
      ENDIF
c
c-----create array of powers of r and e^(i*m*phi).
c
      DD = -DSQRT(3.0D0)/4.0D0
      POWERS(0) = 1.0D0
      DO L = 1,NTERMS+1
        POWERS(L) = POWERS(L-1)*DD
        EPHI(L+1) = EPHI(L)*EPHI(1)
      ENDDO
c
c-----a rotation of phi radians about the z-axis in the
c     original coordinate system.
c
      DO L=0,NTERMS
        DO M=0,L
          LOCALN(L,M)=DCONJG(EPHI(M))*LOCAL(L,M)
        ENDDO
      ENDDO
c
c-----a rotation about the y'-axis to align z' axis.
c
      DO L=0,NTERMS
        DO M=0,L
          MARRAY(L,M)=LOCALN(L,0)*RD(L,0,M)
          DO MP=1,L
            MARRAY(L,M)=MARRAY(L,M)+LOCALN(L,MP)*RD(L,MP,M)+
     1        DCONJG(LOCALN(L,MP))*RD(L,MP,-M)
          ENDDO
        ENDDO
      ENDDO
c
c-----shift along z'-axis.
c
      DO JNEW= 0,NTERMS
        DO KNEW=0,JNEW
          LOCALN(JNEW,KNEW) = MARRAY(JNEW,KNEW)
          DO L=1,NTERMS-JNEW
            LL=L+JNEW
            LOCALN(JNEW,KNEW)=LOCALN(JNEW,KNEW)+MARRAY(LL,KNEW)*
     1        POWERS(L)*DC(LL+KNEW,L)*DC(LL-KNEW,L)
          ENDDO
        ENDDO
      ENDDO
c
c-----rotate back about the y'-axis.
c
      DO L=0,NTERMS
        DO M=0,L,2
          MARRAY(L,M)=LOCALN(L,0)*RD(L,0,M)
          DO MP=1,L,2
            MARRAY(L,M)=MARRAY(L,M)-(LOCALN(L,MP)*RD(L,MP,M)+
     1        DCONJG(LOCALN(L,MP))*RD(L,MP,-M))
          ENDDO
          DO MP=2,L,2
            MARRAY(L,M)=MARRAY(L,M)+(LOCALN(L,MP)*RD(L,MP,M)+
     1        DCONJG(LOCALN(L,MP))*RD(L,MP,-M))
          ENDDO
        ENDDO
        DO M=1,L,2
          MARRAY(L,M)=-LOCALN(L,0)*RD(L,0,M)
          DO MP=1,L,2
            MARRAY(L,M)=MARRAY(L,M)+(LOCALN(L,MP)*RD(L,MP,M)+
     1        DCONJG(LOCALN(L,MP))*RD(L,MP,-M))
          ENDDO
          DO MP=2,L,2
            MARRAY(L,M)=MARRAY(L,M)-(LOCALN(L,MP)*RD(L,MP,M)+
     1        DCONJG(LOCALN(L,MP))*RD(L,MP,-M))
          ENDDO
        ENDDO
      ENDDO
c
c-----rotate back about the z-axis.
c
      DO L=0,NTERMS
        DO M=0,L
          LOCALN(L,M)=EPHI(M)*MARRAY(L,M)
        ENDDO
      ENDDO
c
c-----rescale expansion from sc1 to sc2
c
      POWERS(0) = 1.0D0
      DD = SC1/SC2
      DO L = 1,NTERMS+1
        POWERS(L) = POWERS(L-1)*DD
      ENDDO
      DO L=0,NTERMS
        DO M=0,L
          LOCALN(L,M)=LOCALN(L,M)*POWERS(L)
        ENDDO
      ENDDO
c
      RETURN
      END
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
c
c     this file contains all of the expansion creation routines for
c     a parent box from its four children.
c
c     mkexps creates the table of shifting coefficients in the physical
c     domain ( e^{lambda z}, e^{ i x cos u}, e^{ i y sin  u} ) for
c     a given lambda discretization.
c
c     mkfexp creates the table of exponentials needed for mapping from
c     fourier to physical domain
c
c     mkudexp creates all up and down expansions centered at child 1
c
c     mknsexp creates all north and south expansions centered at child 1
c
c     mkewexp creates all east and west expansions centered at child 1
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
      SUBROUTINE MKEXPS(RLAMS,NLAMBS,NUMPHYS,NEXPTOTP,XS,YS,ZS)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
c
c  purpose:
c
c     this subroutine computes the tables of exponentials needed
c     for translating exponential representations of harmonic
c     functions, discretized via norman's quadratures.
c
c     u   = \int_0^\infty e^{-(lambda+beta) z}
c     \int_0^{2\pi} e^{i \dsqrt(lambda*lambda+2*beta*lambda)
c                         (x cos(u)+y sin(u))}
c           mexpphys(lambda,u) du dlambda
c
c     mexpphys(*):  discrete values of the moment function
c                   m(\lambda,u), ordered as follows.
c
c         mexpphys(1),...,mexpphys(numphys(1)) = m(\lambda_1,0),...,
c              m(\lambda_1, 2*pi*(numphys(1)-1)/numphys(1)).
c         mexpphys(numphys(1)+1),...,mexpphys(numphys(2)) =
c              m(\lambda_2,0),...,
c                  m(\lambda_2, 2*pi*(numphys(2)-1)/numphys(2)).
c         etc.
c         note : in the current version, only half of the modes are
c                stored because of the symmetry of u and u+pi.
c
c  on input:
c
c     beta : the scaled frequency.
c     rlams(nlambs)   discretization points in lambda integral
c     nlambs          number of discret. pts. in lambda integral
c     numphys(j)     number of nodes in u integral needed
c                    for corresponding lambda =  lambda_j.
c     nexptotp        sum_j numphys(j)
c
c  on output:
c
c        define w1=\lambda_j+beta, and w2= sqrt(lambda_j**2+2*beta*lambda_j)
c     xs(1,nexptotp)   e^{i w2 (cos(u_k)}  in above ordering
c     xs(2,nexptotp)   e^{i w2 (2 cos(u_k)}  in above ordering.
c     xs(3,nexptotp)   e^{i w2 (3 cos(u_k)}  in above ordering.
c     ys(1,nexptotp)   e^{i w2 (sin(u_k)}  in above ordering.
c     ys(2,nexptotp)   e^{i w2 (2 sin(u_k)}  in above ordering.
c     ys(3,nexptotp)   e^{i w2 (3 sin(u_k)}  in above ordering.
c     zs(1,nexptotp)   e^{-w1}     in above ordering.
c     zs(2,nexptotp)    e^{-2 w1}   in above ordering.
c     zs(3,nexptotp)    e^{-3 w1}   in above ordering.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
c
      IMPLICIT NONE
c
      INTEGER *4 NLAMBS,NUMPHYS(NLAMBS),NEXPTOTP
c
      REAL *8 RLAMS(NLAMBS)
      REAL *8 ZS(3,NEXPTOTP)
c
      COMPLEX *16 XS(3,NEXPTOTP)
      COMPLEX *16 YS(3,NEXPTOTP)
c
c-----local variables.
c
      INTEGER *4 NTOT,NL,MTH,NCURRENT
      REAL *8 PI,HU,U
      COMPLEX *16 IMA
c
c-----functions.
c
      REAL *8 DATAN,DSQRT,DEXP,DCOS,DSIN
      COMPLEX *16 CDEXP
      DATA IMA/(0.0D0,1.0D0)/
c
c-----loop over each lambda value
c
      PI = 4*DATAN(1.0D0)
      NTOT = 0
      DO NL = 1,NLAMBS
        HU=2D0*PI/NUMPHYS(NL)
        DO MTH = 1,NUMPHYS(NL)/2
          U = (MTH-1)*HU
          NCURRENT = NTOT+MTH
          ZS(1,NCURRENT)=DEXP( -RLAMS(NL) )
          ZS(2,NCURRENT)=ZS(1,NCURRENT)*ZS(1,NCURRENT)
          ZS(3,NCURRENT)=ZS(2,NCURRENT)*ZS(1,NCURRENT)
          XS(1,NCURRENT)=CDEXP(IMA*RLAMS(NL)*DCOS(U))
          XS(2,NCURRENT)=XS(1,NCURRENT)*XS(1,NCURRENT)
          XS(3,NCURRENT)=XS(2,NCURRENT)*XS(1,NCURRENT)
          YS(1,NCURRENT)=CDEXP(IMA*RLAMS(NL)*DSIN(U))
          YS(2,NCURRENT)=YS(1,NCURRENT)*YS(1,NCURRENT)
          YS(3,NCURRENT)=YS(2,NCURRENT)*YS(1,NCURRENT)
        ENDDO
        NTOT = NTOT+NUMPHYS(NL)/2
      ENDDO
c
      RETURN
      END
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
      SUBROUTINE MKFEXP(NLAMBS,NUMFOUR,NUMPHYS,FEXPE,FEXPO,FEXPBACK)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
c
c  purpose:
c
c     precomputes the e^(im*alpha) needed in mp->pw->local.
c
c     this subroutine computes the tables of exponentials needed
c     for mapping from fourier to physical domain.
c     in order to minimize storage, they are organized in a
c     one-dimenional array corresponding to the order in which they
c     are accessed by subroutine ftophys.
c
c     size of fexpe, fexpo =          40000   for nlambs = 39
c     size of fexpe, fexpo =          15000   for nlambs = 30
c     size of fexpe, fexpo =           4000   for nlambs = 20
c     size of fexpe, fexpo =            400   for nlambs = 10
c
c
c  on input :
c
c       nlambs : the total number of nodes for the outer integral.
c       numfour : contains the number of nodes for the fourier
c                 representation.
c       numphys : contains the number of nodes for the inner integral
c                 for the plane wave expansion.
c
c  on output :
c
c       fexpe : the exponentials for the fourier modes. e^(im*alpha)
c               where m is all the fourier modes and alpha comes
c               from the physical modes. odd terms.
c       fexpo : the even terms.
c       fexpback : the exponentials used for the translation
c                  from plane wave to local.
c
c     functions called :
c
c     called from :
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
c
      IMPLICIT NONE
c
      INTEGER *4  NLAMBS,NUMPHYS(NLAMBS),NUMFOUR(NLAMBS)
c
      COMPLEX *16 FEXPE(1)
      COMPLEX *16 FEXPO(1)
      COMPLEX *16 FEXPBACK(1)
c
c-----local variables.
c
      INTEGER *4 I,J,MM,NEXTE,NEXTO,NEXT,NALPHA,NALPHA2
      REAL *8 HALPHA,ALPHA,PI
      COMPLEX *16 IMA
      DATA IMA/(0.0D0,1.0D0)/
c
      PI = 4D0*DATAN(1.0D0)
      NEXTE = 1
      NEXTO = 1
      DO I=1,NLAMBS
        NALPHA = NUMPHYS(I)
        NALPHA2 = NALPHA/2
        HALPHA=2*PI/NALPHA
        DO J=1,NALPHA2
          ALPHA=(J-1)*HALPHA
          DO MM = 2,NUMFOUR(I),2
            FEXPE(NEXTE)  = CDEXP(IMA*(MM-1)*ALPHA)
            NEXTE = NEXTE + 1
          ENDDO
          DO MM = 3,NUMFOUR(I),2
            FEXPO(NEXTO)  = CDEXP(IMA*(MM-1)*ALPHA)
            NEXTO = NEXTO + 1
          ENDDO
        ENDDO
      ENDDO
c
      NEXT = 1
      DO I=1,NLAMBS
        NALPHA = NUMPHYS(I)
        NALPHA2 = NALPHA/2
        HALPHA=2*PI/NALPHA
        DO MM = 3,NUMFOUR(I),2
          DO J=1,NALPHA2
            ALPHA=(J-1)*HALPHA
            FEXPBACK(NEXT)  = CDEXP(-IMA*(MM-1)*ALPHA)
            NEXT = NEXT + 1
          ENDDO
        ENDDO
        DO MM = 2,NUMFOUR(I),2
          DO J=1,NALPHA2
            ALPHA=(J-1)*HALPHA
            FEXPBACK(NEXT)  = CDEXP(-IMA*(MM-1)*ALPHA)
            NEXT = NEXT + 1
          ENDDO
        ENDDO
      ENDDO
c
      RETURN
      END
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
      SUBROUTINE FORMLC(X0Y0Z0,ZPARTS,CHARGE,NPARTS,
     1  LOCAL,NTERMS,SCALE,P,C)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
c
c  purpose:
c    this subroutine forms the local expansion caused by the nparts
c      particles in the box.
c
c  on input:
c    x0y0z0: center of the expansion
c    nparts: number of sources
c    zparts(3,nparts): array of coordinates of sources
c    charge(nparts): array of strengths of sources
c    nparts: the total number of particles.
c    nterms: order of desired expansion
c    scale: the scaling factor.
c
c  on output:
c    local: coefficients of local expansion
c
c  working space :
c    p : used for storing the associate legendre polynomials.
c
c  subroutine called : dsqrt(), in(), lgndr()
c
c  called from : ladapfmm()
c
c  note 1: this subroutine needs the precomputed variables c(,)
c          derived from entry frmini()
c
c       2: the multipole expansion is scaled to avoid over- and
c          under-flow.
c
c       3: only the n=0, ... ,nterms, m=0, ..., nterms coefficients
c          are calculated.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
c
      IMPLICIT NONE
c
      INTEGER *4 NTERMS,NPARTS
      REAL *8 X0Y0Z0(3),ZPARTS(3,NPARTS),CHARGE(NPARTS)
      REAL *8 P(0:NTERMS,0:NTERMS)
      REAL *8 C(0:60,0:60),SCALE
      COMPLEX *16 LOCAL(0:NTERMS,0:NTERMS)
c
c-----local variables.
c
      INTEGER *4 I,L,M,NCALC
      REAL *8 PRECIS,D,CP,PROJ
      REAL *8 RX,RY,RZ,RR,CTHETA
      REAL *8 POWERS(0:60)
      COMPLEX *16 IMAG,EPHI(1:60)
c
      DATA IMAG/(0.0D0,1.0D0)/
      DATA PRECIS/1.0D-14/
c
      DO 4900 I = 1, NPARTS
        RX = ZPARTS(1,I) - X0Y0Z0(1)
        RY = ZPARTS(2,I) - X0Y0Z0(2)
        RZ = ZPARTS(3,I) - X0Y0Z0(3)
c
c-------compute  distance rr
c         ctheta = cos(theta)
c         ephi(1)=e^(i*phi)
c
        PROJ = RX*RX+RY*RY
        RR = PROJ+RZ*RZ
        PROJ = DSQRT(PROJ)
        D = DSQRT(RR)
c
c-------note: here is a hack. when computing cos(theta) as
c         rz/d, we have to be careful about the possibility
c         of d  being 0 (in which case ctheta is not well
c         defined - we arbitrarily set it to 1.)
c
        IF ( D .LE. PRECIS ) THEN
          CTHETA = 1.0D0
        ELSE
          CTHETA = RZ/D
        ENDIF
c
        IF ( PROJ .LE. PRECIS*D ) THEN
          EPHI(1) = 1.0D0
        ELSE
          EPHI(1) = RX/PROJ - IMAG*RY/PROJ
        ENDIF
c
c-------create array of powers of e^(-i*phi).
c
        D=1D0/D
        POWERS(0)=1.0D0
        POWERS(1)=D
        D=D/SCALE
        DO 4100 L = 2,NTERMS+2
          EPHI(L) = EPHI(L-1)*EPHI(1)
          POWERS(L) = POWERS(L-1)*D
4100    CONTINUE
c
        LOCAL(0,0) = LOCAL(0,0) + CHARGE(I)*POWERS(1)
c
c-------compute legendre polynomials of argument cos(theta) = ctheta
c         and add contributions from legendre polynomials
c
        CALL LGNDR(NTERMS,CTHETA,P)
c
        DO 4300 L = 1,NTERMS
          CP = CHARGE(I)*P(L,0)*POWERS(L+1)
          LOCAL(L,0) = LOCAL(L,0) + CP
4300    CONTINUE
c
c-------add contributions from associated legendre functions.
c
        DO 4500 L = 1,NTERMS
          DO 4400 M=1,L
            CP = CHARGE(I)*POWERS(L+1)*C(L,M)*P(L,M)
            LOCAL(L,M) = LOCAL(L,M) + CP*EPHI(M)
4400      CONTINUE
4500    CONTINUE
4900  CONTINUE
c
      RETURN
      END
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
      SUBROUTINE BRMPEV(MPOLE,X0Y0Z0,POINT,NTERMS,RPOT,FIELD,SCALE,P,
     1  C,CS,CSINV)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
c
c  purpose:
c    evaluates multipole expansion at arbitrary point.
c
c  on input:
c    mpole: coefficients of multipole expansion(scaled)
c    x0y0z0: the center of the expansion
c    point: point of evaluation
c    nterms: order of expansion
c    p: work arrays to hold legendre polynomials
c       and associated legendre functions, respectively.
c
c  on output:
c    rpot: computed potential
c    field: computed field.
c
c  note: the current version will compute the potential only,
c            and the computation of the field will be added later.
c
c************************************************************************
c
      IMPLICIT NONE
c
      INTEGER *4 NTERMS
      REAL *8 X0Y0Z0(3),POINT(3),RPOT
      REAL *8 FIELD(3)
      REAL *8 P(0:NTERMS,0:NTERMS)
      REAL *8 C(0:60,0:60),CS(0:60,0:60),CSINV(0:60,0:60)
      REAL *8 SCALE
      COMPLEX *16 MPOLE(0:NTERMS,0:NTERMS)
c
c-----local variables.
c
      INTEGER *4 NTP1,I,L,M,NCALC
      REAL *8 PRECIS,D,CP,PROJ
      REAL *8 POWERS(0:60)
      REAL *8 RX,RY,RZ,RR,CTHETA
      REAL *8 FIELD0(3),RMP,RTEMP,RPOTZ
      COMPLEX *16 ZS1,ZS2,ZS3
      COMPLEX *16 IMAG,EPHI(1:60),CPZ
      DATA IMAG/(0.0D0,1.0D0)/
      DATA PRECIS/1.0D-14/
c
c-----functions called.
c
      REAL *8 DREAL
c
      NTP1 = NTERMS+1
      RPOTZ = 0.0D0
      ZS1 = 0.0D0
      ZS2 = 0.0D0
      ZS3 = 0.0D0
c
c-----compute relevant functions of spherical coordinates
c       d = distance*beta, ctheta = cos(theta), ephi = exp(i*phi)
c
      RX = POINT(1) - X0Y0Z0(1)
      RY = POINT(2) - X0Y0Z0(2)
      RZ = POINT(3) - X0Y0Z0(3)
      PROJ = RX*RX+RY*RY
      RR = PROJ+RZ*RZ
      PROJ = DSQRT(PROJ)
      D = DSQRT(RR)
c
      IF (D .LE. PRECIS) THEN
        CTHETA = 0.0D0
      ELSE
        CTHETA = RZ/D
      ENDIF
      IF ( PROJ .LE. PRECIS*D ) THEN
        EPHI(1) = 1.0D0
      ELSE
        EPHI(1) = RX/PROJ + IMAG*RY/PROJ
      ENDIF
c
      D=1.0D0/D
      POWERS(0) = 1.0D0
      POWERS(1) = D
      D=D/SCALE
c
c-----create array of powers of e^(i*m*phi) and inverse powers of r.
c
      DO 5600 L = 1,NTERMS+2
        POWERS(L+1) = POWERS(L)*D
        EPHI(L+1) = EPHI(L)*EPHI(1)
5600  CONTINUE
c
c-----compute values of legendre functions
c
      CALL LGNDR(NTERMS,CTHETA,P)
c
c----- add contribution from monopole moment
c
      RMP=DREAL(MPOLE(0,0))
      RPOT=RPOT+RMP*POWERS(1)
      CPZ=EPHI(1)*(RMP*POWERS(2)*C(1,1)*P(1,1)*CSINV(1,1))
      ZS1=ZS1+CPZ
      CP=RMP*POWERS(2)*P(1,0)*CS(0,0)*CSINV(1,0)
      FIELD0(3)=CP
c
c-----add contributions from legendre polynomials
c
      DO 5300 L = 1,NTERMS
        RMP=DREAL(MPOLE(L,0))
        CP=RMP*POWERS(L+1)*P(L,0)
        RPOT=RPOT+CP
        ZS1=ZS1+EPHI(1)*
     1    (RMP*POWERS(L+2)*C(L+1,1)*P(L+1,1)*CS(L,0)*CSINV(L+1,1))
        CPZ = MPOLE(L,1)
        RTEMP=POWERS(L+2)*P(L+1,0)*CSINV(L+1,0)
        ZS2=ZS2+CPZ*(RTEMP*CS(L,1))
        CP=RMP*RTEMP*CS(L,0)
        FIELD0(3)=FIELD0(3)+CP
5300  CONTINUE
c
c-----add contributions from associated legendre functions.
c
      DO 5500 L=1,NTERMS
        DO 5400 M=1,L
          CPZ = MPOLE(L,M)*(POWERS(L+1)*C(L,M)*P(L,M))
          RPOTZ = RPOTZ + DREAL(CPZ*EPHI(M))
          CPZ=MPOLE(L,M)*(CS(L,M)*POWERS(L+2))
          ZS1 = ZS1 + CPZ*EPHI(M+1)*(
     1      CSINV(L+1,M+1)*C(L+1,M+1)*P(L+1,M+1))
          IF (M .GT. 1) THEN
            ZS2 =ZS2+ CPZ*EPHI(M-1)*(
     1      CSINV(L+1,M-1)*C(L+1,M-1)*P(L+1,M-1))
          ENDIF
          ZS3 = ZS3 + CPZ*EPHI(M)*
     1      (C(L+1,M)*P(L+1,M)*CSINV(L+1,M))
5400    CONTINUE
5500  CONTINUE
c
      RPOT = RPOT + 2D0*RPOTZ
      FIELD0(1) = DREAL( ZS2 - ZS1 )
      FIELD0(2) = -DIMAG( ZS2 + ZS1 )
      FIELD0(3) = FIELD0(3) + 2D0*DREAL(ZS3)
      FIELD(1) = FIELD(1)+FIELD0(1)*SCALE
      FIELD(2) = FIELD(2)+FIELD0(2)*SCALE
      FIELD(3) = FIELD(3)+FIELD0(3)*SCALE
c
      RETURN
      END
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
      SUBROUTINE DNFORMMP(X0Y0Z0,ZPARTS,CHARGE,DN,NPARTS,
     1  MPOLE,NTERMS,SCALE,P,C,CS,CSINV)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
c
c  purpose:
c
c    this subroutine forms the multipole expansion cause by the nparts
c      particles in the box.
c
c  on input:
c
c    x0y0z0: center of the expansion
c    nparts: number of sources
c    zparts(3,nparts): array of coordinates of sources
c    charge(nparts): array of strengths of sources
c    dn(3,nparts): the normal direction.
c    nparts: the total number of particles.
c    nterms: order of desired expansion
c    scale: the scaling factor.
c    c: precomputed numbers.
c
c  on output:
c
c    mpole: coefficients of multipole expansion
c
c  working space :
c
c    p: used for storing the associate legendre polynomials.
c
c  subroutine called : dsqrt(), lgndr()
c  called from : ladapfmm()
c
c  note 1: this subroutine needs the precomputed variables c(,)
c          derived from entry frmini()
c
c       2: the multipole expansion is scaled to avoid over- and
c          under-flow.
c
c       3: only the n=0, ... ,nterms, m=0, ..., nterms coefficients
c          are calculated.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
c
      IMPLICIT NONE
c
      INTEGER *4 NTERMS,NPARTS
      REAL *8 SCALE
      REAL *8 X0Y0Z0(3),ZPARTS(3,NPARTS),CHARGE(NPARTS)
      REAL *8 DN(3,NPARTS)
      REAL *8 P(0:NTERMS,0:NTERMS),C(0:60,0:60)
      REAL *8 CS(0:60,0:60),CSINV(0:60,0:60)
      COMPLEX *16 MPOLE(0:NTERMS,0:NTERMS)
c
c-----local variables.
c
      INTEGER *4 I,L,M,NCALC
      REAL *8 POWERS(0:60)
      REAL *8 PRECIS,D,CP,PROJ
      REAL *8 RX,RY,RZ,RR,CTHETA,SR2
      COMPLEX *16 IMAG,EPHI(1:60),CPZ
      COMPLEX *16 DNZ1,DNZ2
c
      DATA IMAG/(0.0D0,1.0D0)/
      DATA PRECIS/1.0D-14/
c
      DO 4900 I = 1, NPARTS
        RX = ZPARTS(1,I) - X0Y0Z0(1)
        RY = ZPARTS(2,I) - X0Y0Z0(2)
        RZ = ZPARTS(3,I) - X0Y0Z0(3)
c
c-------compute  distance rr
c         ctheta = cos(theta)
c         ephi(1)=e^(i*phi)
c
        PROJ = RX*RX+RY*RY
        RR = PROJ+RZ*RZ
        PROJ = DSQRT(PROJ)
        D = DSQRT(RR)
c
c-------note: here is a hack. when computing cos(theta) as
c             rz/d, we have to be careful about the possibility
c             of d  being 0 (in which case ctheta is not well
c             defined - we arbitrarily set it to 1.)
c
        IF ( D .LE. PRECIS ) THEN
          CTHETA = 1.0D0
        ELSE
          CTHETA = RZ/D
        ENDIF
        IF ( PROJ .LE. PRECIS*D ) THEN
          EPHI(1) = 1.0D0
        ELSE
          EPHI(1) = RX/PROJ - IMAG*RY/PROJ
        ENDIF
c
c-------create array of powers of r and powers of e^(i*phi).
c
        D = D*SCALE
        POWERS(0) = 1.0D0
c
        DO 4100 L = 1,NTERMS+1
          POWERS(L) = POWERS(L-1)*D
          EPHI(L+1) = EPHI(L)*EPHI(1)
4100    CONTINUE
c
c-------contribution from (0,0) term.
c
        CP=CHARGE(I)*SCALE
        DNZ1=DN(1,I)-DN(2,I)*IMAG
        DNZ2=DCONJG(DNZ1)
        MPOLE(1,1)=MPOLE(1,1)+DNZ1*(CP*CS(1,1))
        MPOLE(1,0)=MPOLE(1,0)-CP*DN(3,I)
c
c-------compute legendre polynomials of argument cos(theta) = ctheta
c         and add contributions from legendre polynomials
c
        CALL LGNDR(NTERMS,CTHETA,P)
c
c-------contributions from (l,0) terms.
c
        DO 4300 L = 1,NTERMS-1
          CP=CHARGE(I)*POWERS(L)*P(L,0)*SCALE
          CPZ=CP*CSINV(L+1,1)*CS(L,0)/2D0
          MPOLE(L+1,1)=MPOLE(L+1,1)+CPZ*DNZ1
          MPOLE(L+1,0)=MPOLE(L+1,0)-CP*DBLE(L+1)*DN(3,I)
4300    CONTINUE
c
c-------add contributions from associated legendre functions.
c
        DO 4500 L = 1,NTERMS-1
          DO 4400 M=1,L
            CP=CHARGE(I)*POWERS(L)*C(L,M)*P(L,M)*SCALE
            CPZ=CP*EPHI(M)
c
            SR2=CSINV(L+1,M-1)*CS(L,M)
            IF (M.EQ.1) THEN
              CP=SR2*DREAL(CPZ*DNZ2)
              MPOLE(L+1,0)=MPOLE(L+1,0)-CP
            ELSE
              MPOLE(L+1,M-1)=MPOLE(L+1,M-1)-
     1          SR2/2D0*CPZ*DNZ2
            ENDIF
c
            SR2=CSINV(L+1,M+1)*CS(L,M)/2D0
            MPOLE(L+1,M+1)=MPOLE(L+1,M+1)+SR2*CPZ*DNZ1
            MPOLE(L+1,M)=MPOLE(L+1,M)-(CSINV(L+1,M)*CS(L,M))*
     1        CPZ*DN(3,I)
4400      CONTINUE
4500    CONTINUE
4900  CONTINUE
c
      RETURN
      END
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
      SUBROUTINE DNLC(X0Y0Z0,ZPARTS,CHARGE,DN,NPARTS,
     1  LOCAL,NTERMS,SCALE,P,C,CS,CSINV)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
c
c  purpose:
c
c    this subroutine forms the multipole expansion cause by the nparts
c      particles in the box.
c
c  on input:
c
c    x0y0z0: center of the expansion
c    nparts: number of sources
c    zparts(3,nparts): array of coordinates of sources
c    charge(nparts): array of strengths of sources
c    dn(3,nparts): the normal direction.
c    nparts: the total number of particles.
c    nterms: order of desired expansion
c    scale: the scaling factor.
c    c,cs,csinv: precomputed numbers.
c
c  on output:
c
c    local: coefficients of multipole expansion
c
c  working space :
c
c    p: used for storing the associate legendre polynomials.
c
c  subroutine called : dsqrt(), lgndr()
c  called from : ladapfmm()
c
c  note 1: this subroutine needs the precomputed variables c(,)
c          derived from entry frmini()
c
c       2: the multipole expansion is scaled to avoid over- and
c          under-flow.
c
c       3: only the n=0, ... ,nterms, m=0, ..., nterms coefficients
c          are calculated.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
c
      IMPLICIT NONE
c
      INTEGER *4 NTERMS,NPARTS
      REAL *8 SCALE
      REAL *8 X0Y0Z0(3),ZPARTS(3,NPARTS),CHARGE(NPARTS)
      REAL *8 DN(3,NPARTS)
      REAL *8 P(0:NTERMS,0:NTERMS),C(0:60,0:60)
      REAL *8 CS(0:60,0:60),CSINV(0:60,0:60)
      COMPLEX *16 LOCAL(0:NTERMS,0:NTERMS)
c
c-----local variables.
c
      INTEGER *4 I,L,M,NCALC
      REAL *8 POWERS(0:60)
      REAL *8 PRECIS,D,CP,PROJ,CP1
      REAL *8 RX,RY,RZ,RR,CTHETA,SR2
      COMPLEX *16 IMAG,EPHI(1:60),CPZ,DNZ1,DNZ2
c
      DATA IMAG/(0.0D0,1.0D0)/
      DATA PRECIS/1.0D-14/
c
      DO 4900 I = 1, NPARTS
        RX = ZPARTS(1,I) - X0Y0Z0(1)
        RY = ZPARTS(2,I) - X0Y0Z0(2)
        RZ = ZPARTS(3,I) - X0Y0Z0(3)
c
c-------compute  distance rr
c         ctheta = cos(theta)
c         ephi(1)=e^(i*phi)
c
        PROJ = RX*RX+RY*RY
        RR = PROJ+RZ*RZ
        PROJ = DSQRT(PROJ)
        D = DSQRT(RR)
c
c-------note: here is a hack. when computing cos(theta) as
c             rz/d, we have to be careful about the possibility
c             of d  being 0 (in which case ctheta is not well
c             defined - we arbitrarily set it to 1.)
c
        IF ( D .LE. PRECIS ) THEN
          CTHETA = 1.0D0
        ELSE
          CTHETA = RZ/D
        ENDIF
c
        IF ( PROJ .LE. PRECIS*D ) THEN
          EPHI(1) = 1.0D0
        ELSE
          EPHI(1) = RX/PROJ - IMAG*RY/PROJ
        ENDIF
c
c-------create array of powers of r and powers of e^(-i*phi).
c
        D=1D0/D
        POWERS(0)=1.0D0
        POWERS(1)=D
        D=D/SCALE
        DO 4100 L = 2,NTERMS+2
          EPHI(L) = EPHI(L-1)*EPHI(1)
          POWERS(L) = POWERS(L-1)*D
4100    CONTINUE
c
c-------no contribution from (0,0) term.
c         compute legendre polynomials of argument cos(theta) = ctheta
c         and add contributions from legendre polynomials
c
        CALL LGNDR(NTERMS,CTHETA,P)
c
c-------contribution from (1,0) term.
c
        CP=CHARGE(I)*P(1,0)*POWERS(2)*SCALE
        LOCAL(0,0)=LOCAL(0,0)+CP*DN(3,I)
        DNZ1=DN(1,I)-DN(2,I)*IMAG
        DNZ2=DCONJG(DNZ1)
c
c-------contributions from (l,0), l>1 terms.
c
        DO 4300 L = 2,NTERMS
          CP=CHARGE(I)*P(L,0)*POWERS(L+1)*SCALE
          CP1=CP*CSINV(L-1,1)*CS(L-2,0)/2D0
          LOCAL(L-1,1)=LOCAL(L-1,1)+CP1*DNZ1
          LOCAL(L-1,0)=LOCAL(L-1,0)+CP*DBLE(L)*DN(3,I)
4300    CONTINUE
c
c-------add contributions from associated legendre functions.
c
        DO 4500 L = 1,NTERMS
c
c---------contributions from (l,m), m=1,l-2.
c
          DO 4400 M=1,L-2
            CP=CHARGE(I)*POWERS(L+1)*C(L,M)*P(L,M)*SCALE
            CPZ=CP*EPHI(M)
c
            SR2=CSINV(L,M)*CS(L-1,M-1)
            IF (M.EQ.1) THEN
              CP=SR2*DREAL(CPZ*DNZ2)
              LOCAL(L-1,0)=LOCAL(L-1,0)-CP
            ELSE
              LOCAL(L-1,M-1)=LOCAL(L-1,M-1)-
     1          (SR2/2D0)*CPZ*DNZ2
            ENDIF
c
            SR2=CSINV(L,M)*CS(L-1,M+1)/2D0
            LOCAL(L-1,M+1)=LOCAL(L-1,M+1)+SR2*CPZ*DNZ1
            LOCAL(L-1,M)=LOCAL(L-1,M)+(CS(L-1,M)*CSINV(L,M)*DN(3,I))*
     1        CPZ
4400      CONTINUE
c
c---------m=l-1
c
          M=L-1
          CP=CHARGE(I)*POWERS(L+1)*C(L,M)*P(L,M)*SCALE
          CPZ=CP*EPHI(M)
c
          SR2=CSINV(L,M)*CS(L-1,M-1)
          IF (M.EQ.1) THEN
            CP=SR2*DREAL(CPZ*DNZ2)
            LOCAL(L-1,0)=LOCAL(L-1,0)-CP
          ELSE
            LOCAL(L-1,M-1)=LOCAL(L-1,M-1)-(SR2/2D0)*CPZ*DNZ2
          ENDIF
c
          LOCAL(L-1,M)=LOCAL(L-1,M)+(CS(L-1,M)*CSINV(L,M)*DN(3,I))*
     1      CPZ
c
c---------m=l.
c
          M=L
          CP=CHARGE(I)*POWERS(L+1)*C(L,M)*P(L,M)*SCALE
          CPZ=CP*EPHI(M)
c
          SR2=CSINV(L,M)*CS(L-1,M-1)
          IF (M.EQ.1) THEN
            CP=SR2*DREAL(CPZ*DNZ2)
            LOCAL(L-1,0)=LOCAL(L-1,0)-CP
          ELSE
            LOCAL(L-1,M-1)=LOCAL(L-1,M-1)-(SR2/2D0)*CPZ*DNZ2
          ENDIF
4500    CONTINUE
4900  CONTINUE
c
      RETURN
      END
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
      SUBROUTINE DIRECI(NPARTS,ZPARTS,CHARGE,I,FIELD,POT)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
c
c  purpose:
c
c    calculates directly the potential and field for particle i.
c
c  on input :
c
c    nparts : the total number of particles.
c    zparts(3,:) : the location of the particles.
c    charge : the charge of the particle.
c    i : the location of the target particle.
c
c  on output :
c    rpot : the potential.
c    ptfrc: the field (gradient of potential).
c
c  subroutine called : dsqrt, dexp,
c
c  called from : main()
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
      IMPLICIT NONE
c
      INTEGER *4 NPARTS,I
      REAL *8 ZPARTS(3,NPARTS),CHARGE(NPARTS)
      REAL *8 FIELD(3),POT
c
      INTEGER *4 J
      REAL *8 PI,RX,RY,RZ,RR,RDIS,RMUL
c
      POT = 0.0D0
      FIELD(1) = 0.0D0
      FIELD(2) = 0.0D0
      FIELD(3) = 0.0D0
c
      DO 1000 J = 1,NPARTS
        IF (J .EQ. I) GOTO 1000
        RX = ZPARTS(1,I) - ZPARTS(1,J)
        RY = ZPARTS(2,I) - ZPARTS(2,J)
        RZ = ZPARTS(3,I) - ZPARTS(3,J)
        RR = RX*RX + RY*RY + RZ*RZ
        RDIS = DSQRT(RR)
c
c-------the potential.
c
        POT = POT + CHARGE(J)/RDIS
c
c-------the field.
c
        RMUL = CHARGE(J)/(RDIS*RR)
        FIELD(1) = FIELD(1)+RMUL*RX
        FIELD(2) = FIELD(2)+RMUL*RY
        FIELD(3) = FIELD(3)+RMUL*RZ
1000  CONTINUE
c
      RETURN
      END
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
      SUBROUTINE DNDIRECI(NPARTS,ZPARTS,CHARGE,DN,I,FIELD,POT)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
c
c  purpose:
c
c    calculates directly the potential and field for particle i.
c
c  on input :
c
c    nparts : the total number of particles.
c    zparts(3,:) : the location of the particles.
c    charge : the charge of the particle.
c    dn(3,natoms): the normal direction.
c    i : the location of the target particle.
c
c  on output :
c    rpot : the potential.
c    ptfrc: the field (gradient of potential).
c
c  subroutine called : dsqrt, dexp,
c
c  called from : main()
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
c
      IMPLICIT NONE
c
      INTEGER *4 NPARTS,I
      REAL *8 ZPARTS(3,NPARTS),CHARGE(NPARTS)
      REAL *8 DN(3,NPARTS)
      REAL *8 FIELD(3),POT
c
      INTEGER *4 J
      REAL *8 RX,RY,RZ,RR,RDIS,TERM1,TERM2,TERM3
c
      POT = 0.0D0
      FIELD(1) = 0.0D0
      FIELD(2) = 0.0D0
      FIELD(3) = 0.0D0
c
      DO 1000 J = 1,NPARTS
        IF (J .EQ. I) GOTO 1000
        RX = ZPARTS(1,I) - ZPARTS(1,J)
        RY = ZPARTS(2,I) - ZPARTS(2,J)
        RZ = ZPARTS(3,I) - ZPARTS(3,J)
        RR = RX*RX + RY*RY + RZ*RZ
        RDIS = DSQRT(RR)
c
        TERM1=-CHARGE(J)/(RDIS*RR)
        TERM3=RX*DN(1,J)+RY*DN(2,J)+RZ*DN(3,J)
        TERM2=3D0/RR*TERM1*TERM3
c
c-------the potential.
c
        POT=POT+TERM1*TERM3
c
c-------the field.
c
        FIELD(1)=FIELD(1)-DN(1,J)*TERM1+RX*TERM2
        FIELD(2)=FIELD(2)-DN(2,J)*TERM1+RY*TERM2
        FIELD(3)=FIELD(3)-DN(3,J)*TERM1+RZ*TERM2
1000  CONTINUE
c
      RETURN
      END
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
      SUBROUTINE DIRECIST(NSOU,ZSOU,CHARGE,I,ZTAR,FIELD,POT)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
c
c  purpose:
c
c    calculates directly the potential and field for target particle i.
c
c  on input :
c
c    nsou : the total number of source particles.
c    zsou(3,:) : the location of the source particles.
c    charge : the charge of the particle.
c    i : the index of the target particle.
c    ztar: the location of the target particles.
c
c  on output :
c    rpot : the potential.
c    ptfrc: the field (gradient of potential).
c
c  subroutine called : dsqrt, dexp,
c
c  called from : main()
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
      IMPLICIT NONE
c
      INTEGER *4 NSOU,I
      REAL *8 ZSOU(3,NSOU),CHARGE(NSOU),ZTAR(3,1)
      REAL *8 FIELD(3),POT
c
      INTEGER *4 J
      REAL *8 PI,RX,RY,RZ,RR,RDIS,RMUL
c
      POT = 0.0D0
      FIELD(1) = 0.0D0
      FIELD(2) = 0.0D0
      FIELD(3) = 0.0D0
c
      DO 1000 J = 1,NSOU
        RX = ZTAR(1,I) - ZSOU(1,J)
        RY = ZTAR(2,I) - ZSOU(2,J)
        RZ = ZTAR(3,I) - ZSOU(3,J)
        RR = RX*RX + RY*RY + RZ*RZ
c
c-------stop if two particles are too close.
c
        IF (RR.LE.1D-60) THEN
          PRINT *, 'source and target toooo close, stop'
          STOP
        ENDIF
c
        RDIS = DSQRT(RR)
c
c-------the potential.
c
        POT = POT + CHARGE(J)/RDIS
c
c-------the field.
c
        RMUL = CHARGE(J)/(RDIS*RR)
        FIELD(1) = FIELD(1)+RMUL*RX
        FIELD(2) = FIELD(2)+RMUL*RY
        FIELD(3) = FIELD(3)+RMUL*RZ
1000  CONTINUE
c
      RETURN
      END
c
