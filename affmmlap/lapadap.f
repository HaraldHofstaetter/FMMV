cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
      SUBROUTINE EXPTOLOCAL(LOCAL,NTERMS,RLAMS,WHTS,NLAMBS,
     1  NUMTETS,NTHMAX,NEXPTOT,IEXPU,MEXPUP,IEXPD,MEXPDOWN,CS)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
c
c  purpose:
c    this subroutine converts the fourier representation of two
c    exponential moment functions into a local multipole expansion
c    (with respect to the same box center).
c      l_n^m= (see reference). and is scaled.
c
c    u(x,y,z) = \int_0^\infty e^{-\lambda z}
c                \int_0^{2\pi} e^{i\lambda(xcos(alpha)+ysin(alpha))}
c                mexpup(lambda,alpha) dalpha dlambda
c            +
c                \int_0^\infty e^{\lambda z}
c                \int_0^{2\pi} e^{i\lambda(xcos(alpha)+ysin(alpha))}
c                mexpdown(lambda,alpha) dalpha dlambda
c
c             = \sum_{n=0}^{nterms} \sum_{m=-n,n}
c                local(n,m) y_n^m(cos theta) e^{i m \phi} r^{n}
c
c  on input:
c    nterms : the total number of expansions.
c    iexpu: indicator if mexpup=0 or not.
c    mexpup(nexptot): fourier coefficients of the function
c                    mexpup for discrete lambda
c                    values. they are ordered as follows:
c
c                 mexpup(1,...,numtets(1)) = fourier modes
c                             for lambda_1
c                 mexpup(numtets(1)+1,...,numtets(2)) = fourier modes
c                             for lambda_2
c                 etc.
c
c    iexpd: indicator if mexpdown=0
c    mexpdown(nexptot): as above for down expansion
c
c    rlams(nlambs): discretization points in lambda integral
c    whts(nlambs): quadrature weights in lambda integral
c    nlambs:      number of discretization pts. in lambda integral
c    numtets(j): number of fourier modes in expansion of alpha
c                variable for lambda_j.
c    nthmax:     max_j numtets(j)
c    nexptot:    sum_j numtets(j)
c
c  on output:
c    local(0:nterms,0:nterms): output multipole expansion of order
c                              nterms.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
c
      IMPLICIT NONE
c
      INTEGER *4 NTERMS,NLAMBS,NUMTETS(NLAMBS),NEXPTOT,NTHMAX
      INTEGER *4 IEXPU,IEXPD
c
      REAL *8 RLAMS(NLAMBS),WHTS(NLAMBS)
      REAL *8 CS(0:60,0:60)
c
      COMPLEX *16 LOCAL(0:NTERMS,0:NTERMS)
      COMPLEX *16 MEXPUP(NEXPTOT)
      COMPLEX *16 MEXPDOWN(NEXPTOT)
c
c-----local variables.
c
      INTEGER *4 I,J,NM,NTOT,NL,MMAX,MTH,NCURRENT
      REAL *8 RMUL,RLAMPOW(0:100)
      COMPLEX *16 ZEYE(0:200)
      COMPLEX *16 MEXPPLUS(2000)
      COMPLEX *16 MEXPMINUS(2000)
c
c-----local functions.
c
      COMPLEX *16 DCMPLX
c
c-----compute necessary powers of -i
c
      ZEYE(0) = 1.0D0
      DO I = 1,NTHMAX
        ZEYE(I) = ZEYE(I-1)*DCMPLX(0.0D0,1.0D0)
      ENDDO
c
c-----initialize local expansion
c
      DO NM = 0,NTERMS
        DO MTH = 0,NTERMS
          LOCAL(NM,MTH) = 0.0D0
        ENDDO
      ENDDO
c
c-----compute sum and difference of mexpup and mexpdown
c
      DO NM = 1,NEXPTOT
        IF (IEXPU.LE.0) THEN
          MEXPPLUS(NM) = MEXPDOWN(NM)
          MEXPMINUS(NM) = MEXPDOWN(NM)
        ELSEIF (IEXPD.LE.0) THEN
          MEXPPLUS(NM) =  MEXPUP(NM)
          MEXPMINUS(NM) = -MEXPUP(NM)
        ELSE
          MEXPPLUS(NM) = MEXPDOWN(NM) + MEXPUP(NM)
          MEXPMINUS(NM) = MEXPDOWN(NM) - MEXPUP(NM)
        ENDIF
      ENDDO
c
c-----loop over multipole order to generate mexp values.
c
      NTOT = 1
      DO NL = 1,NLAMBS
c
c-------compute powers of lambda_nl
c
        RLAMPOW(0) = WHTS(NL)
        RMUL = RLAMS(NL)
        DO J = 1,NTERMS
          RLAMPOW(J) = RLAMPOW(J-1)*RMUL
        ENDDO
c
c-------add contributions to local expansion.
c
        DO NM = 0,NTERMS,2
          MMAX = NUMTETS(NL)-1
          IF (MMAX.GT.NM) MMAX = NM
          RMUL = RLAMPOW(NM)
          DO MTH = 0,MMAX
            NCURRENT = NTOT+MTH
            LOCAL(NM,MTH) = LOCAL(NM,MTH) + RMUL*
     1        MEXPPLUS(NCURRENT)
          ENDDO
        ENDDO
        DO NM = 1,NTERMS,2
          MMAX = NUMTETS(NL)-1
          IF (MMAX.GT.NM) MMAX = NM
          RMUL = RLAMPOW(NM)
          DO MTH = 0,MMAX
            NCURRENT = NTOT+MTH
            LOCAL(NM,MTH) = LOCAL(NM,MTH) + RMUL*
     1        MEXPMINUS(NCURRENT)
          ENDDO
        ENDDO
        NTOT = NTOT+NUMTETS(NL)
      ENDDO
c
c-----scale the expansions according to formula
c
      DO NM = 0,NTERMS
        DO MTH = 0,NM
          LOCAL(NM,MTH)=LOCAL(NM,MTH)*ZEYE(MTH)*CS(NM,MTH)
        ENDDO
      ENDDO
c
      RETURN
      END
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
      SUBROUTINE PROCESSUP(SCALE,LEXP1,IUALL,NUALL,IXUALL,
     1  IYUALL,IU1234,NU1234,IX1234,IY1234,MEXUALL,
     2  MEXU1234,XS,YS,ZS,NEXPTOTP,IEXP,IEXP1)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
c
c  purpose:
c
c     this subroutine processes the up interaction lists.
c
c  on input:
c
c     iuall(nuall), ixuall(nuall), iyuall(nuall) are the boxes
c          receiving all child box data and the x and y offsets from
c          child 1, respectively.
c     iu1234(nu1234), ix1234(nu1234), iy1234(nu1234) are the boxes
c          receiving data from child boxes 1,2,3,4 and the x and y
c          offsets from child 1, respectively.
c
c     mexuall is the exponential expansion for all eight children.
c     mexu1234(nexptotp) is the exponential expansion for
c          children 1,2,3,4.
c     xs,ys,zs are the shift coefficients computed by subroutine
c          mkexps.
c     iexp(1-4): indicating if certain source is empty or not.
c     iexp1: indicating if target point is empty or not.
c
c  on output:
c
c     lexp1, which contains the local up expansion information for
c            all boxes, is incremented for each box in the
c            interaction lists.
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
c
      IMPLICIT NONE
c
      INTEGER *4 NEXPTOTP
      INTEGER *4 NUALL, NU1234
      INTEGER *4 IUALL(NUALL),IXUALL(NUALL),IYUALL(NUALL)
      INTEGER *4 IU1234(NU1234),IX1234(NU1234),IY1234(NU1234)
      INTEGER *4 IEXP(4),IEXP1(1)
c
      REAL *8 SCALE
      REAL *8 ZS(3,NEXPTOTP)
c
      COMPLEX *16 LEXP1(NEXPTOTP,1)
      COMPLEX *16 MEXUALL(NEXPTOTP)
      COMPLEX *16 MEXU1234(NEXPTOTP)
      COMPLEX *16 XS(3,NEXPTOTP)
      COMPLEX *16 YS(3,NEXPTOTP)
c
c-----local variables.
c
      INTEGER *4 I,JJ
      COMPLEX *16 ZMUL
c
c-----check to see if source is empty or not.
c
      IF (IEXP(1).GT.0) THEN
        DO I = 1,NUALL
          IEXP1(IUALL(I))=IEXP1(IUALL(I))+1
          DO JJ = 1,NEXPTOTP
            ZMUL = ZS(3,JJ)*SCALE
            IF (IXUALL(I).GT.0)
     1        ZMUL = ZMUL*XS(IXUALL(I),JJ)
            IF (IXUALL(I).LT.0)
     1        ZMUL = ZMUL*DCONJG(XS(-IXUALL(I),JJ))
            IF (IYUALL(I).GT.0)
     1        ZMUL = ZMUL*YS(IYUALL(I),JJ)
            IF (IYUALL(I).LT.0)
     1        ZMUL = ZMUL*DCONJG(YS(-IYUALL(I),JJ))
            LEXP1(JJ,IUALL(I)) = LEXP1(JJ,IUALL(I)) +
     1        MEXUALL(JJ)*ZMUL
          ENDDO
        ENDDO
      ENDIF
c
      IF (IEXP(2).GT.0) THEN
        DO I = 1,NU1234
          IEXP1(IU1234(I))=IEXP1(IU1234(I))+1
          DO JJ = 1,NEXPTOTP
            ZMUL = ZS(2,JJ)*SCALE
            IF (IX1234(I).GT.0)
     1        ZMUL = ZMUL*XS(IX1234(I),JJ)
            IF (IX1234(I).LT.0)
     1        ZMUL = ZMUL*DCONJG(XS(-IX1234(I),JJ))
            IF (IY1234(I).GT.0)
     1        ZMUL = ZMUL*YS(IY1234(I),JJ)
            IF (IY1234(I).LT.0)
     1        ZMUL = ZMUL*DCONJG(YS(-IY1234(I),JJ))
            LEXP1(JJ,IU1234(I)) = LEXP1(JJ,IU1234(I)) +
     1        MEXU1234(JJ)*ZMUL
          ENDDO
        ENDDO
      ENDIF
c
      RETURN
      END
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
      SUBROUTINE PROCESSDN(SCALE,LEXP2,IDALL,NDALL,IXDALL,
     1  IYDALL,ID5678,ND5678,IX5678,IY5678,MEXDALL,
     2  MEXD5678,XS,YS,ZS,NEXPTOTP,IEXP,IEXP1)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
c
c  purpose:
c
c     this subroutine processes the down interaction lists.
c
c  on input:
c
c     idall(ndall), ixdall(ndall), iydall(ndall) are the boxes
c          receiving all child box data and the x and y offsets from
c          child 1, respectively.
c     id5678(nd5678), ix5678(nd5678), iy5678(nd5678) are the boxes
c          receiving data from child boxes 5,6,7,8 and the x and y
c          offsets from child 1, respectively.
c
c     mexdall is the exponential expansion for all eight children.
c     mexd5678(nexptotp) is the exponential expansion for
c          children 5,6,7,8.
c     xs,ys,zs are the shift coefficients computed by subroutine
c          mkexps.
c     iexp: indicates if sources are empty.
c
c  on output:
c
c     lexp2: which contains the local down expansion information for
c            all boxes, is incremented for each box in the
c            interaction lists.
c     iexp1: indicate if lexp2 has been changed.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
c
      IMPLICIT NONE
c
      INTEGER *4 NEXPTOTP
      INTEGER *4 NDALL,ND5678
      INTEGER *4 IDALL(NDALL),IXDALL(NDALL),IYDALL(NDALL)
      INTEGER *4 ID5678(ND5678),IX5678(ND5678),IY5678(ND5678)
      INTEGER *4 IEXP(4),IEXP1(1)
c
      REAL *8 SCALE
      REAL *8 ZS(3,NEXPTOTP)
c
      COMPLEX *16 LEXP2(NEXPTOTP,1)
      COMPLEX *16 MEXDALL(NEXPTOTP)
      COMPLEX *16 MEXD5678(NEXPTOTP)
      COMPLEX *16 XS(3,NEXPTOTP)
      COMPLEX *16 YS(3,NEXPTOTP)
c
c-----local variables.
c
      INTEGER *4 I, JJ
      COMPLEX *16 ZMUL
c
c-----check to see if source is empty or not.
c
      IF (IEXP(3).GT.0) THEN
        DO I = 1,NDALL
          IEXP1(IDALL(I))=IEXP1(IDALL(I))+1
          DO JJ = 1,NEXPTOTP
            ZMUL = ZS(2,JJ)*SCALE
            IF (IXDALL(I).GT.0)
     1        ZMUL = ZMUL*DCONJG(XS(IXDALL(I),JJ))
            IF (IXDALL(I).LT.0)
     1        ZMUL = ZMUL*XS(-IXDALL(I),JJ)
            IF (IYDALL(I).GT.0)
     1        ZMUL = ZMUL*DCONJG(YS(IYDALL(I),JJ))
            IF (IYDALL(I).LT.0)
     1        ZMUL = ZMUL*YS(-IYDALL(I),JJ)
            LEXP2(JJ,IDALL(I)) = LEXP2(JJ,IDALL(I)) +
     1        MEXDALL(JJ)*ZMUL
          ENDDO
        ENDDO
      ENDIF
c
      IF (IEXP(4).GT.0) THEN
        DO I = 1,ND5678
          IEXP1(ID5678(I))=IEXP1(ID5678(I))+1
          DO JJ = 1,NEXPTOTP
            ZMUL = ZS(1,JJ)*SCALE
            IF (IX5678(I).GT.0)
     1        ZMUL = ZMUL*DCONJG(XS(IX5678(I),JJ))
            IF (IX5678(I).LT.0)
     1        ZMUL = ZMUL*XS(-IX5678(I),JJ)
            IF (IY5678(I).GT.0)
     1        ZMUL = ZMUL*DCONJG(YS(IY5678(I),JJ))
            IF (IY5678(I).LT.0)
     1        ZMUL = ZMUL*YS(-IY5678(I),JJ)
            LEXP2(JJ,ID5678(I)) = LEXP2(JJ,ID5678(I)) +
     1        MEXD5678(JJ)*ZMUL
          ENDDO
        ENDDO
      ENDIF
c
      RETURN
      END
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
      SUBROUTINE PROCESSNO(SCALE,LEXP1,INALL,NNALL,IXNALL,IYNALL,
     1  IN1256,NN1256,IX1256,IY1256,
     1  IN12,NN12,IX12,IY12,IN56,NN56,IX56,IY56,
     2  MEXNALL,MEXN1256,MEXN12,MEXN56,XS,YS,ZS,NEXPTOTP,
     3  IEXP,IEXP1)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
c
c  purpose:
c
c     this subroutine processes the north interaction lists.
c
c  on input:
c
c     inall(nnall), ixnall(nnall), iynall(nnall) are the boxes
c          receiving all child box data and the x and y offsets from
c          child 1, respectively.
c     the other north lists are similarly defined (see ymknolist).
c
c     mexnall is the exponential expansion for all eight children, etc.
c     xs,ys,zs are the shift coefficients computed by subroutine
c          mkexps.
c     iexp: indicate if certain sources is empty or not.
c     1: nall; 2: n1256; 3: n12; 4: n56;
c
c  on output:
c
c     lexp1, which contains the local north expansion information for
c            all boxes, is incremented for each box in the
c            interaction lists.
c     iexp1: indicates if the target receives any information.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
c
      IMPLICIT NONE
c
      INTEGER *4 NEXPTOTP
      INTEGER *4 NNALL,NN1256,NN12,NN56
      INTEGER *4 INALL(NNALL),IXNALL(NNALL),IYNALL(NNALL)
      INTEGER *4 IN1256(NN1256),IX1256(NN1256),IY1256(NN1256)
      INTEGER *4 IN12(NN12),IX12(NN12),IY12(NN12)
      INTEGER *4 IN56(NN56),IX56(NN56),IY56(NN56)
      INTEGER *4 IEXP(4),IEXP1(1)
c
      REAL *8 SCALE
      REAL *8 ZS(3,NEXPTOTP)
c
      COMPLEX *16 LEXP1(NEXPTOTP,*)
      COMPLEX *16 MEXNALL(NEXPTOTP)
      COMPLEX *16 MEXN1256(NEXPTOTP)
      COMPLEX *16 MEXN12(NEXPTOTP)
      COMPLEX *16 MEXN56(NEXPTOTP)
      COMPLEX *16 XS(3,NEXPTOTP)
      COMPLEX *16 YS(3,NEXPTOTP)
c
c-----local variables.
c
      INTEGER *4 I, JJ
      COMPLEX *16 ZMUL
c
      IF (IEXP(1).GT.0) THEN
        DO I = 1,NNALL
          IEXP1(INALL(I))=IEXP1(INALL(I))+1
          DO JJ = 1,NEXPTOTP
            ZMUL = ZS(3,JJ)*SCALE
            IF (IXNALL(I).GT.0)
     1        ZMUL = ZMUL*XS(IXNALL(I),JJ)
            IF (IXNALL(I).LT.0)
     1        ZMUL = ZMUL*DCONJG(XS(-IXNALL(I),JJ))
            IF (IYNALL(I).GT.0)
     1        ZMUL = ZMUL*YS(IYNALL(I),JJ)
            IF (IYNALL(I).LT.0)
     1        ZMUL = ZMUL*DCONJG(YS(-IYNALL(I),JJ))
            LEXP1(JJ,INALL(I)) = LEXP1(JJ,INALL(I)) +
     1        MEXNALL(JJ)*ZMUL
          ENDDO
        ENDDO
      ENDIF
c
      IF (IEXP(2).GT.0) THEN
        DO I = 1,NN1256
          IEXP1(IN1256(I))=IEXP1(IN1256(I))+1
          DO JJ = 1,NEXPTOTP
            ZMUL = ZS(2,JJ)*SCALE
            IF (IX1256(I).GT.0)
     1        ZMUL = ZMUL*XS(IX1256(I),JJ)
            IF (IX1256(I).LT.0)
     1        ZMUL = ZMUL*DCONJG(XS(-IX1256(I),JJ))
            IF (IY1256(I).GT.0)
     1        ZMUL = ZMUL*YS(IY1256(I),JJ)
            IF (IY1256(I).LT.0)
     1        ZMUL = ZMUL*DCONJG(YS(-IY1256(I),JJ))
            LEXP1(JJ,IN1256(I)) = LEXP1(JJ,IN1256(I)) +
     1        MEXN1256(JJ)*ZMUL
          ENDDO
        ENDDO
      ENDIF
c
      IF (IEXP(3).GT.0) THEN
        DO I = 1,NN12
          IEXP1(IN12(I))=IEXP1(IN12(I))+1
          DO JJ = 1,NEXPTOTP
            ZMUL = ZS(2,JJ)*SCALE
            IF (IX12(I).GT.0)
     1        ZMUL = ZMUL*XS(IX12(I),JJ)
            IF (IX12(I).LT.0)
     1        ZMUL = ZMUL*DCONJG(XS(-IX12(I),JJ))
            IF (IY12(I).GT.0)
     1        ZMUL = ZMUL*YS(IY12(I),JJ)
            IF (IY12(I).LT.0)
     1        ZMUL = ZMUL*DCONJG(YS(-IY12(I),JJ))
            LEXP1(JJ,IN12(I)) = LEXP1(JJ,IN12(I)) +
     1        MEXN12(JJ)*ZMUL
          ENDDO
        ENDDO
      ENDIF
c
      IF (IEXP(4).GT.0) THEN
        DO I = 1,NN56
          IEXP1(IN56(I))=IEXP1(IN56(I))+1
          DO JJ = 1,NEXPTOTP
            ZMUL = ZS(2,JJ)*SCALE
            IF (IX56(I).GT.0)
     1        ZMUL = ZMUL*XS(IX56(I),JJ)
            IF (IX56(I).LT.0)
     1        ZMUL = ZMUL*DCONJG(XS(-IX56(I),JJ))
            IF (IY56(I).GT.0)
     1        ZMUL = ZMUL*YS(IY56(I),JJ)
            IF (IY56(I).LT.0)
     1        ZMUL = ZMUL*DCONJG(YS(-IY56(I),JJ))
            LEXP1(JJ,IN56(I)) = LEXP1(JJ,IN56(I)) +
     1        MEXN56(JJ)*ZMUL
          ENDDO
        ENDDO
      ENDIF
c
      RETURN
      END
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
      SUBROUTINE PROCESSSO(SCALE,LEXP2,ISALL,NSALL,IXSALL,IYSALL,
     1  IS3478,NS3478,IX3478,IY3478,
     1  IS34,NS34,IX34,IY34,IS78,NS78,IX78,IY78,
     2  MEXSALL,MEXS3478,MEXS34,MEXS78,XS,YS,ZS,NEXPTOTP,
     3  IEXP,IEXP1)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
c
c  purpose:
c
c     this subroutine processes the south interaction lists.
c
c  on input:
c
c     isall(nsall), ixsall(nsall), iysall(nsall) are the boxes
c          receiving all child box data and the x and y offsets from
c          child 1, respectively.
c     the other south lists are similarly defined (see mksolist).
c
c     mexsall is the exponential expansion for all eight children, etc.
c     xs,ys,zs are the shift coefficients computed by subroutine
c          mkexps.
c     iexp: indicate if certain sources is empty or not.
c     5: sall; 6: s3478; 7: s34; 8: s78
c
c  on output:
c
c     lexp2, which contains the local north expansion information for
c            all boxes, is incremented for each box in the
c            interaction lists.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
c
      IMPLICIT NONE
c
      INTEGER *4 NEXPTOTP
      INTEGER *4 NSALL,NS3478,NS34,NS78
      INTEGER *4 ISALL(NSALL),IXSALL(NSALL),IYSALL(NSALL)
      INTEGER *4 IS3478(NS3478),IX3478(NS3478),IY3478(NS3478)
      INTEGER *4 IS34(NS34),IX34(NS34),IY34(NS34)
      INTEGER *4 IS78(NS78),IX78(NS78),IY78(NS78)
      INTEGER *4 IEXP(8),IEXP1(1)
c
      REAL *8 SCALE
      REAL *8 ZS(3,NEXPTOTP)
c
      COMPLEX *16 LEXP2(NEXPTOTP,*)
      COMPLEX *16 MEXSALL(NEXPTOTP)
      COMPLEX *16 MEXS3478(NEXPTOTP)
      COMPLEX *16 MEXS34(NEXPTOTP)
      COMPLEX *16 MEXS78(NEXPTOTP)
      COMPLEX *16 XS(3,NEXPTOTP)
      COMPLEX *16 YS(3,NEXPTOTP)
c
c-----local variables.
c
      INTEGER *4 I,JJ
      COMPLEX *16 ZMUL
c
      IF (IEXP(5).GT.0) THEN
        DO I = 1,NSALL
          IEXP1(ISALL(I))=IEXP1(ISALL(I))+1
          DO JJ = 1,NEXPTOTP
            ZMUL = ZS(2,JJ)*SCALE
            IF (IXSALL(I).GT.0)
     1        ZMUL = ZMUL*DCONJG(XS(IXSALL(I),JJ))
            IF (IXSALL(I).LT.0)
     1        ZMUL = ZMUL*XS(-IXSALL(I),JJ)
            IF (IYSALL(I).GT.0)
     1        ZMUL = ZMUL*DCONJG(YS(IYSALL(I),JJ))
            IF (IYSALL(I).LT.0)
     1        ZMUL = ZMUL*YS(-IYSALL(I),JJ)
            LEXP2(JJ,ISALL(I)) = LEXP2(JJ,ISALL(I)) +
     1        MEXSALL(JJ)*ZMUL
          ENDDO
        ENDDO
      ENDIF
c
      IF (IEXP(6).GT.0) THEN
        DO I = 1,NS3478
          IEXP1(IS3478(I))=IEXP1(IS3478(I))+1
          DO JJ = 1,NEXPTOTP
            ZMUL = ZS(1,JJ)*SCALE
            IF (IX3478(I).GT.0)
     1        ZMUL = ZMUL*DCONJG(XS(IX3478(I),JJ))
            IF (IX3478(I).LT.0)
     1        ZMUL = ZMUL*XS(-IX3478(I),JJ)
            IF (IY3478(I).GT.0)
     1        ZMUL = ZMUL*DCONJG(YS(IY3478(I),JJ))
            IF (IY3478(I).LT.0)
     1        ZMUL = ZMUL*YS(-IY3478(I),JJ)
            LEXP2(JJ,IS3478(I)) = LEXP2(JJ,IS3478(I)) +
     1        MEXS3478(JJ)*ZMUL
          ENDDO
        ENDDO
      ENDIF
c
      IF (IEXP(7).GT.0) THEN
        DO I = 1,NS34
          IEXP1(IS34(I))=IEXP1(IS34(I))+1
          DO JJ = 1,NEXPTOTP
            ZMUL = ZS(1,JJ)*SCALE
            IF (IX34(I).GT.0)
     1        ZMUL = ZMUL*DCONJG(XS(IX34(I),JJ))
            IF (IX34(I).LT.0)
     1        ZMUL = ZMUL*XS(-IX34(I),JJ)
            IF (IY34(I).GT.0)
     1        ZMUL = ZMUL*DCONJG(YS(IY34(I),JJ))
            IF (IY34(I).LT.0)
     1        ZMUL = ZMUL*YS(-IY34(I),JJ)
            LEXP2(JJ,IS34(I)) = LEXP2(JJ,IS34(I)) +
     1        MEXS34(JJ)*ZMUL
          ENDDO
        ENDDO
      ENDIF
c
      IF (IEXP(8).GT.0) THEN
        DO I = 1,NS78
          IEXP1(IS78(I))=IEXP1(IS78(I))+1
          DO JJ = 1,NEXPTOTP
            ZMUL = ZS(1,JJ)*SCALE
            IF (IX78(I).GT.0)
     1        ZMUL = ZMUL*DCONJG(XS(IX78(I),JJ))
            IF (IX78(I).LT.0)
     1        ZMUL = ZMUL*XS(-IX78(I),JJ)
            IF (IY78(I).GT.0)
     1        ZMUL = ZMUL*DCONJG(YS(IY78(I),JJ))
            IF (IY78(I).LT.0)
     1        ZMUL = ZMUL*YS(-IY78(I),JJ)
            LEXP2(JJ,IS78(I)) = LEXP2(JJ,IS78(I)) +
     1        MEXS78(JJ)*ZMUL
          ENDDO
        ENDDO
      ENDIF
c
      RETURN
      END
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
      SUBROUTINE PROCESSEA(SCALE,LEXP1,IEALL,NEALL,IXEALL,IYEALL,
     1  IE1357,NE1357,IX1357,IY1357,IE13,NE13,IX13,IY13,IE57,NE57,
     2  IX57,IY57,IE1,NE1,IX1,IY1,IE3,NE3,IX3,IY3,IE5,NE5,IX5,IY5,
     3  IE7,NE7,IX7,IY7,MEXEALL,MEXE1357,MEXE13,MEXE57,MEXE1,MEXE3,
     4  MEXE5,MEXE7,XS,YS,ZS,NEXPTOTP,IEXP,IEXP1)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
c
c  purpose:
c
c     this subroutine processes the east interaction lists.
c
c  on input:
c
c     ieall(neall), ixeall(neall), iyeall(neall) are the boxes
c          receiving all child box data and the x and y offsets from
c          child 1, respectively.
c     the other east lists are similarly defined (see mkealist).
c
c     mexeall is the exponential expansion for all eight children, etc.
c     xs,ys,zs are the shift coefficients computed by subroutine
c          mkexps.
c     iexp: indicate if certain sources is empty or not.
c     1: eall; 2: e1357; 3: e13; 4: e57; 5: e1; 6: e3; 7: e5; 8: e7;
c     9: wall; 10:w2468; 11:w24;12: w68; 13:w2; 14:w4; 15:w6; 16:w8
c
c  on output:
c
c     lexp1, which contains the local east expansion information for
c            all boxes, is incremented for each box in the
c            interaction lists.
c     iexp1: indicates if the target receives any information.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
c
      IMPLICIT NONE
c
      INTEGER *4 NEXPTOTP
      INTEGER *4 NEALL,NE1357,NE13,NE57,NE1,NE3,NE5,NE7
      INTEGER *4 IEALL(NEALL),IXEALL(NEALL),IYEALL(NEALL)
      INTEGER *4 IE1357(NE1357),IX1357(NE1357),IY1357(NE1357)
      INTEGER *4 IE13(NE13),IX13(NE13),IY13(NE13)
      INTEGER *4 IE57(NE57),IX57(NE57),IY57(NE57)
      INTEGER *4 IE1(NE1),IX1(NE1),IY1(NE1)
      INTEGER *4 IE3(NE3),IX3(NE3),IY3(NE3)
      INTEGER *4 IE5(NE5),IX5(NE5),IY5(NE5)
      INTEGER *4 IE7(NE7),IX7(NE7),IY7(NE7)
      INTEGER *4 IEXP(16),IEXP1(1)
c
      REAL *8 SCALE
      REAL *8 ZS(3,NEXPTOTP)
c
      COMPLEX *16 LEXP1(NEXPTOTP,*)
      COMPLEX *16 MEXEALL(NEXPTOTP)
      COMPLEX *16 MEXE1357(NEXPTOTP)
      COMPLEX *16 MEXE13(NEXPTOTP)
      COMPLEX *16 MEXE57(NEXPTOTP)
      COMPLEX *16 MEXE1(NEXPTOTP)
      COMPLEX *16 MEXE3(NEXPTOTP)
      COMPLEX *16 MEXE5(NEXPTOTP)
      COMPLEX *16 MEXE7(NEXPTOTP)
      COMPLEX *16 XS(3,NEXPTOTP)
      COMPLEX *16 YS(3,NEXPTOTP)
c
c-----local variables.
c
      INTEGER *4 I, JJ
      COMPLEX *16 ZMUL
c
      IF (IEXP(1).GT.0) THEN
        DO I = 1,NEALL
          IEXP1(IEALL(I))=IEXP1(IEALL(I))+1
          DO JJ = 1,NEXPTOTP
            ZMUL = ZS(3,JJ)*SCALE
            IF (IXEALL(I).GT.0)
     1        ZMUL = ZMUL*XS(IXEALL(I),JJ)
            IF (IXEALL(I).LT.0)
     1        ZMUL = ZMUL*DCONJG(XS(-IXEALL(I),JJ))
            IF (IYEALL(I).GT.0)
     1        ZMUL = ZMUL*YS(IYEALL(I),JJ)
            IF (IYEALL(I).LT.0)
     1        ZMUL = ZMUL*DCONJG(YS(-IYEALL(I),JJ))
            LEXP1(JJ,IEALL(I)) = LEXP1(JJ,IEALL(I)) +
     1        MEXEALL(JJ)*ZMUL
          ENDDO
        ENDDO
      ENDIF
c
      IF (IEXP(2).GT.0) THEN
        DO I = 1,NE1357
          IEXP1(IE1357(I))=IEXP1(IE1357(I))+1
          DO JJ = 1,NEXPTOTP
            ZMUL = ZS(2,JJ)*SCALE
            IF (IX1357(I).GT.0)
     1        ZMUL = ZMUL*XS(IX1357(I),JJ)
            IF (IX1357(I).LT.0)
     1        ZMUL = ZMUL*DCONJG(XS(-IX1357(I),JJ))
            IF (IY1357(I).GT.0)
     1        ZMUL = ZMUL*YS(IY1357(I),JJ)
            IF (IY1357(I).LT.0)
     1        ZMUL = ZMUL*DCONJG(YS(-IY1357(I),JJ))
            LEXP1(JJ,IE1357(I)) = LEXP1(JJ,IE1357(I)) +
     1        MEXE1357(JJ)*ZMUL
          ENDDO
        ENDDO
      ENDIF
c
      IF (IEXP(3).GT.0) THEN
        DO I = 1,NE13
          IEXP1(IE13(I))=IEXP1(IE13(I))+1
          DO JJ = 1,NEXPTOTP
            ZMUL = ZS(2,JJ)*SCALE
            IF (IX13(I).GT.0)
     1        ZMUL = ZMUL*XS(IX13(I),JJ)
            IF (IX13(I).LT.0)
     1        ZMUL = ZMUL*DCONJG(XS(-IX13(I),JJ))
            IF (IY13(I).GT.0)
     1        ZMUL = ZMUL*YS(IY13(I),JJ)
            IF (IY13(I).LT.0)
     1        ZMUL = ZMUL*DCONJG(YS(-IY13(I),JJ))
            LEXP1(JJ,IE13(I)) = LEXP1(JJ,IE13(I)) +
     1        MEXE13(JJ)*ZMUL
          ENDDO
        ENDDO
      ENDIF
c
      IF (IEXP(4).GT.0) THEN
        DO I = 1,NE57
          IEXP1(IE57(I))=IEXP1(IE57(I))+1
          DO JJ = 1,NEXPTOTP
            ZMUL = ZS(2,JJ)*SCALE
            IF (IX57(I).GT.0)
     1        ZMUL = ZMUL*XS(IX57(I),JJ)
            IF (IX57(I).LT.0)
     1        ZMUL = ZMUL*DCONJG(XS(-IX57(I),JJ))
            IF (IY57(I).GT.0)
     1        ZMUL = ZMUL*YS(IY57(I),JJ)
            IF (IY57(I).LT.0)
     1        ZMUL = ZMUL*DCONJG(YS(-IY57(I),JJ))
            LEXP1(JJ,IE57(I)) = LEXP1(JJ,IE57(I)) +
     1        MEXE57(JJ)*ZMUL
          ENDDO
        ENDDO
      ENDIF
c
      IF (IEXP(5).GT.0) THEN
        DO I = 1,NE1
          IEXP1(IE1(I))=IEXP1(IE1(I))+1
          DO JJ = 1,NEXPTOTP
            ZMUL = ZS(2,JJ)*SCALE
            IF (IX1(I).GT.0)
     1        ZMUL = ZMUL*XS(IX1(I),JJ)
            IF (IX1(I).LT.0)
     1        ZMUL = ZMUL*DCONJG(XS(-IX1(I),JJ))
            IF (IY1(I).GT.0)
     1        ZMUL = ZMUL*YS(IY1(I),JJ)
            IF (IY1(I).LT.0)
     1        ZMUL = ZMUL*DCONJG(YS(-IY1(I),JJ))
            LEXP1(JJ,IE1(I)) = LEXP1(JJ,IE1(I)) +
     1        MEXE1(JJ)*ZMUL
          ENDDO
        ENDDO
      ENDIF
c
      IF (IEXP(6).GT.0) THEN
        DO I = 1,NE3
          IEXP1(IE3(I))=IEXP1(IE3(I))+1
          DO JJ = 1,NEXPTOTP
            ZMUL = ZS(2,JJ)*SCALE
            IF (IX3(I).GT.0)
     1        ZMUL = ZMUL*XS(IX3(I),JJ)
            IF (IX3(I).LT.0)
     1        ZMUL = ZMUL*DCONJG(XS(-IX3(I),JJ))
            IF (IY3(I).GT.0)
     1        ZMUL = ZMUL*YS(IY3(I),JJ)
            IF (IY3(I).LT.0)
     1        ZMUL = ZMUL*DCONJG(YS(-IY3(I),JJ))
            LEXP1(JJ,IE3(I)) = LEXP1(JJ,IE3(I)) +
     1        MEXE3(JJ)*ZMUL
          ENDDO
        ENDDO
      ENDIF
c
      IF (IEXP(7).GT.0) THEN
        DO I = 1,NE5
          IEXP1(IE5(I))=IEXP1(IE5(I))+1
          DO JJ = 1,NEXPTOTP
            ZMUL = ZS(2,JJ)*SCALE
            IF (IX5(I).GT.0)
     1        ZMUL = ZMUL*XS(IX5(I),JJ)
            IF (IX5(I).LT.0)
     1        ZMUL = ZMUL*DCONJG(XS(-IX5(I),JJ))
            IF (IY5(I).GT.0)
     1        ZMUL = ZMUL*YS(IY5(I),JJ)
            IF (IY5(I).LT.0)
     1        ZMUL = ZMUL*DCONJG(YS(-IY5(I),JJ))
            LEXP1(JJ,IE5(I)) = LEXP1(JJ,IE5(I)) +
     1        MEXE5(JJ)*ZMUL
          ENDDO
        ENDDO
      ENDIF
c
      IF (IEXP(8).GT.0) THEN
        DO I = 1,NE7
          IEXP1(IE7(I))=IEXP1(IE7(I))+1
          DO JJ = 1,NEXPTOTP
            ZMUL = ZS(2,JJ)*SCALE
            IF (IX7(I).GT.0)
     1        ZMUL = ZMUL*XS(IX7(I),JJ)
            IF (IX7(I).LT.0)
     1        ZMUL = ZMUL*DCONJG(XS(-IX7(I),JJ))
            IF (IY7(I).GT.0)
     1        ZMUL = ZMUL*YS(IY7(I),JJ)
            IF (IY7(I).LT.0)
     1        ZMUL = ZMUL*DCONJG(YS(-IY7(I),JJ))
            LEXP1(JJ,IE7(I)) = LEXP1(JJ,IE7(I)) +
     1        MEXE7(JJ)*ZMUL
          ENDDO
        ENDDO
      ENDIF
c
      RETURN
      END
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
      SUBROUTINE PROCESSWE(SCALE,LEXP2,IWALL,NWALL,IXWALL,IYWALL,
     1  IW2468,NW2468,IX2468,IY2468,IW24,NW24,IX24,IY24,IW68,NW68,
     2  IX68,IY68,IW2,NW2,IX2,IY2,IW4,NW4,IX4,IY4,IW6,NW6,IX6,IY6,
     3  IW8,NW8,IX8,IY8,MEXWALL,MEXW2468,MEXW24,MEXW68,MEXW2,MEXW4,
     4  MEXW6,MEXW8,XS,YS,ZS,NEXPTOTP,IEXP,IEXP1)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
c
c  purpose:
c
c     this subroutine processes the west interaction lists.
c
c  on input:
c
c     iwall(nwall), ixwall(nwall), iywall(nwall) are the boxes
c          receiving all child box data and the x and y offsets from
c          child 1, respectively.
c     the other west lists are similarly defined (see mkwelist).
c
c     mexeall is the exponential expansion for all eight children, etc.
c     xs,ys,zs are the shift coefficients computed by subroutine
c          mkexps.
c     iexp: indicate if certain sources is empty or not.
c     1: eall; 2: e1357; 3: e13; 4: e57; 5: e1; 6: e3; 7: e5; 8: e7;
c     9: wall; 10:w2468; 11:w24;12: w68; 13:w2; 14:w4; 15:w6; 16:w8
c
c     iexp1: indicates if the target receives any information.
c
c  on output:
c
c     lexp1, which contains the local west expansion information for
c            all boxes, is incremented for each box in the
c            interaction lists.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
c
      IMPLICIT NONE
c
      INTEGER *4 NEXPTOTP,NWALL,NW2468,NW24,NW68,NW2
      INTEGER *4 NW4,NW6,NW8
      INTEGER *4 IWALL(NWALL),IXWALL(NWALL),IYWALL(NWALL)
      INTEGER *4 IW2468(NW2468),IX2468(NW2468),IY2468(NW2468)
      INTEGER *4 IW24(NW24),IX24(NW24),IY24(NW24)
      INTEGER *4 IW68(NW68),IX68(NW68),IY68(NW68)
      INTEGER *4 IW2(NW2),IX2(NW2),IY2(NW2)
      INTEGER *4 IW4(NW4),IX4(NW4),IY4(NW4)
      INTEGER *4 IW6(NW6),IX6(NW6),IY6(NW6)
      INTEGER *4 IW8(NW8),IX8(NW8),IY8(NW8)
      INTEGER *4 IEXP(16),IEXP1(1)
c
      REAL *8 SCALE
      REAL *8 ZS(3,NEXPTOTP)
c
      COMPLEX *16 LEXP2(NEXPTOTP,*)
      COMPLEX *16 MEXWALL(NEXPTOTP)
      COMPLEX *16 MEXW2468(NEXPTOTP)
      COMPLEX *16 MEXW24(NEXPTOTP)
      COMPLEX *16 MEXW68(NEXPTOTP)
      COMPLEX *16 MEXW2(NEXPTOTP)
      COMPLEX *16 MEXW4(NEXPTOTP)
      COMPLEX *16 MEXW6(NEXPTOTP)
      COMPLEX *16 MEXW8(NEXPTOTP)
      COMPLEX *16 XS(3,NEXPTOTP)
      COMPLEX *16 YS(3,NEXPTOTP)
c
c-----local variables.
c
      INTEGER *4 I,JJ
      COMPLEX *16 ZMUL
c
      IF (IEXP(9).GT.0) THEN
        DO I = 1,NWALL
          IEXP1(IWALL(I))=IEXP1(IWALL(I))+1
          DO JJ = 1,NEXPTOTP
            ZMUL = ZS(2,JJ)*SCALE
            IF (IXWALL(I).GT.0)
     1        ZMUL = ZMUL*DCONJG(XS(IXWALL(I),JJ))
            IF (IXWALL(I).LT.0)
     1        ZMUL = ZMUL*XS(-IXWALL(I),JJ)
            IF (IYWALL(I).GT.0)
     1        ZMUL = ZMUL*DCONJG(YS(IYWALL(I),JJ))
            IF (IYWALL(I).LT.0)
     1        ZMUL = ZMUL*YS(-IYWALL(I),JJ)
            LEXP2(JJ,IWALL(I)) = LEXP2(JJ,IWALL(I)) +
     1        MEXWALL(JJ)*ZMUL
          ENDDO
        ENDDO
      ENDIF
c
      IF (IEXP(10).GT.0) THEN
        DO I = 1,NW2468
          IEXP1(IW2468(I))=IEXP1(IW2468(I))+1
          DO JJ = 1,NEXPTOTP
            ZMUL = ZS(1,JJ)*SCALE
            IF (IX2468(I).GT.0)
     1        ZMUL = ZMUL*DCONJG(XS(IX2468(I),JJ))
            IF (IX2468(I).LT.0)
     1        ZMUL = ZMUL*XS(-IX2468(I),JJ)
            IF (IY2468(I).GT.0)
     1        ZMUL = ZMUL*DCONJG(YS(IY2468(I),JJ))
            IF (IY2468(I).LT.0)
     1        ZMUL = ZMUL*YS(-IY2468(I),JJ)
            LEXP2(JJ,IW2468(I)) = LEXP2(JJ,IW2468(I)) +
     1        MEXW2468(JJ)*ZMUL
          ENDDO
        ENDDO
      ENDIF
c
      IF (IEXP(11).GT.0) THEN
        DO I = 1,NW24
          IEXP1(IW24(I))=IEXP1(IW24(I))+1
          DO JJ = 1,NEXPTOTP
            ZMUL = ZS(1,JJ)*SCALE
            IF (IX24(I).GT.0)
     1        ZMUL = ZMUL*DCONJG(XS(IX24(I),JJ))
            IF (IX24(I).LT.0)
     1        ZMUL = ZMUL*XS(-IX24(I),JJ)
            IF (IY24(I).GT.0)
     1        ZMUL = ZMUL*DCONJG(YS(IY24(I),JJ))
            IF (IY24(I).LT.0)
     1        ZMUL = ZMUL*YS(-IY24(I),JJ)
            LEXP2(JJ,IW24(I)) = LEXP2(JJ,IW24(I)) +
     1        MEXW24(JJ)*ZMUL
          ENDDO
        ENDDO
      ENDIF
c
      IF (IEXP(12).GT.0) THEN
        DO I = 1,NW68
          IEXP1(IW68(I))=IEXP1(IW68(I))+1
          DO JJ = 1,NEXPTOTP
            ZMUL = ZS(1,JJ)*SCALE
            IF (IX68(I).GT.0)
     1        ZMUL = ZMUL*DCONJG(XS(IX68(I),JJ))
            IF (IX68(I).LT.0)
     1        ZMUL = ZMUL*XS(-IX68(I),JJ)
            IF (IY68(I).GT.0)
     1        ZMUL = ZMUL*DCONJG(YS(IY68(I),JJ))
            IF (IY68(I).LT.0)
     1        ZMUL = ZMUL*YS(-IY68(I),JJ)
            LEXP2(JJ,IW68(I)) = LEXP2(JJ,IW68(I)) +
     1        MEXW68(JJ)*ZMUL
          ENDDO
        ENDDO
      ENDIF
c
      IF (IEXP(13).GT.0) THEN
        DO I = 1,NW2
          IEXP1(IW2(I))=IEXP1(IW2(I))+1
          DO JJ = 1,NEXPTOTP
            ZMUL = ZS(1,JJ)*SCALE
            IF (IX2(I).GT.0)
     1        ZMUL = ZMUL*DCONJG(XS(IX2(I),JJ))
            IF (IX2(I).LT.0)
     1        ZMUL = ZMUL*XS(-IX2(I),JJ)
            IF (IY2(I).GT.0)
     1        ZMUL = ZMUL*DCONJG(YS(IY2(I),JJ))
            IF (IY2(I).LT.0)
     1        ZMUL = ZMUL*YS(-IY2(I),JJ)
            LEXP2(JJ,IW2(I)) = LEXP2(JJ,IW2(I)) +
     1        MEXW2(JJ)*ZMUL
          ENDDO
        ENDDO
      ENDIF
c
      IF (IEXP(14).GT.0) THEN
        DO I = 1,NW4
          IEXP1(IW4(I))=IEXP1(IW4(I))+1
          DO JJ = 1,NEXPTOTP
            ZMUL = ZS(1,JJ)*SCALE
            IF (IX4(I).GT.0)
     1        ZMUL = ZMUL*DCONJG(XS(IX4(I),JJ))
            IF (IX4(I).LT.0)
     1        ZMUL = ZMUL*XS(-IX4(I),JJ)
            IF (IY4(I).GT.0)
     1        ZMUL = ZMUL*DCONJG(YS(IY4(I),JJ))
            IF (IY4(I).LT.0)
     1        ZMUL = ZMUL*YS(-IY4(I),JJ)
            LEXP2(JJ,IW4(I)) = LEXP2(JJ,IW4(I)) +
     1        MEXW4(JJ)*ZMUL
          ENDDO
        ENDDO
      ENDIF
c
      IF (IEXP(15).GT.0) THEN
        DO I = 1,NW6
          IEXP1(IW6(I))=IEXP1(IW6(I))+1
          DO JJ = 1,NEXPTOTP
            ZMUL = ZS(1,JJ)*SCALE
            IF (IX6(I).GT.0)
     1        ZMUL = ZMUL*DCONJG(XS(IX6(I),JJ))
            IF (IX6(I).LT.0)
     1        ZMUL = ZMUL*XS(-IX6(I),JJ)
            IF (IY6(I).GT.0)
     1        ZMUL = ZMUL*DCONJG(YS(IY6(I),JJ))
            IF (IY6(I).LT.0)
     1        ZMUL = ZMUL*YS(-IY6(I),JJ)
            LEXP2(JJ,IW6(I)) = LEXP2(JJ,IW6(I)) +
     1        MEXW6(JJ)*ZMUL
          ENDDO
        ENDDO
      ENDIF
c
      IF (IEXP(16).GT.0) THEN
        DO I = 1,NW8
          IEXP1(IW8(I))=IEXP1(IW8(I))+1
          DO JJ = 1,NEXPTOTP
            ZMUL = ZS(1,JJ)*SCALE
            IF (IX8(I).GT.0)
     1        ZMUL = ZMUL*DCONJG(XS(IX8(I),JJ))
            IF (IX8(I).LT.0)
     1        ZMUL = ZMUL*XS(-IX8(I),JJ)
            IF (IY8(I).GT.0)
     1        ZMUL = ZMUL*DCONJG(YS(IY8(I),JJ))
            IF (IY8(I).LT.0)
     1        ZMUL = ZMUL*YS(-IY8(I),JJ)
            LEXP2(JJ,IW8(I)) = LEXP2(JJ,IW8(I)) +
     1        MEXW8(JJ)*ZMUL
          ENDDO
        ENDDO
      ENDIF
c
      RETURN
      END
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
c
c     this following file contains all of the expansion creation
c     routines for a parent box from its four children.
c
c     mkudexp creates all up and down expansions centered at child 1
c
c     mknsexp creates all north and south expansions centered at child 1
c
c     mkewexp creates all east and west expansions centered at child 1
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
      SUBROUTINE MKUDEXP(IBOX,BOX,NTERMS,MPOLE,RLAMS,NLAMBS,
     1  NUMFOUR,NUMPHYS,NTHMAX,NEXPTOT,NEXPTOTP,MEXPUP,
     2  MEXPDOWN,MEXPUPHYS,MEXPDPHYS,MEXUALL,MEXU1234,
     3  MEXDALL,MEXD5678,IEXP,XS,YS,ZS,FEXPE,FEXPO,RLSC)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
c
c  purpose:
c
c     this subroutine creates the up (+z)  and down (-z) exponential
c     expansions for a parent box due to all eight children, if exist.
c
c  note:
c
c     some intelligence is used in the order of summation. thus
c     mexu1234 is computed first and used to initialize mexuall.
c     the contributions from boxes 5 6 7 8  are then added in
c     separately, etc.
c
c  on input:
c
c     ibox: current box number.
c     box: current box information.
c     nterms: number of terms in the multipole expansion.
c     mpole: the multipole expansion coefficients.
c     rlams: exponential expansion coefficients.
c     nlambs: number of terms in the exponential expansion.
c     numfour: number of fourier modes in the expansion.
c     numphys: number of modes in the plane wave expansion.
c     nthmax: max number of terms in the exponential expansion.
c     nexptot: total number of fourier modes in the expansion.
c     nextpotp: half of the fourier modes.
c
c  precomputed tables:
c
c    xs,ys,zs: stores the diagonal translation operators when shifting exponential
c        expansions.
c    fexpe,fexp0: how exponential expansions will be merged.
c    rlsc: stores p_n^m for different lambda_k.
c
c  on output:
c
c    mexuall: up expansion from all boxes.
c    mexu1234: up expansion from box 1-4.
c    mexdall: down expansion from all box.
c    mexd5678: down expansion from 5-8.
c    iexp(1:4): whether any of the output is empty, so we can
c               take advantage of it in the adaptive code.
c
c  variables:
c    mexpup:
c    mexpdown
c    mexpuphys:
c    mexpdphys:
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
c
      IMPLICIT NONE
c
      INTEGER *4 IBOX,BOX(15),NTERMS
      INTEGER *4 NLAMBS,NUMFOUR(NLAMBS),NEXPTOT,NTHMAX
      INTEGER *4 NUMPHYS(NLAMBS),NEXPTOTP
      INTEGER *4 IEXP(4)
c
      REAL *8 ZS(3,NEXPTOTP)
      REAL *8 RLAMS(NLAMBS)
      REAL *8 RLSC(0:NTERMS,0:NTERMS,NLAMBS)
c
      COMPLEX *16 MPOLE(0:NTERMS,0:NTERMS,*)
      COMPLEX *16 MEXPUP(NEXPTOT),MEXPDOWN(NEXPTOT)
      COMPLEX *16 MEXPUPHYS(NEXPTOTP),MEXPDPHYS(NEXPTOTP)
      COMPLEX *16 MEXUALL(NEXPTOTP),MEXU1234(NEXPTOTP)
      COMPLEX *16 MEXDALL(NEXPTOTP),MEXD5678(NEXPTOTP)
      COMPLEX *16 XS(3,NEXPTOTP)
      COMPLEX *16 YS(3,NEXPTOTP)
      COMPLEX *16 FEXPE(1),FEXPO(1)
c
c-----local varialbles:
c
      INTEGER *4 JJ
      COMPLEX *16 ZTMP
c
c-----set all lists to zero list.
c       1: uall; 2: u1234; 3: dall; 4: d5678
c
      IEXP(1)=0
      IEXP(2)=0
      IEXP(3)=0
      IEXP(4)=0
      DO JJ=1,NEXPTOTP
        MEXU1234(JJ)=0.0D0
        MEXD5678(JJ)=0.0D0
        MEXUALL(JJ)=0.0D0
        MEXDALL(JJ)=0.0D0
      ENDDO
c
c-----process the interaction list.
c
      IF (BOX(6).GT.0) THEN
c
c-------child 1 exists, add contributions from child 1
c
        CALL MPOLETOEXP(MPOLE(0,0,BOX(6)),NTERMS,NLAMBS,
     1    NUMFOUR,NEXPTOT,MEXPUP,MEXPDOWN,RLSC)
        CALL FTOPHYS(MEXPUP,NLAMBS,RLAMS,NUMFOUR,NUMPHYS,
     1    NTHMAX,MEXPUPHYS,FEXPE,FEXPO)
        CALL FTOPHYS(MEXPDOWN,NLAMBS,RLAMS,NUMFOUR,NUMPHYS,
     1    NTHMAX,MEXPDPHYS,FEXPE,FEXPO)
        DO JJ = 1,NEXPTOTP
          MEXU1234(JJ) = MEXPUPHYS(JJ)
          MEXDALL(JJ) =  MEXPDPHYS(JJ)
        ENDDO
c
        IEXP(2)=IEXP(2)+1
        IEXP(3)=IEXP(3)+1
c
      ENDIF
c
      IF (BOX(7).GT.0) THEN
c
c-------child 2 exists, add contributions from child 2
c
        CALL MPOLETOEXP(MPOLE(0,0,BOX(7)),NTERMS,NLAMBS,
     1    NUMFOUR,NEXPTOT,MEXPUP,MEXPDOWN,RLSC)
        CALL FTOPHYS(MEXPUP,NLAMBS,RLAMS,NUMFOUR,NUMPHYS,
     1    NTHMAX,MEXPUPHYS,FEXPE,FEXPO)
        CALL FTOPHYS(MEXPDOWN,NLAMBS,RLAMS,NUMFOUR,NUMPHYS,
     1    NTHMAX,MEXPDPHYS,FEXPE,FEXPO)
        DO JJ = 1,NEXPTOTP
          MEXPUPHYS(JJ) = MEXPUPHYS(JJ)*DCONJG(XS(1,JJ))
          MEXU1234(JJ) = MEXU1234(JJ) + MEXPUPHYS(JJ)
          MEXPDPHYS(JJ) = MEXPDPHYS(JJ)*XS(1,JJ)
          MEXDALL(JJ) = MEXDALL(JJ) + MEXPDPHYS(JJ)
        ENDDO
c
        IEXP(2)=IEXP(2)+1
        IEXP(3)=IEXP(3)+1
c
      ENDIF
c
      IF (BOX(8).GT.0) THEN
c
c-------child 3 exists, add contributions from child 3
c
        CALL MPOLETOEXP(MPOLE(0,0,BOX(8)),NTERMS,NLAMBS,
     1    NUMFOUR,NEXPTOT,MEXPUP,MEXPDOWN,RLSC)
        CALL FTOPHYS(MEXPUP,NLAMBS,RLAMS,NUMFOUR,NUMPHYS,
     1    NTHMAX,MEXPUPHYS,FEXPE,FEXPO)
        CALL FTOPHYS(MEXPDOWN,NLAMBS,RLAMS,NUMFOUR,NUMPHYS,
     1    NTHMAX,MEXPDPHYS,FEXPE,FEXPO)
        DO JJ = 1,NEXPTOTP
          MEXPUPHYS(JJ) = MEXPUPHYS(JJ)*DCONJG(YS(1,JJ))
          MEXU1234(JJ) = MEXU1234(JJ) + MEXPUPHYS(JJ)
          MEXPDPHYS(JJ) = MEXPDPHYS(JJ)*YS(1,JJ)
          MEXDALL(JJ) = MEXDALL(JJ) + MEXPDPHYS(JJ)
        ENDDO
c
        IEXP(2)=IEXP(2)+1
        IEXP(3)=IEXP(3)+1
c
      ENDIF
c
      IF (BOX(9).GT.0) THEN
c
c-------child 4 exists, add contributions from child 4
c
        CALL MPOLETOEXP(MPOLE(0,0,BOX(9)),NTERMS,NLAMBS,
     1    NUMFOUR,NEXPTOT,MEXPUP,MEXPDOWN,RLSC)
        CALL FTOPHYS(MEXPUP,NLAMBS,RLAMS,NUMFOUR,NUMPHYS,
     1    NTHMAX,MEXPUPHYS,FEXPE,FEXPO)
        CALL FTOPHYS(MEXPDOWN,NLAMBS,RLAMS,NUMFOUR,NUMPHYS,
     1    NTHMAX,MEXPDPHYS,FEXPE,FEXPO)
        DO JJ = 1,NEXPTOTP
          ZTMP = XS(1,JJ)*YS(1,JJ)
          MEXPUPHYS(JJ) = MEXPUPHYS(JJ)*DCONJG(ZTMP)
          MEXU1234(JJ) = MEXU1234(JJ) + MEXPUPHYS(JJ)
          MEXPDPHYS(JJ) = MEXPDPHYS(JJ)*ZTMP
          MEXDALL(JJ) = MEXDALL(JJ) + MEXPDPHYS(JJ)
        ENDDO
c
        IEXP(2)=IEXP(2)+1
        IEXP(3)=IEXP(3)+1
c
      ENDIF
c
      IF (BOX(10).GT.0) THEN
c
c-------child 5 exists, add contributions from child 5
c
        CALL MPOLETOEXP(MPOLE(0,0,BOX(10)),NTERMS,NLAMBS,
     1    NUMFOUR,NEXPTOT,MEXPUP,MEXPDOWN,RLSC)
        CALL FTOPHYS(MEXPUP,NLAMBS,RLAMS,NUMFOUR,NUMPHYS,
     1    NTHMAX,MEXPUPHYS,FEXPE,FEXPO)
        CALL FTOPHYS(MEXPDOWN,NLAMBS,RLAMS,NUMFOUR,NUMPHYS,
     1    NTHMAX,MEXPDPHYS,FEXPE,FEXPO)
        DO JJ = 1,NEXPTOTP
          MEXPUPHYS(JJ) = MEXPUPHYS(JJ)/ZS(1,JJ)
          MEXUALL(JJ) = MEXU1234(JJ) + MEXPUPHYS(JJ)
          MEXPDPHYS(JJ) = MEXPDPHYS(JJ)*ZS(1,JJ)
          MEXD5678(JJ) = MEXPDPHYS(JJ)
        ENDDO
c
        IEXP(1)=IEXP(1)+1
        IEXP(4)=IEXP(4)+1
      ELSEIF (BOX(10).EQ.0 .AND. IEXP(2).GT.0) THEN
        DO JJ = 1,NEXPTOTP
          MEXUALL(JJ) = MEXU1234(JJ)
        ENDDO
        IEXP(1)=IEXP(1)+1
      ENDIF
c
      IF (BOX(11).GT.0) THEN
c
c-------child 6 exists, add contributions from child 6.
c
        CALL MPOLETOEXP(MPOLE(0,0,BOX(11)),NTERMS,NLAMBS,
     1    NUMFOUR,NEXPTOT,MEXPUP,MEXPDOWN,RLSC)
        CALL FTOPHYS(MEXPUP,NLAMBS,RLAMS,NUMFOUR,NUMPHYS,
     1    NTHMAX,MEXPUPHYS,FEXPE,FEXPO)
        CALL FTOPHYS(MEXPDOWN,NLAMBS,RLAMS,NUMFOUR,NUMPHYS,
     1    NTHMAX,MEXPDPHYS,FEXPE,FEXPO)
        DO JJ = 1,NEXPTOTP
          ZTMP = XS(1,JJ)*ZS(1,JJ)
          MEXPUPHYS(JJ) = MEXPUPHYS(JJ)/ZTMP
          MEXUALL(JJ) = MEXUALL(JJ) + MEXPUPHYS(JJ)
          MEXPDPHYS(JJ) = MEXPDPHYS(JJ)*ZTMP
          MEXD5678(JJ) = MEXD5678(JJ) + MEXPDPHYS(JJ)
        ENDDO
c
        IEXP(1)=IEXP(1)+1
        IEXP(4)=IEXP(4)+1
c
      ENDIF
c
      IF (BOX(12).GT.0) THEN
c
c-------child 7 exists, add contributions from child 7
c
        CALL MPOLETOEXP(MPOLE(0,0,BOX(12)),NTERMS,NLAMBS,
     1    NUMFOUR,NEXPTOT,MEXPUP,MEXPDOWN,RLSC)
        CALL FTOPHYS(MEXPUP,NLAMBS,RLAMS,NUMFOUR,NUMPHYS,
     1    NTHMAX,MEXPUPHYS,FEXPE,FEXPO)
        CALL FTOPHYS(MEXPDOWN,NLAMBS,RLAMS,NUMFOUR,NUMPHYS,
     1    NTHMAX,MEXPDPHYS,FEXPE,FEXPO)
        DO JJ = 1,NEXPTOTP
          ZTMP = YS(1,JJ)*ZS(1,JJ)
          MEXPUPHYS(JJ) = MEXPUPHYS(JJ)/ZTMP
          MEXUALL(JJ) = MEXUALL(JJ) + MEXPUPHYS(JJ)
          MEXPDPHYS(JJ) = MEXPDPHYS(JJ)*ZTMP
          MEXD5678(JJ) = MEXD5678(JJ) + MEXPDPHYS(JJ)
        ENDDO
c
        IEXP(1)=IEXP(1)+1
        IEXP(4)=IEXP(4)+1
c
      ENDIF
c
      IF (BOX(13).GT.0) THEN
c
c-------child 8 exists, add contributions from child 8
c
        CALL MPOLETOEXP(MPOLE(0,0,BOX(13)),NTERMS,NLAMBS,
     1    NUMFOUR,NEXPTOT,MEXPUP,MEXPDOWN,RLSC)
        CALL FTOPHYS(MEXPUP,NLAMBS,RLAMS,NUMFOUR,NUMPHYS,
     1    NTHMAX,MEXPUPHYS,FEXPE,FEXPO)
        CALL FTOPHYS(MEXPDOWN,NLAMBS,RLAMS,NUMFOUR,NUMPHYS,
     1    NTHMAX,MEXPDPHYS,FEXPE,FEXPO)
        DO JJ = 1,NEXPTOTP
          ZTMP = XS(1,JJ)*YS(1,JJ)*ZS(1,JJ)
          MEXPUPHYS(JJ) = MEXPUPHYS(JJ)/ZTMP
          MEXUALL(JJ) = MEXUALL(JJ) + MEXPUPHYS(JJ)
          MEXPDPHYS(JJ) = MEXPDPHYS(JJ)*ZTMP
          MEXD5678(JJ) = MEXD5678(JJ) + MEXPDPHYS(JJ)
          MEXDALL(JJ) = MEXDALL(JJ) + MEXD5678(JJ)
        ENDDO
c
        IEXP(1)=IEXP(1)+1
        IEXP(4)=IEXP(4)+1
        IEXP(3)=IEXP(3)+1
      ELSEIF (BOX(13).EQ.0 .AND. IEXP(4).GT.0) THEN
        DO JJ = 1,NEXPTOTP
          MEXDALL(JJ) = MEXDALL(JJ) + MEXD5678(JJ)
        ENDDO
        IEXP(3)=IEXP(3)+1
      ENDIF
c
      RETURN
      END
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
      SUBROUTINE MKNSEXP(IBOX,BOX,NTERMS,MPOLE,MROTATE,MWORK,
     1  RLAMS,NLAMBS,NUMFOUR,NUMPHYS,NTHMAX,NEXPTOT,NEXPTOTP,
     2  MEXPNOF,MEXPSOF,MEXPNPHYS,MEXPSPHYS,RDMINUS,
     3  MEXNALL,MEXN1256,MEXN12,MEXN56,MEXSALL,MEXS3478,
     4  MEXS34,MEXS78,IEXP,XS,YS,ZS,FEXPE,FEXPO,RLSC)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
c
c  purpose:
c
c     this subroutine creates the north (+y)  and south (-y) exponential
c     expansions for a parent box due to all eight children.
c
c  note:
c
c     some intelligence is used in the order of summation. thus
c     mexn12 and mexn56 are computed separately. mexn1256 is then
c     obtained by adding these two expansions together, etc.
c
c  on input:
c
c     ibox: current box number.
c     box: current box information.
c     nterms: number of terms in the multipole expansion.
c     mpole: the multipole expansion coefficients.
c     rlams: exponential expansion coefficients.
c     nlambs: number of terms in the exponential expansion.
c     numfour: number of fourier modes in the expansion.
c     numphys: number of modes in the plane wave expansion.
c     nthmax: max number of terms in the exponential expansion.
c     nexptot: total number of fourier modes in the expansion.
c     nextpotp: half of the fourier modes.
c
c  precomputed tables:
c    mrotate: rotation matrix so we shift along z-axis.
c    rdminus:
c    xs,ys,zs: stores the diagonal translation operators when shifting exponential
c        expansions.
c    fexpe,fexp0: how exponential expansions will be merged.
c    rlsc: stores p_n^m for different lambda_k.
c
c  on output:
c    mexnall: up expansion from all boxes.
c    mexn1256: up expansion from box 1-4.
c    mexn12:
c    mexn56:
c    mexsall:
c    mexs3478:
c    mexs34:
c    mexs78:
c    iexp(1:8): whether any of the output is empty, so we can
c               take advantage of it in the adaptive code.
c
c  variables:
c    mwork:
c    mexpnof:
c    mexpsof:
c    mexpnphys:
c    mexpsphys:
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
c
      IMPLICIT NONE
c
      INTEGER *4 IBOX,BOX(15),NTERMS
      INTEGER *4 NLAMBS,NUMFOUR(NLAMBS),NEXPTOT,NTHMAX
      INTEGER *4 NUMPHYS(NLAMBS),NEXPTOTP
      INTEGER *4 IEXP(8)
c
      REAL *8 ZS(3,NEXPTOTP)
      REAL *8 RDMINUS(0:NTERMS,0:NTERMS,-NTERMS:NTERMS)
      REAL *8 RLAMS(NLAMBS)
      REAL *8 RLSC(0:NTERMS,0:NTERMS,NLAMBS)
c
      COMPLEX *16 MPOLE(0:NTERMS,0:NTERMS,*)
      COMPLEX *16 MROTATE(0:NTERMS,0:NTERMS)
      COMPLEX *16 MWORK(0:NTERMS,0:NTERMS)
      COMPLEX *16 MEXPNOF(NEXPTOT),MEXPSOF(NEXPTOT)
      COMPLEX *16 MEXPNPHYS(NEXPTOTP),MEXPSPHYS(NEXPTOTP)
      COMPLEX *16 MEXNALL(NEXPTOTP),MEXN1256(NEXPTOTP)
      COMPLEX *16 MEXN12(NEXPTOTP),MEXN56(NEXPTOTP)
      COMPLEX *16 MEXSALL(NEXPTOTP),MEXS3478(NEXPTOTP)
      COMPLEX *16 MEXS34(NEXPTOTP),MEXS78(NEXPTOTP)
      COMPLEX *16 XS(3,NEXPTOTP),YS(3,NEXPTOTP)
      COMPLEX *16 FEXPE(1),FEXPO(1)
c
c-----local variables.
c
      INTEGER *4 JJ
      COMPLEX *16 ZTMP
c
c-----set all lists to zero list.
c     1: nall; 2: n1256; 3: n12; 4: n56;
c     5: sall; 6: s3478; 7: s34; 8: s78
c
      DO JJ=1,8
        IEXP(JJ)=0
      ENDDO
c
      DO JJ=1,NEXPTOTP
        MEXNALL(JJ)=0.0D0
        MEXN1256(JJ)=0.0D0
        MEXN12(JJ)=0.0D0
        MEXN56(JJ)=0.0D0
        MEXSALL(JJ)=0.0D0
        MEXS3478(JJ)=0.0D0
        MEXS34(JJ)=0.0D0
        MEXS78(JJ)=0.0D0
      ENDDO
c
c-----process the interaction list.
c
      IF (BOX(6).GT.0) THEN
c
c-------child 1 exists, add contributions from child 1
c
        CALL ROTZTOY(NTERMS,MPOLE(0,0,BOX(6)),MWORK,MROTATE,RDMINUS)
        CALL MPOLETOEXP(MROTATE,NTERMS,NLAMBS,NUMFOUR,NEXPTOT,
     1    MEXPNOF,MEXPSOF,RLSC)
        CALL FTOPHYS(MEXPNOF,NLAMBS,RLAMS,NUMFOUR,NUMPHYS,
     1    NTHMAX,MEXPNPHYS,FEXPE,FEXPO)
        CALL FTOPHYS(MEXPSOF,NLAMBS,RLAMS,NUMFOUR,NUMPHYS,
     1    NTHMAX,MEXPSPHYS,FEXPE,FEXPO)
        DO JJ = 1,NEXPTOTP
          MEXN12(JJ) = MEXPNPHYS(JJ)
          MEXSALL(JJ) = MEXPSPHYS(JJ)
        ENDDO
        IEXP(3)=IEXP(3)+1
        IEXP(5)=IEXP(5)+1
      ENDIF
c
      IF (BOX(7).GT.0) THEN
c
c-------child 2 exists, add contributions from child 2
c
        CALL ROTZTOY(NTERMS,MPOLE(0,0,BOX(7)),MWORK,MROTATE,RDMINUS)
        CALL MPOLETOEXP(MROTATE,NTERMS,NLAMBS,
     1    NUMFOUR,NEXPTOT,MEXPNOF,MEXPSOF,RLSC)
        CALL FTOPHYS(MEXPNOF,NLAMBS,RLAMS,NUMFOUR,NUMPHYS,
     1    NTHMAX,MEXPNPHYS,FEXPE,FEXPO)
        CALL FTOPHYS(MEXPSOF,NLAMBS,RLAMS,NUMFOUR,NUMPHYS,
     1    NTHMAX,MEXPSPHYS,FEXPE,FEXPO)
        DO JJ = 1,NEXPTOTP
          MEXPNPHYS(JJ) = MEXPNPHYS(JJ)*DCONJG(YS(1,JJ))
          MEXN12(JJ) = MEXN12(JJ) + MEXPNPHYS(JJ)
          MEXPSPHYS(JJ) = MEXPSPHYS(JJ)*YS(1,JJ)
          MEXSALL(JJ) = MEXSALL(JJ) + MEXPSPHYS(JJ)
        ENDDO
        IEXP(3)=IEXP(3)+1
        IEXP(5)=IEXP(5)+1
      ENDIF
c
      IF (BOX(8).GT.0) THEN
c
c-------child 3 exists, add contributions from child 3
c
        CALL ROTZTOY(NTERMS,MPOLE(0,0,BOX(8)),MWORK,MROTATE,RDMINUS)
        CALL MPOLETOEXP(MROTATE,NTERMS,NLAMBS,
     1    NUMFOUR,NEXPTOT,MEXPNOF,MEXPSOF,RLSC)
        CALL FTOPHYS(MEXPNOF,NLAMBS,RLAMS,NUMFOUR,NUMPHYS,
     1    NTHMAX,MEXPNPHYS,FEXPE,FEXPO)
        CALL FTOPHYS(MEXPSOF,NLAMBS,RLAMS,NUMFOUR,NUMPHYS,
     1    NTHMAX,MEXPSPHYS,FEXPE,FEXPO)
        DO JJ = 1,NEXPTOTP
          MEXPNPHYS(JJ) = MEXPNPHYS(JJ)/ZS(1,JJ)
          MEXNALL(JJ) = MEXPNPHYS(JJ)
          MEXPSPHYS(JJ) = MEXPSPHYS(JJ)*ZS(1,JJ)
          MEXS34(JJ) = MEXPSPHYS(JJ)
        ENDDO
        IEXP(1)=IEXP(1)+1
        IEXP(7)=IEXP(7)+1
      ENDIF
c
      IF (BOX(9).GT.0) THEN
c
c-------child 4 exists, add contributions from child 4
c
        CALL ROTZTOY(NTERMS,MPOLE(0,0,BOX(9)),MWORK,MROTATE,RDMINUS)
        CALL MPOLETOEXP(MROTATE,NTERMS,NLAMBS,
     1    NUMFOUR,NEXPTOT,MEXPNOF,MEXPSOF,RLSC)
        CALL FTOPHYS(MEXPNOF,NLAMBS,RLAMS,NUMFOUR,NUMPHYS,
     1    NTHMAX,MEXPNPHYS,FEXPE,FEXPO)
        CALL FTOPHYS(MEXPSOF,NLAMBS,RLAMS,NUMFOUR,NUMPHYS,
     1    NTHMAX,MEXPSPHYS,FEXPE,FEXPO)
        DO JJ = 1,NEXPTOTP
          ZTMP = YS(1,JJ)*ZS(1,JJ)
          MEXPNPHYS(JJ) = MEXPNPHYS(JJ)/ZTMP
          MEXNALL(JJ) = MEXNALL(JJ) + MEXPNPHYS(JJ)
          MEXPSPHYS(JJ) = MEXPSPHYS(JJ)*ZTMP
          MEXS34(JJ) = MEXS34(JJ) + MEXPSPHYS(JJ)
        ENDDO
        IEXP(1)=IEXP(1)+1
        IEXP(7)=IEXP(7)+1
      ENDIF
c
      IF (BOX(10).GT.0) THEN
c
c-------child 5 exists, add contributions from child 5
c
        CALL ROTZTOY(NTERMS,MPOLE(0,0,BOX(10)),MWORK,MROTATE,RDMINUS)
        CALL MPOLETOEXP(MROTATE,NTERMS,NLAMBS,
     1    NUMFOUR,NEXPTOT,MEXPNOF,MEXPSOF,RLSC)
        CALL FTOPHYS(MEXPNOF,NLAMBS,RLAMS,NUMFOUR,NUMPHYS,
     1    NTHMAX,MEXPNPHYS,FEXPE,FEXPO)
        CALL FTOPHYS(MEXPSOF,NLAMBS,RLAMS,NUMFOUR,NUMPHYS,
     1    NTHMAX,MEXPSPHYS,FEXPE,FEXPO)
        DO JJ = 1,NEXPTOTP
          MEXPNPHYS(JJ) = MEXPNPHYS(JJ)*DCONJG(XS(1,JJ))
          MEXN56(JJ) = MEXPNPHYS(JJ)
          MEXPSPHYS(JJ) = MEXPSPHYS(JJ)*XS(1,JJ)
          MEXSALL(JJ) = MEXSALL(JJ) + MEXPSPHYS(JJ)
        ENDDO
        IEXP(4)=IEXP(4)+1
        IEXP(5)=IEXP(5)+1
      ENDIF
c
      IF (BOX(11).GT.0) THEN
c
c-------child 6 exists, add contributions from child 6.
c
        CALL ROTZTOY(NTERMS,MPOLE(0,0,BOX(11)),MWORK,MROTATE,RDMINUS)
        CALL MPOLETOEXP(MROTATE,NTERMS,NLAMBS,
     1    NUMFOUR,NEXPTOT,MEXPNOF,MEXPSOF,RLSC)
        CALL FTOPHYS(MEXPNOF,NLAMBS,RLAMS,NUMFOUR,NUMPHYS,
     1    NTHMAX,MEXPNPHYS,FEXPE,FEXPO)
        CALL FTOPHYS(MEXPSOF,NLAMBS,RLAMS,NUMFOUR,NUMPHYS,
     1    NTHMAX,MEXPSPHYS,FEXPE,FEXPO)
        DO JJ = 1,NEXPTOTP
          ZTMP = YS(1,JJ)*XS(1,JJ)
          MEXPNPHYS(JJ) = MEXPNPHYS(JJ)*DCONJG(ZTMP)
          MEXN56(JJ) = MEXN56(JJ) + MEXPNPHYS(JJ)
          MEXN1256(JJ) = MEXN56(JJ) + MEXN12(JJ)
          MEXPSPHYS(JJ) = MEXPSPHYS(JJ)*ZTMP
          MEXSALL(JJ) = MEXSALL(JJ) + MEXPSPHYS(JJ)
        ENDDO
        IEXP(4)=IEXP(4)+1
        IEXP(2)=IEXP(2)+1
        IEXP(5)=IEXP(5)+1
      ELSEIF (BOX(11).EQ.0 .AND. (IEXP(3)+IEXP(4)).GT.0) THEN
         DO JJ = 1,NEXPTOTP
          MEXN1256(JJ) = MEXN56(JJ) + MEXN12(JJ)
        ENDDO
        IEXP(2)=IEXP(2)+1
      ENDIF
c
      IF (BOX(12).GT.0) THEN
c
c-------child 7 exists, add contributions from child 7
c
        CALL ROTZTOY(NTERMS,MPOLE(0,0,BOX(12)),MWORK,MROTATE,RDMINUS)
        CALL MPOLETOEXP(MROTATE,NTERMS,NLAMBS,
     1    NUMFOUR,NEXPTOT,MEXPNOF,MEXPSOF,RLSC)
        CALL FTOPHYS(MEXPNOF,NLAMBS,RLAMS,NUMFOUR,NUMPHYS,
     1    NTHMAX,MEXPNPHYS,FEXPE,FEXPO)
        CALL FTOPHYS(MEXPSOF,NLAMBS,RLAMS,NUMFOUR,NUMPHYS,
     1    NTHMAX,MEXPSPHYS,FEXPE,FEXPO)
        DO JJ = 1,NEXPTOTP
          ZTMP = ZS(1,JJ)*XS(1,JJ)
          MEXPNPHYS(JJ) = MEXPNPHYS(JJ)/ZTMP
          MEXNALL(JJ) = MEXNALL(JJ) + MEXPNPHYS(JJ)
          MEXPSPHYS(JJ) = MEXPSPHYS(JJ)*ZTMP
          MEXS78(JJ) = MEXPSPHYS(JJ)
        ENDDO
        IEXP(1)=IEXP(1)+1
        IEXP(8)=IEXP(8)+1
      ENDIF
c
      IF (BOX(13).GT.0) THEN
c
c-------child 8 exists, add contributions from child 8
c
        CALL ROTZTOY(NTERMS,MPOLE(0,0,BOX(13)),MWORK,MROTATE,RDMINUS)
        CALL MPOLETOEXP(MROTATE,NTERMS,NLAMBS,
     1    NUMFOUR,NEXPTOT,MEXPNOF,MEXPSOF,RLSC)
        CALL FTOPHYS(MEXPNOF,NLAMBS,RLAMS,NUMFOUR,NUMPHYS,
     1    NTHMAX,MEXPNPHYS,FEXPE,FEXPO)
        CALL FTOPHYS(MEXPSOF,NLAMBS,RLAMS,NUMFOUR,NUMPHYS,
     1    NTHMAX,MEXPSPHYS,FEXPE,FEXPO)
        DO JJ = 1,NEXPTOTP
          ZTMP = ZS(1,JJ)*YS(1,JJ)*XS(1,JJ)
          MEXPNPHYS(JJ) = MEXPNPHYS(JJ)/ZTMP
          MEXNALL(JJ) = MEXNALL(JJ) + MEXPNPHYS(JJ) + MEXN1256(JJ)
          MEXPSPHYS(JJ) = MEXPSPHYS(JJ)*ZTMP
          MEXS78(JJ) = MEXS78(JJ) + MEXPSPHYS(JJ)
          MEXS3478(JJ) = MEXS78(JJ) + MEXS34(JJ)
          MEXSALL(JJ) = MEXSALL(JJ) + MEXS3478(JJ)
        ENDDO
        IEXP(1)=IEXP(1)+1
        IEXP(8)=IEXP(8)+1
        IEXP(6)=IEXP(6)+1
        IEXP(5)=IEXP(5)+1
      ELSEIF (BOX(13).EQ.0 .AND. IEXP(2)+IEXP(7)+IEXP(8).GT.0) THEN
        DO JJ = 1,NEXPTOTP
          MEXNALL(JJ) = MEXNALL(JJ) + MEXN1256(JJ)
          MEXS3478(JJ) = MEXS78(JJ) + MEXS34(JJ)
          MEXSALL(JJ) = MEXSALL(JJ) + MEXS3478(JJ)
        ENDDO
        IEXP(1)=IEXP(1)+1
        IEXP(5)=IEXP(5)+1
        IEXP(6)=IEXP(6)+1
      ENDIF
c
      RETURN
      END
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
      SUBROUTINE MKEWEXP(IBOX,BOX,NTERMS,MPOLE,MROTATE,
     1  RLAMS,NLAMBS,NUMFOUR,NUMPHYS,NTHMAX,NEXPTOT,NEXPTOTP,
     2  MEXPEF,MEXPWF,MEXPEPHYS,MEXPWPHYS,RDPLUS,
     3  MEXEALL,MEXE1357,MEXE13,MEXE57,MEXE1,MEXE3,MEXE5,
     4  MEXE7,MEXWALL,MEXW2468,MEXW24,MEXW68,MEXW2,MEXW4,
     5  MEXW6,MEXW8,IEXP,XS,YS,ZS,FEXPE,FEXPO,RLSC)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
c
c  purpose:
c
c     this subroutine creates the east (+x)  and west (-x) exponential
c     expansions for a parent box due to all eight children.
c
c  note:
c
c     some intelligence is used in the order of summation. thus
c     mexe1, mexe3, mexe5, mexe7 are computed separately. mexe13, mexe57
c     mex1357 are then obtained by adding these two expansions
c     together, etc.
c
c  on input:
c     ibox: current box number.
c     box: current box information.
c     nterms: number of terms in the multipole expansion.
c     mpole: the multipole expansion coefficients.
c     rlams: exponential expansion coefficients.
c     nlambs: number of terms in the exponential expansion.
c     numfour: number of fourier modes in the expansion.
c     numphys: number of modes in the plane wave expansion.
c     nthmax: max number of terms in the exponential expansion.
c     nexptot: total number of fourier modes in the expansion.
c     nextpotp: half of the fourier modes.
c
c  precomputed tables:
c    mrotate: rotation matrix so we shift along z-axis.
c    rdminus:
c    xs,ys,zs: stores the diagonal translation operators when shifting exponential
c        expansions.
c    fexpe,fexp0: how exponential expansions will be merged.
c    rlsc: stores p_n^m for different lambda_k.
c
c  on output:
c    mexeall: east expansion from all boxes.
c    mexe1357: east expansion from box 1-4.
c    mexe13:
c    mexe57:
c    mexe1:
c    mexe3:
c    mexe5:
c    mexe7:
c    mexwall: west expansion from all boxes.
c    mexw2468: west expansion from box 1-4.
c    mexw24:
c    mexw68:
c    mexw2:
c    mexw4:
c    mexw6:
c    mexw8:
c    iexp(1:16): whether any of the output is empty, so we can
c               take advantage of it in the adaptive code.
c
c  variables:
c    mwork:
c    mexpnof:
c    mexpsof:
c    mexpnphys:
c    mexpsphys:
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
c
      IMPLICIT NONE
c
      INTEGER *4 IBOX,BOX(15),NTERMS
      INTEGER *4 NLAMBS,NUMFOUR(NLAMBS),NEXPTOT,NTHMAX
      INTEGER *4 NUMPHYS(NLAMBS),NEXPTOTP
      INTEGER *4 IEXP(16)
c
      REAL *8 RDPLUS(0:NTERMS,0:NTERMS,-NTERMS:NTERMS)
      REAL *8 ZS(3,NEXPTOTP)
      REAL *8 RLAMS(NLAMBS)
      REAL *8 RLSC(0:NTERMS,0:NTERMS,NLAMBS)
c
      COMPLEX *16 MPOLE(0:NTERMS,0:NTERMS,*)
      COMPLEX *16 MROTATE(0:NTERMS,0:NTERMS)
      COMPLEX *16 MEXPEF(NEXPTOT),MEXPWF(NEXPTOT)
      COMPLEX *16 MEXPEPHYS(NEXPTOTP),MEXPWPHYS(NEXPTOTP)
      COMPLEX *16 MEXEALL(NEXPTOTP),MEXE1357(NEXPTOTP)
      COMPLEX *16 MEXE13(NEXPTOTP),MEXE57(NEXPTOTP)
      COMPLEX *16 MEXE1(NEXPTOTP),MEXE3(NEXPTOTP)
      COMPLEX *16 MEXE5(NEXPTOTP),MEXE7(NEXPTOTP)
      COMPLEX *16 MEXWALL(NEXPTOTP),MEXW2468(NEXPTOTP)
      COMPLEX *16 MEXW24(NEXPTOTP),MEXW68(NEXPTOTP)
      COMPLEX *16 MEXW2(NEXPTOTP),MEXW4(NEXPTOTP)
      COMPLEX *16 MEXW6(NEXPTOTP),MEXW8(NEXPTOTP)
      COMPLEX *16 XS(3,NEXPTOTP),YS(3,NEXPTOTP)
      COMPLEX *16 FEXPE(1),FEXPO(1)
c
c-----local variables.
c
      INTEGER *4 JJ
      COMPLEX *16 ZTMP
c
c-----set all lists to zero list.
c     1: eall; 2: e1357; 3: e13; 4: e57; 5: e1; 6: e3; 7: e5; 8: e7;
c     9: wall; 10:w2468; 11:w24;12: w68; 13:w2; 14:w4; 15:w6; 16:w8
c
      DO JJ=1,16
        IEXP(JJ)=0
      ENDDO
c
      DO JJ=1,NEXPTOTP
        MEXEALL(JJ)=0.0D0
        MEXE1357(JJ)=0.0D0
        MEXE13(JJ)=0.0D0
        MEXE57(JJ)=0.0D0
        MEXE1(JJ)=0.0D0
        MEXE3(JJ)=0.0D0
        MEXE5(JJ)=0.0D0
        MEXE7(JJ)=0.0D0
        MEXWALL(JJ)=0.0D0
        MEXW2468(JJ)=0.0D0
        MEXW24(JJ)=0.0D0
        MEXW68(JJ)=0.0D0
        MEXW2(JJ)=0.0D0
        MEXW4(JJ)=0.0D0
        MEXW6(JJ)=0.0D0
        MEXW8(JJ)=0.0D0
      ENDDO
c
c-----process the interaction list.
c
      IF (BOX(6).GT.0) THEN
c
c-------child 1 exists, add contributions from child 1
c
        CALL ROTZTOX(NTERMS,MPOLE(0,0,BOX(6)),MROTATE,RDPLUS)
        CALL MPOLETOEXP(MROTATE,NTERMS,NLAMBS,
     1    NUMFOUR,NEXPTOT,MEXPEF,MEXPWF,RLSC)
        CALL FTOPHYS(MEXPEF,NLAMBS,RLAMS,NUMFOUR,NUMPHYS,
     1    NTHMAX,MEXPEPHYS,FEXPE,FEXPO)
        CALL FTOPHYS(MEXPWF,NLAMBS,RLAMS,NUMFOUR,NUMPHYS,
     1    NTHMAX,MEXPWPHYS,FEXPE,FEXPO)
        DO JJ = 1,NEXPTOTP
          MEXE1(JJ) = MEXPEPHYS(JJ)
          MEXWALL(JJ) = MEXPWPHYS(JJ)
        ENDDO
        IEXP(5)=IEXP(5)+1
        IEXP(9)=IEXP(9)+1
      ENDIF
c
      IF (BOX(7).GT.0) THEN
c
c-------child 2 exists, add contributions from child 2
c
        CALL ROTZTOX(NTERMS,MPOLE(0,0,BOX(7)),MROTATE,RDPLUS)
        CALL MPOLETOEXP(MROTATE,NTERMS,NLAMBS,
     1    NUMFOUR,NEXPTOT,MEXPEF,MEXPWF,RLSC)
        CALL FTOPHYS(MEXPEF,NLAMBS,RLAMS,NUMFOUR,NUMPHYS,
     1    NTHMAX,MEXPEPHYS,FEXPE,FEXPO)
        CALL FTOPHYS(MEXPWF,NLAMBS,RLAMS,NUMFOUR,NUMPHYS,
     1    NTHMAX,MEXPWPHYS,FEXPE,FEXPO)
        DO JJ = 1,NEXPTOTP
          MEXEALL(JJ) = MEXPEPHYS(JJ)/ZS(1,JJ)
          MEXW2(JJ) = MEXPWPHYS(JJ)*ZS(1,JJ)
        ENDDO
        IEXP(1)=IEXP(1)+1
        IEXP(13)=IEXP(13)+1
      ENDIF
c
      IF (BOX(8).GT.0) THEN
c
c-------child 3 exists, add contributions from child 3
c
        CALL ROTZTOX(NTERMS,MPOLE(0,0,BOX(8)),MROTATE,RDPLUS)
        CALL MPOLETOEXP(MROTATE,NTERMS,NLAMBS,
     1    NUMFOUR,NEXPTOT,MEXPEF,MEXPWF,RLSC)
        CALL FTOPHYS(MEXPEF,NLAMBS,RLAMS,NUMFOUR,NUMPHYS,
     1    NTHMAX,MEXPEPHYS,FEXPE,FEXPO)
        CALL FTOPHYS(MEXPWF,NLAMBS,RLAMS,NUMFOUR,NUMPHYS,
     1    NTHMAX,MEXPWPHYS,FEXPE,FEXPO)
        DO JJ = 1,NEXPTOTP
          MEXE3(JJ) = MEXPEPHYS(JJ)*DCONJG(YS(1,JJ))
          MEXE13(JJ) = MEXE1(JJ) + MEXE3(JJ)
          MEXPWPHYS(JJ) = MEXPWPHYS(JJ)*YS(1,JJ)
          MEXWALL(JJ) = MEXWALL(JJ) + MEXPWPHYS(JJ)
        ENDDO
        IEXP(6)=IEXP(6)+1
        IEXP(3)=IEXP(3)+1
        IEXP(9)=IEXP(9)+1
      ELSEIF (BOX(8).EQ.0 .AND. IEXP(5)+IEXP(6).GT.0 ) THEN
        DO JJ = 1,NEXPTOTP
          MEXE13(JJ) = MEXE1(JJ) + MEXE3(JJ)
        ENDDO
        IEXP(3)=IEXP(3)+1
      ENDIF
c
      IF (BOX(9).GT.0) THEN
c
c-------child 4 exists, add contributions from child 4
c
        CALL ROTZTOX(NTERMS,MPOLE(0,0,BOX(9)),MROTATE,RDPLUS)
        CALL MPOLETOEXP(MROTATE,NTERMS,NLAMBS,
     1    NUMFOUR,NEXPTOT,MEXPEF,MEXPWF,RLSC)
        CALL FTOPHYS(MEXPEF,NLAMBS,RLAMS,NUMFOUR,NUMPHYS,
     1    NTHMAX,MEXPEPHYS,FEXPE,FEXPO)
        CALL FTOPHYS(MEXPWF,NLAMBS,RLAMS,NUMFOUR,NUMPHYS,
     1    NTHMAX,MEXPWPHYS,FEXPE,FEXPO)
        DO JJ = 1,NEXPTOTP
          ZTMP = YS(1,JJ)*ZS(1,JJ)
          MEXPEPHYS(JJ) = MEXPEPHYS(JJ)/ZTMP
          MEXEALL(JJ) = MEXEALL(JJ) + MEXPEPHYS(JJ)
          MEXW4(JJ) = MEXPWPHYS(JJ)*ZTMP
          MEXW24(JJ) = MEXW2(JJ) + MEXW4(JJ)
        ENDDO
        IEXP(1)=IEXP(1)+1
        IEXP(11)=IEXP(11)+1
        IEXP(14)=IEXP(14)+1
      ELSEIF (BOX(9).EQ.0 .AND. IEXP(13)+IEXP(14).GT.0) THEN
        DO JJ = 1,NEXPTOTP
          MEXW24(JJ) = MEXW2(JJ) + MEXW4(JJ)
        ENDDO
        IEXP(11)=IEXP(11)+1
      ENDIF
c
      IF (BOX(10).GT.0) THEN
c
c-------child 5 exists, add contributions from child 5
c
        CALL ROTZTOX(NTERMS,MPOLE(0,0,BOX(10)),MROTATE,RDPLUS)
        CALL MPOLETOEXP(MROTATE,NTERMS,NLAMBS,
     1    NUMFOUR,NEXPTOT,MEXPEF,MEXPWF,RLSC)
        CALL FTOPHYS(MEXPEF,NLAMBS,RLAMS,NUMFOUR,NUMPHYS,
     1    NTHMAX,MEXPEPHYS,FEXPE,FEXPO)
        CALL FTOPHYS(MEXPWF,NLAMBS,RLAMS,NUMFOUR,NUMPHYS,
     1    NTHMAX,MEXPWPHYS,FEXPE,FEXPO)
        DO JJ = 1,NEXPTOTP
          MEXE5(JJ) = MEXPEPHYS(JJ)*XS(1,JJ)
          MEXPWPHYS(JJ) = MEXPWPHYS(JJ)*DCONJG(XS(1,JJ))
          MEXWALL(JJ) = MEXWALL(JJ) + MEXPWPHYS(JJ)
        ENDDO
        IEXP(7)=IEXP(7)+1
        IEXP(9)=IEXP(9)+1
      ENDIF
c
      IF (BOX(11).GT.0) THEN
c
c-------child 6 exists, add contributions from child 6.
c
        CALL ROTZTOX(NTERMS,MPOLE(0,0,BOX(11)),MROTATE,RDPLUS)
        CALL MPOLETOEXP(MROTATE,NTERMS,NLAMBS,
     1    NUMFOUR,NEXPTOT,MEXPEF,MEXPWF,RLSC)
        CALL FTOPHYS(MEXPEF,NLAMBS,RLAMS,NUMFOUR,NUMPHYS,
     1    NTHMAX,MEXPEPHYS,FEXPE,FEXPO)
        CALL FTOPHYS(MEXPWF,NLAMBS,RLAMS,NUMFOUR,NUMPHYS,
     1    NTHMAX,MEXPWPHYS,FEXPE,FEXPO)
        DO JJ = 1,NEXPTOTP
          ZTMP = XS(1,JJ)/ZS(1,JJ)
          MEXPEPHYS(JJ) = MEXPEPHYS(JJ)*ZTMP
          MEXEALL(JJ) = MEXEALL(JJ) + MEXPEPHYS(JJ)
          MEXW6(JJ) = MEXPWPHYS(JJ)/ZTMP
        ENDDO
        IEXP(1)=IEXP(1)+1
        IEXP(15)=IEXP(15)+1
      ENDIF
c
      IF (BOX(12).GT.0) THEN
c
c-------child 7 exists, add contributions from child 7
c
        CALL ROTZTOX(NTERMS,MPOLE(0,0,BOX(12)),MROTATE,RDPLUS)
        CALL MPOLETOEXP(MROTATE,NTERMS,NLAMBS,
     1    NUMFOUR,NEXPTOT,MEXPEF,MEXPWF,RLSC)
        CALL FTOPHYS(MEXPEF,NLAMBS,RLAMS,NUMFOUR,NUMPHYS,
     1    NTHMAX,MEXPEPHYS,FEXPE,FEXPO)
        CALL FTOPHYS(MEXPWF,NLAMBS,RLAMS,NUMFOUR,NUMPHYS,
     1    NTHMAX,MEXPWPHYS,FEXPE,FEXPO)
        DO JJ = 1,NEXPTOTP
          ZTMP = XS(1,JJ)*DCONJG(YS(1,JJ))
          MEXE7(JJ) = MEXPEPHYS(JJ)*ZTMP
          MEXE57(JJ) = MEXE5(JJ) + MEXE7(JJ)
          MEXE1357(JJ) = MEXE13(JJ) + MEXE57(JJ)
          MEXPWPHYS(JJ) = MEXPWPHYS(JJ)*DCONJG(ZTMP)
          MEXWALL(JJ) = MEXWALL(JJ) + MEXPWPHYS(JJ)
        ENDDO
        IEXP(2)=IEXP(2)+1
        IEXP(4)=IEXP(4)+1
        IEXP(8)=IEXP(8)+1
        IEXP(9)=IEXP(9)+1
      ELSEIF (BOX(12).EQ.0 .AND. IEXP(7)+IEXP(8)+IEXP(3)+IEXP(4).GT.0)
     1  THEN
        DO JJ = 1,NEXPTOTP
          MEXE57(JJ) = MEXE5(JJ) + MEXE7(JJ)
          MEXE1357(JJ) = MEXE13(JJ) + MEXE57(JJ)
        ENDDO
        IEXP(2)=IEXP(2)+1
        IEXP(4)=IEXP(4)+1
      ENDIF
c
      IF (BOX(13).GT.0) THEN
c
c-------child 8 exists, add contributions from child 8
c
        CALL ROTZTOX(NTERMS,MPOLE(0,0,BOX(13)),MROTATE,RDPLUS)
        CALL MPOLETOEXP(MROTATE,NTERMS,NLAMBS,
     1    NUMFOUR,NEXPTOT,MEXPEF,MEXPWF,RLSC)
        CALL FTOPHYS(MEXPEF,NLAMBS,RLAMS,NUMFOUR,NUMPHYS,
     1    NTHMAX,MEXPEPHYS,FEXPE,FEXPO)
        CALL FTOPHYS(MEXPWF,NLAMBS,RLAMS,NUMFOUR,NUMPHYS,
     1    NTHMAX,MEXPWPHYS,FEXPE,FEXPO)
        DO JJ = 1,NEXPTOTP
          ZTMP = XS(1,JJ)*DCONJG(YS(1,JJ))/ZS(1,JJ)
          MEXPEPHYS(JJ) = MEXPEPHYS(JJ)*ZTMP
          MEXEALL(JJ) = MEXEALL(JJ) + MEXPEPHYS(JJ) + MEXE1357(JJ)
          MEXW8(JJ) = MEXPWPHYS(JJ)/ZTMP
          MEXW68(JJ) = MEXW8(JJ) + MEXW6(JJ)
          MEXW2468(JJ) = MEXW68(JJ) + MEXW24(JJ)
          MEXWALL(JJ) = MEXWALL(JJ) + MEXW2468(JJ)
        ENDDO
        IEXP(1)=IEXP(1)+1
        IEXP(9)=IEXP(9)+1
        IEXP(10)=IEXP(10)+1
        IEXP(12)=IEXP(12)+1
        IEXP(16)=IEXP(16)+1
      ELSEIF (BOX(13).EQ.0 .AND. IEXP(15)+IEXP(16)+IEXP(11)+
     1  IEXP(12)+IEXP(2).GT.0) THEN
        DO JJ = 1,NEXPTOTP
          MEXEALL(JJ) = MEXEALL(JJ) + MEXE1357(JJ)
          MEXW68(JJ) = MEXW8(JJ) + MEXW6(JJ)
          MEXW2468(JJ) = MEXW68(JJ) + MEXW24(JJ)
          MEXWALL(JJ) = MEXWALL(JJ) + MEXW2468(JJ)
        ENDDO
        IEXP(1)=IEXP(1)+1
        IEXP(9)=IEXP(9)+1
        IEXP(10)=IEXP(10)+1
        IEXP(12)=IEXP(12)+1
      ENDIF
c
      RETURN
      END
