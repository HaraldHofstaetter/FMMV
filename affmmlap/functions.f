cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
      SUBROUTINE DUMMY(NMOLS,ZAT,CHARGE)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
c
c  purpose: sets the particle distributions, either random, or on a
c    sphere, or other types for adaptive fmm testing.
c
c  subroutine called : rand()
c
c  called from : main()
c
c  on input :
c
c    nmols : the number of random particles to be generated.
c
c  on output :
c
c    zat(3,:) : the location of the particles.
c    charge : the charge of the particles.
c
c  note: for adaptive code, particles can be located arbitrary.
c        for uniform, only unit box is allowed.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
c
      IMPLICIT NONE
c
      INTEGER *4 NMOLS
      REAL *8 VAL
      REAL *8 ZAT(3,NMOLS),CHARGE(NMOLS)
c
c-----local variables
c
      INTEGER *4 J
      REAL *4 RAND
      REAL *8 THETA, PHI
c
c-----initial seed.
c
      VAL = RAND(1)
c
c-----distribute molecules randomly in box and assign random charges.
c
c     note: for uniform code, the particles are inside the unit box.
c       for nonuniform code, particles can be arbitrary.
c
      DO 1000 J = 1,NMOLS
c
c-------the following generates a distribution on a sphere.
c
ccc        theta=rand(0)*3.1415926d0
ccc        phi=rand(0)*2.0d0*3.1415926d0
ccc        zat(1,j) =1.0d0*dsin(theta)*dsin(phi)
ccc        zat(2,j) =1.0d0*dsin(theta)*dcos(phi)
ccc        zat(3,j) =1.0d0*dcos(theta)
c
c-------the following generates a distribution inside a cylinder.
c
        ZAT(1,J) =0.10D0*( RAND(0)-0.5D0)
        ZAT(2,J) =( RAND(0)-0.5D0)
        ZAT(3,J) =( RAND(0)-0.5D0)
c
c-------generate the charge distribution.
c
        CHARGE(J) = RAND(0) - 0.5D0
1000  CONTINUE
c
      RETURN
      END
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
      SUBROUTINE DNDUMMY(NMOLS,ZAT,CHARGE,DN)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
c
c  purpose: sets the particle distributions, either random, or on a
c    sphere, or other types for adaptive fmm testing.
c
c  subroutine called : rand()
c
c  called from : main()
c
c  on input :
c
c    nmols : the number of random particles to be generated.
c
c  on output :
c
c    zat(3,:) : the location of the particles.
c    charge : the charge of the particles.
c    dn(3,:) : the normal direction.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
c
      IMPLICIT NONE
c
      INTEGER *4 NMOLS
      REAL *8 VAL
      REAL *8 ZAT(3,NMOLS),CHARGE(NMOLS),DN(3,NMOLS)
c
c-----local variables
c
      INTEGER *4 J
      REAL *4 RAND
      REAL *8 THETA,PHI,XX,YY,ZZ,RR
c
c-----initial seed.
c
      VAL = RAND(1)
c
c-----distribute molecules randomly in box and assign random charges.
c
c     note: for uniform code, the particles are inside the unit box.
c       for nonuniform code, particles can be arbitrary.
c
      DO 1000 J = 1,NMOLS
c
c-------the following generates a distribution on a sphere.
c
ccc        theta=rand(0)*3.1415926d0
ccc        phi=rand(0)*2.0d0*3.1415926d0
ccc        zat(1,j) =1.0d0*dsin(theta)*dsin(phi)
ccc        zat(2,j) =1.0d0*dsin(theta)*dcos(phi)
ccc        zat(3,j) =1.0d0*dcos(theta)
c
c-------the following generates a distribution inside a cylinder.
c
        ZAT(1,J) =0.1D0*( RAND(0)-0.5D0)
        ZAT(2,J) =0.1D0*( RAND(0)-0.5D0)
        ZAT(3,J) =( RAND(0)-0.5D0)
c
c-------generate the charge distribution.
c
        CHARGE(J) = RAND(0) - 0.5D0
c
c-------generate the normal derivative vector.
c
        XX =( RAND(0)-0.5D0)
        YY =( RAND(0)-0.5D0)
        ZZ =( RAND(0)-0.5D0)
        RR=DSQRT(XX*XX+YY*YY+ZZ*ZZ)
        DN(1,J)=XX/RR
        DN(2,J)=YY/RR
        DN(3,J)=ZZ/RR
ccc        dn(1,j)=0.0d0
ccc        dn(2,j)=0.0d0
ccc        dn(3,j)=zz/rr
c
1000  CONTINUE
c
      RETURN
      END
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
      SUBROUTINE DUMMYST(NSOU,ZSOU,CHARGE,NTAR,ZTAR)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
c
c  purpose: sets the source and target particle distributions, either
c    random, or on a sphere, or other types for adaptive fmm testing.
c
c  subroutine called : rand()
c
c  called from : main()
c
c  on input :
c
c    nsou : the number of sourc3 particles to be generated.
c    ntar : the number of target particles to be generated.
c
c  on output :
c    zsou(3,:) : the location of source particles.
c    charge : the charge of the source particles.
c    ztar(3,:) : the location of source particles.
c
c  note: for adaptive code, particles can be located arbitrary.
c        for uniform, only unit box is allowed.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
c
      IMPLICIT NONE
c
      INTEGER *4 NSOU,NTAR
      REAL *8 VAL
      REAL *8 ZSOU(3,NSOU),CHARGE(NSOU),ZTAR(3,NTAR)
c
c-----local variables
c
      INTEGER *4 J
      REAL *4 RAND
      REAL *8 THETA, PHI
c
c-----initial seed.
c
      VAL = RAND(1)
c
c-----distribute molecules randomly in box and assign random charges.
c
c     note: for uniform code, the particles are inside the unit box.
c       for nonuniform code, particles can be arbitrary.
c
      DO 1000 J = 1,NSOU
c
c-------the following generates a distribution on a sphere.
c
ccc        theta=rand(0)*3.1415926d0
ccc        phi=rand(0)*2.0d0*3.1415926d0
ccc        zat(1,j) =1.0d0*dsin(theta)*dsin(phi)
ccc        zat(2,j) =1.0d0*dsin(theta)*dcos(phi)
ccc        zat(3,j) =1.0d0*dcos(theta)
c
c-------the following generates a distribution inside a cylinder.
c
        ZSOU(1,J) =1D0*( RAND(0)-0.5D0)
        ZSOU(2,J) =( RAND(0)-0.5D0)
        ZSOU(3,J) =( RAND(0)-0.5D0)
c
c-------generate the charge distribution.
c
        CHARGE(J) = RAND(0) - 0.5D0
1000  CONTINUE
c
      DO 2000 J = 1,NTAR
c
c-------the following generates a distribution on a sphere.
c
ccc        theta=rand(0)*3.1415926d0
ccc        phi=rand(0)*2.0d0*3.1415926d0
ccc        zat(1,j) =1.0d0*dsin(theta)*dsin(phi)
ccc        zat(2,j) =1.0d0*dsin(theta)*dcos(phi)
ccc        zat(3,j) =1.0d0*dcos(theta)
c
c-------the following generates a distribution inside a cylinder.
c
        ZTAR(1,J) =0.1D0*( RAND(0)-0.5D0)
        ZTAR(2,J) =( RAND(0)-0.5D0)
        ZTAR(3,J) =( RAND(0)-0.5D0)
2000  CONTINUE
cccc---debug
c       CALL PRINF('nsou is *', NSOU,1)
c       CALL PRINF('ntar is *', NTAR,1)
c       ZSOU(1,1)=0D0
c       ZSOU(2,1)=0D0
c       ZSOU(3,1)=0D0
c       ZTAR(1,1)=120D0
c       ZTAR(2,1)=0D0
c       ZTAR(3,1)=0D0
c       CHARGE(1)=50D0
cccc----debug
c
      RETURN
      END
c


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
c
c  this file contains common files for both uniform and  adaptive fmm
c    poisson and yukawa solvers.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
      SUBROUTINE ADDEXP(B,A,NTERMS)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
c
c  purpose: adds one vector to another.
c
c  note: this can be replaced by future blast subroutines.
c
c  on input: expansions a and b of size (0:nterms,0:nterms)
c
c  on output: a is over written by (a+b).
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
c
      IMPLICIT NONE
c
      INTEGER *4  NTERMS
c
      COMPLEX *16 B(0:NTERMS,0:NTERMS),A(0:NTERMS,0:NTERMS)
c
c-----local variables.
c
      INTEGER *4 I,J
c
      DO I = 0,NTERMS
        DO J = 0,NTERMS
          A(I,J) = A(I,J) + B(I,J)
        ENDDO
      ENDDO
c
      RETURN
      END
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
      SUBROUTINE PRINM(MPOLE,NTERMS)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
c
c  print out coefficients of multipole expansion
c
c  note: this subroutine is mostly used for debugging purposes.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
c
      IMPLICIT NONE
c
      INTEGER *4 NTERMS
      COMPLEX *16 MPOLE(0:NTERMS,0:NTERMS)
c
      INTEGER *4 L,M
c
      DO 100 L = 0,NTERMS
         WRITE(6,1000)(MPOLE(L,M),M=0,L)
         WRITE(13,1000)(MPOLE(L,M),M=0,L)
         WRITE(6,1001)
         WRITE(13,1001)
100   CONTINUE
1000  FORMAT(6D12.5)
1001  FORMAT(/)
c
      RETURN
      END
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
      SUBROUTINE WRITETIMES(NATOMS,NLEV,NTERMS,TIME_UPPASS,TMKEXPS,
     1  TSHIFTS,TIME_EXPTOMP,TIME_TATA,TIME_LCEVAL)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
c
c  purpose: output the cpu time information for the fmm.
c
c  note: this subroutine is only used by the uniform code.
c    therefore it is ok to put it here.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
c
      IMPLICIT NONE
c
      INTEGER *4 NATOMS,NLEV,NTERMS
      REAL *8 TIME_UPPASS,TMKEXPS,TSHIFTS,TIME_EXPTOMP,TIME_TATA
      REAL *8 TIME_LCEVAL
c
      WRITE(11,*)'n     nlev  nterms up pass mkexps  shifts  ',
     1              'exptomp tata    lceval'
      WRITE(11,1001)NATOMS,NLEV,NTERMS,TIME_UPPASS,TMKEXPS,
     1                TSHIFTS,TIME_EXPTOMP,TIME_TATA,TIME_LCEVAL
1001  FORMAT(I7,2X,I2,4X,I2,1X,6F8.2)
c
      RETURN
      END
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
      SUBROUTINE FSTRTN(NTERMS,D,SQC,THETA)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
c
c  purpose:
c
c    implement the fast version of rotation matrices from
c      the recurrences formulas.
c
c  on input:
c    nterms: an integer indicates the dimension of d.
c    sqc: an array contains the square root of the
c       binormial coefficients.
c    theta:  the rotate angle about the y-axis.
c
c  on output:
c    d: an array which contains the rotation matrix.
c
c  note: only half of d are evaluated, the other
c    half can be obtained by using the symmetricity.
c
c  subroutine called : dabs(), dcos(), dsin(), dsqrt()
c
c  called from : rotgen()
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
c
      IMPLICIT NONE
c
      INTEGER *4 NTERMS
c
      REAL *8 D(0:NTERMS,0:NTERMS,-NTERMS:NTERMS)
      REAL *8 SQC(0:2*NTERMS,0:2*NTERMS)
      REAL *8 THETA
c
c-----local variables.
c
      INTEGER *4 N,M,IJ,IM,IMP,MP,MPABS
      REAL *8 CTHETA, STHETA, HSTHTA
      REAL *8 CTHTAP, CTHTAN, PRECIS
      REAL *8 WW
c
      DATA PRECIS/1.0D-19/
      DATA WW/0.7071067811865476D+00/
c
      CTHETA=DCOS(THETA)
      IF (DABS(CTHETA).LE.PRECIS) CTHETA=0.0D0
      STHETA=DSIN(-THETA)
      IF (DABS(STHETA).LE.PRECIS) STHETA=0.0D0
      HSTHTA=WW*STHETA
      CTHTAP=WW*(1.0D0+CTHETA)
      CTHTAN=-WW*(1.0D0-CTHETA)
c
c-----initial setup for some coefficient matrix.
c
      D(0,0,0)=1.0D0
c
      DO IJ=1,NTERMS
c
c-------compute the result for m'=0 case, use formula (1).
c
        DO IM=-IJ,-1
          D(IJ,0,IM)=-SQC(IJ-IM,2)*D(IJ-1,0,IM+1)
          IF (IM.GT.(1-IJ)) THEN
            D(IJ,0,IM)=D(IJ,0,IM)+SQC(IJ+IM,2)*D(IJ-1,0,IM-1)
          ENDIF
          D(IJ,0,IM)=D(IJ,0,IM)*HSTHTA
          IF (IM.GT.-IJ) THEN
            D(IJ,0,IM)=D(IJ,0,IM)+
     1        D(IJ-1,0,IM)*CTHETA*SQC(IJ+IM,1)*SQC(IJ-IM,1)
          ENDIF
          D(IJ,0,IM)=D(IJ,0,IM)/IJ
        ENDDO
c
        D(IJ,0,0)=D(IJ-1,0,0)*CTHETA
        IF (IJ.GT.1) THEN
          D(IJ,0,0)=D(IJ,0,0)+HSTHTA*SQC(IJ,2)*(D(IJ-1,0,-1)+
     1      D(IJ-1,0,1))/IJ
        ENDIF
c
        DO IM=1,IJ
          D(IJ,0,IM)=-SQC(IJ+IM,2)*D(IJ-1,0,IM-1)
          IF (IM.LT.(IJ-1)) THEN
            D(IJ,0,IM)=D(IJ,0,IM)+SQC(IJ-IM,2)*D(IJ-1,0,IM+1)
          ENDIF
          D(IJ,0,IM)=D(IJ,0,IM)*HSTHTA
          IF (IM.LT.IJ) THEN
            D(IJ,0,IM)=D(IJ,0,IM)+
     1        D(IJ-1,0,IM)*CTHETA*SQC(IJ+IM,1)*SQC(IJ-IM,1)
          ENDIF
          D(IJ,0,IM)=D(IJ,0,IM)/IJ
        ENDDO
c
c-------compute the result for 0<m'<=j case, use formula (2).
c
        DO IMP=1,IJ
          DO IM=-IJ,-1
            D(IJ,IMP,IM)=D(IJ-1,IMP-1,IM+1)*CTHTAN*SQC(IJ-IM,2)
            IF (IM.GT.(1-IJ)) THEN
              D(IJ,IMP,IM)=D(IJ,IMP,IM)-
     1          D(IJ-1,IMP-1,IM-1)*CTHTAP*SQC(IJ+IM,2)
            ENDIF
            IF (IM.GT.-IJ) THEN
              D(IJ,IMP,IM)=D(IJ,IMP,IM)+
     1          D(IJ-1,IMP-1,IM)*STHETA*SQC(IJ+IM,1)*SQC(IJ-IM,1)
            ENDIF
            D(IJ,IMP,IM)=D(IJ,IMP,IM)*WW/SQC(IJ+IMP,2)
          ENDDO
c
          D(IJ,IMP,0)=IJ*STHETA*D(IJ-1,IMP-1,0)
          IF (IJ.GT.1) THEN
            D(IJ,IMP,0)=D(IJ,IMP,0)-SQC(IJ,2)*(
     1        D(IJ-1,IMP-1,-1)*CTHTAP+D(IJ-1,IMP-1,1)*CTHTAN)
          ENDIF
          D(IJ,IMP,0)=D(IJ,IMP,0)*WW/SQC(IJ+IMP,2)
c
          DO IM=1,IJ
            D(IJ,IMP,IM)=D(IJ-1,IMP-1,IM-1)*CTHTAP*SQC(IJ+IM,2)
            IF (IM.LT.(IJ-1)) THEN
              D(IJ,IMP,IM)=D(IJ,IMP,IM)-
     1          D(IJ-1,IMP-1,IM+1)*CTHTAN*SQC(IJ-IM,2)
            ENDIF
            IF (IM.LT.IJ) THEN
              D(IJ,IMP,IM)=D(IJ,IMP,IM)+
     1          D(IJ-1,IMP-1,IM)*STHETA*SQC(IJ+IM,1)*SQC(IJ-IM,1)
            ENDIF
            D(IJ,IMP,IM)=D(IJ,IMP,IM)*WW/SQC(IJ+IMP,2)
          ENDDO
c
c---------use symmetricity, i.e. formula (3.80) in biedenharn & louck's
c           book, to compute the lower part of the matrix
c
        ENDDO
      ENDDO
c
      RETURN
      END
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
      SUBROUTINE YHFSTRTN(NTERMS,D,SQC,THETA)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
c
c  purpose:
c
c    implement the fast version of rotation matrices from
c      the recurrences formulas.
c
c  on input:
c    nterms: an integer indicates the dimension of d.
c    sqc: an array contains the square root of the
c       binormial coefficients.
c    theta:  the rotate angle about the y-axis.
c
c  on output:
c    d: an array which contains the rotation matrix.
c
c  note: only half of d are evaluated, the other
c    half can be obtained by using the symmetricity.
c
c  subroutine called : dabs(), dcos(), dsin(), dsqrt()
c
c  called from : rotgen()
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
c
      IMPLICIT NONE
c
      INTEGER *4 NTERMS
c
      REAL *8 D(0:NTERMS,0:NTERMS,-NTERMS:NTERMS)
      REAL *8 SQC(0:4*NTERMS, 0:4*NTERMS)
      REAL *8 THETA
c
c-----local variables.
c
      INTEGER *4 N,M,IJ,IM,IMP,MP,MPABS
      REAL *8 CTHETA, STHETA, HSTHTA
      REAL *8 CTHTAP, CTHTAN, PRECIS
      REAL *8 WW
      REAL *8 FACTS(0:200)
c
      DATA PRECIS/1.0D-19/
      DATA WW/0.7071067811865476D+00/
c
      CTHETA=DCOS(THETA)
      IF (DABS(CTHETA).LE.PRECIS) CTHETA=0.0D0
      STHETA=DSIN(-THETA)
      IF (DABS(STHETA).LE.PRECIS) STHETA=0.0D0
      HSTHTA=WW*STHETA
      CTHTAP=WW*(1.0D0+CTHETA)
      CTHTAN=-WW*(1.0D0-CTHETA)
c
c-----initial setup for some coefficient matrix.
c
      D(0,0,0)=1.0D0
c
      DO IJ=1,NTERMS
c
c-------compute the result for m'=0 case, use formula (1).
c
        DO IM=-IJ,-1
          D(IJ,0,IM)=-SQC(IJ-IM,2)*D(IJ-1,0,IM+1)
          IF (IM.GT.(1-IJ)) THEN
            D(IJ,0,IM)=D(IJ,0,IM)+SQC(IJ+IM,2)*D(IJ-1,0,IM-1)
          ENDIF
          D(IJ,0,IM)=D(IJ,0,IM)*HSTHTA
          IF (IM.GT.-IJ) THEN
            D(IJ,0,IM)=D(IJ,0,IM)+
     1        D(IJ-1,0,IM)*CTHETA*SQC(IJ+IM,1)*SQC(IJ-IM,1)
          ENDIF
          D(IJ,0,IM)=D(IJ,0,IM)/IJ
        ENDDO
c
        D(IJ,0,0)=D(IJ-1,0,0)*CTHETA
        IF (IJ.GT.1) THEN
          D(IJ,0,0)=D(IJ,0,0)+HSTHTA*SQC(IJ,2)*(D(IJ-1,0,-1)+
     1      D(IJ-1,0,1))/IJ
        ENDIF
c
        DO IM=1,IJ
          D(IJ,0,IM)=-SQC(IJ+IM,2)*D(IJ-1,0,IM-1)
          IF (IM.LT.(IJ-1)) THEN
            D(IJ,0,IM)=D(IJ,0,IM)+SQC(IJ-IM,2)*D(IJ-1,0,IM+1)
          ENDIF
          D(IJ,0,IM)=D(IJ,0,IM)*HSTHTA
          IF (IM.LT.IJ) THEN
            D(IJ,0,IM)=D(IJ,0,IM)+
     1        D(IJ-1,0,IM)*CTHETA*SQC(IJ+IM,1)*SQC(IJ-IM,1)
          ENDIF
          D(IJ,0,IM)=D(IJ,0,IM)/IJ
        ENDDO
c
c-------compute the result for 0<m'<=j case, use formula (2).
c
        DO IMP=1,IJ
          DO IM=-IJ,-1
            D(IJ,IMP,IM)=D(IJ-1,IMP-1,IM+1)*CTHTAN*SQC(IJ-IM,2)
            IF (IM.GT.(1-IJ)) THEN
              D(IJ,IMP,IM)=D(IJ,IMP,IM)-
     1          D(IJ-1,IMP-1,IM-1)*CTHTAP*SQC(IJ+IM,2)
            ENDIF
            IF (IM.GT.-IJ) THEN
              D(IJ,IMP,IM)=D(IJ,IMP,IM)+
     1          D(IJ-1,IMP-1,IM)*STHETA*SQC(IJ+IM,1)*SQC(IJ-IM,1)
            ENDIF
            D(IJ,IMP,IM)=D(IJ,IMP,IM)*WW/SQC(IJ+IMP,2)
          ENDDO
c
          D(IJ,IMP,0)=IJ*STHETA*D(IJ-1,IMP-1,0)
          IF (IJ.GT.1) THEN
            D(IJ,IMP,0)=D(IJ,IMP,0)-SQC(IJ,2)*(
     1        D(IJ-1,IMP-1,-1)*CTHTAP+D(IJ-1,IMP-1,1)*CTHTAN)
          ENDIF
          D(IJ,IMP,0)=D(IJ,IMP,0)*WW/SQC(IJ+IMP,2)
c
          DO IM=1,IJ
            D(IJ,IMP,IM)=D(IJ-1,IMP-1,IM-1)*CTHTAP*SQC(IJ+IM,2)
            IF (IM.LT.(IJ-1)) THEN
              D(IJ,IMP,IM)=D(IJ,IMP,IM)-
     1          D(IJ-1,IMP-1,IM+1)*CTHTAN*SQC(IJ-IM,2)
            ENDIF
            IF (IM.LT.IJ) THEN
              D(IJ,IMP,IM)=D(IJ,IMP,IM)+
     1          D(IJ-1,IMP-1,IM)*STHETA*SQC(IJ+IM,1)*SQC(IJ-IM,1)
            ENDIF
            D(IJ,IMP,IM)=D(IJ,IMP,IM)*WW/SQC(IJ+IMP,2)
          ENDDO
c
c---------note: the lower part of the matrix can be computed using
c           symmetricity, i.e. formula (3.80) in biedenharn & louck's
c           book.
c
        ENDDO
      ENDDO
c
c-----now scale the rotation matrix to avoid y_n^m
c       note : since in yukawa, i will use p_n^m instead of
c            y_n^m.
c
      FACTS(0)=1.0D0
      DO 1 N=1, 2*NTERMS
        FACTS(N)=FACTS(N-1)*DBLE(N)
1     CONTINUE
      DO N=0, NTERMS
        DO M=0, N
          DO MP=-N, N
            MPABS=IABS(MP)
            D(N,M,MP)=D(N,M,MP)*DSQRT(FACTS(N+M)/FACTS(N+MPABS)
     1        *FACTS(N-MPABS)/FACTS(N-M) )
          ENDDO
        ENDDO
      ENDDO
c
      RETURN
      END
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
      SUBROUTINE BNLCFT(C, NTERMS)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
c
c  purpose:
c
c    computes the binomial coefficients c_nterms^n, where n=0,1,2,...,nterms.
c
c  on input:
c
c    nterms: an integer indicates the number we are going to choose from.
c
c  on output:
c
c    c:    an array consists of the squre root of the
c                 binomial coefficients.
c
c  note : this is a different version from the laplace equation.
c    since the binomial coefficients are not needed, but
c    only the square root is needed, so we will not store
c    the binomial coefficients. and we will use that space
c    for some must-be computed coefficients at different levels.
c
c  subroutine called : dsqrt()
c
c  called from : rotgen()
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
c
      IMPLICIT NONE
c
      INTEGER *4 NTERMS
      REAL *8 C(0:NTERMS,0:NTERMS)
c
      INTEGER *4 N,M
c
c-----compute c(n,m)
c
      DO N=0,NTERMS
        C(N,0)=1.0D0
      ENDDO
c
      DO M=1,NTERMS
        C(M,M)=1.0D0
        DO N=M+1,NTERMS
          C(N,M)=C(N-1,M)+C(N-1,M-1)
        ENDDO
      ENDDO
c
c-----compute the square root of c(n,m)
c
      DO M=1,NTERMS
        DO N=M+1,NTERMS
          C(N,M)=DSQRT(C(N,M))
        ENDDO
      ENDDO
c
      RETURN
      END
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
      SUBROUTINE LGNDR(NMAX, X, Y)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
c
c  purpose:
c
c    this subroutine computes the lengre polynomial expansion using
c    a recursive expansion.
c
c  on input:
c    nmax: the max number of terms in the expansion.
c    x: where we want to evaluate the expansion.
c
c  on output:
c    y: the function value at x.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
c
      IMPLICIT NONE
c
      INTEGER *4 NMAX
      REAL *8 X,Y(0:NMAX,0:NMAX)
c
c-----local variables.
c
      INTEGER *4 M,N
      REAL *8 U,DSQRT,DBLE
c
      U=-DSQRT(1.0D0-X*X)
      Y(0,0)=1.0D0
      DO 10 M=0, NMAX
        IF (M.GT.0)  Y(M,M)=Y(M-1,M-1)*U*DBLE(2*M-1)
        IF (M.LT.NMAX)  Y(M+1,M)=DBLE(2*M+1)*X*Y(M,M)
        DO 20 N=M+2, NMAX
          Y(N,M)=((2.0D0*DBLE(N)-1.0D0)*X*Y(N-1,M)-DBLE(N+M-1)
     1      *Y(N-2,M)) / DBLE(N-M)
 20     CONTINUE
 10   CONTINUE
c
      RETURN
      END
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
      SUBROUTINE ROTZTOY(NTERMS,MPOLE,MWORK,MROTATE,RDMINUS)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
c
c  purpose:
c
c    the rotation matrix used in the subroutine mknsexp()
c    so the north_south expansions are made the same way
c    as the up-down expansions.
c
c  on input :
c
c    nterms : number of terms in the multipole expansion.
c    mpole : the multipole expansion.
c    rdminus : the rotation matrix generated in subroutine
c      rotgen<- fstrtn()
c
c  output :
c    mrotate : the rotated multiple expansion coefficients.
c
c  working space : mwork().
c
c  called from :  mknsexp()
c
c  subroutine called : none.
c
c  note : this can be further simplified?
c
c     end result         z_new <- y_old
c                        y_new <- x_old
c                        x_new <- z_old
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
c
      IMPLICIT NONE
c
      INTEGER *4 NTERMS
      REAL *8 RDMINUS(0:NTERMS,0:NTERMS,-NTERMS:NTERMS)
      COMPLEX *16 MPOLE(0:NTERMS,0:NTERMS)
      COMPLEX *16 MWORK(0:NTERMS,0:NTERMS)
      COMPLEX *16 MROTATE(0:NTERMS,0:NTERMS)
c
c-----local varialbles.
c
      INTEGER *4 M,L,MP
      COMPLEX *16 EPHI(0:100)
c
c-----a rotation of -pi/2 radians about the z-axis in the
c       original coordinate system.
c
      EPHI(0) = 1.0D0
      DO M=1,NTERMS
        EPHI(M)=EPHI(M-1)*DCMPLX(0.0D0,-1.0D0)
      ENDDO
c
      DO L=0,NTERMS
        DO M=0,L
          MWORK(L,M)=EPHI(M)*MPOLE(L,M)
        ENDDO
      ENDDO
c
c-----a rotation of -pi/2 radians about the y'-axis in the
c       new coordinate system, bringing the +z-axis in line
c       with the original +y axis.
c
      DO L=0,NTERMS
        DO M=0,L
          MROTATE(L,M)=MWORK(L,0)*RDMINUS(L,0,M)
          DO MP=1,L
            MROTATE(L,M)=MROTATE(L,M)+MWORK(L,MP)*RDMINUS(L,MP,M)+
     1        DCONJG(MWORK(L,MP))*RDMINUS(L,MP,-M)
          ENDDO
        ENDDO
      ENDDO
c
      RETURN
      END
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
      SUBROUTINE ROTYTOZ(NTERMS,MPOLE,MWORK,MROTATE,RDPLUS)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
c
c  purpose:
c
c    the rotation matrix used in the subroutine ladapfmm(), yadapfmm()
c      rotate the totated north_south expansions back to the
c      normal direction.
c
c  on input :
c    nterms : number of terms in the multipole expansion.
c    mpole : the local expansion.
c    rdplus : the rotation matrix generated in subroutine
c      rotgen<- fstrtn()
c
c  on output :
c    mrotate : the rotated local expansion coefficients.
c
c  working space : mwork().
c
c  called from :  ladapfmm(), yadapfmm()
c  subroutine called : none.
c
c  note : this can be further simplified?
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
c
      IMPLICIT NONE
c
      INTEGER *4 NTERMS
      REAL *8 RDPLUS(0:NTERMS,0:NTERMS,-NTERMS:NTERMS)
      COMPLEX *16 MPOLE(0:NTERMS,0:NTERMS)
      COMPLEX *16 MWORK(0:NTERMS,0:NTERMS)
      COMPLEX *16 MROTATE(0:NTERMS,0:NTERMS)
c
c-----local variables.
c
      INTEGER *4 L,M,MP
      COMPLEX *16 EPHI(0:100)
c
c----a rotation of pi/2 radians about the y'-axis in the
c      new coordinate system, bringing the +z-axis in line
c      with the original +z axis.
c
      DO L=0,NTERMS
        DO M=0,L
          MWORK(L,M)=MPOLE(L,0)*RDPLUS(L,0,M)
          DO MP=1,L
            MWORK(L,M)=MWORK(L,M)+MPOLE(L,MP)*RDPLUS(L,MP,M)+
     1        DCONJG(MPOLE(L,MP))*RDPLUS(L,MP,-M)
          ENDDO
        ENDDO
      ENDDO
c
c-----a rotation of pi/2 radians about the z-axis in the
c       original coordinate system.
c
      EPHI(0) = 1.0D0
      DO M=1,NTERMS
        EPHI(M)=EPHI(M-1)*DCMPLX(0.0D0,1.0D0)
      ENDDO
c
      DO L=0,NTERMS
        DO M=0,L
          MROTATE(L,M)=EPHI(M)*MWORK(L,M)
        ENDDO
      ENDDO
c
      RETURN
      END
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
      SUBROUTINE ROTZTOX(NTERMS,MPOLE,MROTATE,RD)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
c
c  purpose:
c
c    the rotation matrix used in the subroutine ladapfmm(), yadapfmm(),
c    and mkewexp(). rotate the east-west expansions from or to the
c    normal direction, depending on rd.
c
c  on input :
c
c    nterms : number of terms in the multipole/local expansion.
c    mpole : the multipole/local expansion.
c    rd : the rotation matrix generated in subroutine
c         rotgen<- fstrtn(), it can either be rdplus or rdminus
c
c  on output :
c    mrotate : the rotated local expansion coefficients.
c
c  called from :  mkewexp(), ladapfmm(), yadapfmm()
c
c  subroutine called : none.
c
c  note : this can be further simplified?
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
c
      IMPLICIT NONE
c
      INTEGER *4 NTERMS
      REAL *8 RD(0:NTERMS,0:NTERMS,-NTERMS:NTERMS)
      COMPLEX *16 MPOLE(0:NTERMS,0:NTERMS)
      COMPLEX *16 MROTATE(0:NTERMS,0:NTERMS)
c
c-----local variables.
c
      INTEGER *4 L,M,MP
c
      DO L=0,NTERMS
        DO M=0,L
          MROTATE(L,M)=MPOLE(L,0)*RD(L,0,M)
          DO MP=1,L
            MROTATE(L,M)=MROTATE(L,M)+MPOLE(L,MP)*RD(L,MP,M)+
     1        DCONJG(MPOLE(L,MP))*RD(L,MP,-M)
          ENDDO
        ENDDO
      ENDDO
c
      RETURN
      END
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
c
c  the following files are for the case when both negative and positive
c    terms are stored for the expansion.
c
c  this is important for helmholtz kernels because half the terms are
c    used for the multipole expansion, but all are used for the local
c    expansion.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
         SUBROUTINE LOCADD(B,A,NTERMS)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
c
c  purpose: adds one vector to another.
c
c  note: this subroutine is used by the helmholtz equation solver.
c        as the m in the local expansion is now from -nterms to nterms
c
c  input :  local expansions a and b of size (0:nterms,-nterms:nterms)
c
c  output:     a is over written by (a+b).
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
c
      IMPLICIT NONE
c
      INTEGER *4  NTERMS
      COMPLEX *16 B(0:NTERMS,-NTERMS:NTERMS),A(0:NTERMS,-NTERMS:NTERMS)
c
      INTEGER *4 I,J
c
      DO I = 0,NTERMS
        DO J = -NTERMS,NTERMS
          A(I,J) = A(I,J) + B(I,J)
        ENDDO
      ENDDO
c
      RETURN
      END
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
      SUBROUTINE LOCROTYTOZ(NTERMS,MPOLE,MWORK,MROTATE,RDPLUS)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
c
c  purpose:
c
c    the rotation matrix used in the subroutine hbrfrc()
c      rotate the totated north_south expansions back to the
c      normal direction.
c
c  input :
c    nterms : number of terms in the multipole expansion.
c    mpole : the local expansion.
c    rdplus : the rotation matrix generated in subroutine
c             hrotgen<- hfstrtn()
c
c  output :
c    mrotate : the rotated local expansion coefficients.
c
c  working space : mwork().
c
c  called from :  hbrfrc()
c  subroutine called : none.
c
c  note : this can be further simplified?
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
c
      IMPLICIT NONE
c
      INTEGER *4 NTERMS
      REAL *8 RDPLUS(0:NTERMS,0:NTERMS,-NTERMS:NTERMS)
      COMPLEX *16 MPOLE(0:NTERMS,-NTERMS:NTERMS)
      COMPLEX *16 MWORK(0:NTERMS,-NTERMS:NTERMS)
      COMPLEX *16 MROTATE(0:NTERMS,-NTERMS:NTERMS)
c
c-----local variables.
c
      INTEGER *4 L,M,MP
      COMPLEX *16 EPHI(-100:100)
c
c-----a rotation of pi/2 radians about the y'-axis in the
c       new coordinate system, bringing the +z-axis in line
c       with the original +z axis.
c
      DO L=0,NTERMS
        DO M=-L,L
          MWORK(L,M)=MPOLE(L,0)*RDPLUS(L,0,M)
          DO MP=1,L
            MWORK(L,M)=MWORK(L,M)+MPOLE(L,MP)*RDPLUS(L,MP,M)+
     1        MPOLE(L,-MP)*RDPLUS(L,MP,-M)
          ENDDO
        ENDDO
      ENDDO
c
c-----a rotation of pi/2 radians about the z-axis in the
c       original coordinate system.
c
      EPHI(0) = 1.0D0
      DO M=1,NTERMS
        EPHI(M)=EPHI(M-1)*DCMPLX(0.0D0,1.0D0)
        EPHI(-M)=DCONJG(EPHI(M))
      ENDDO
c
      DO L=0,NTERMS
        DO M=-L,L
          MROTATE(L,M)=EPHI(M)*MWORK(L,M)
        ENDDO
      ENDDO
c
      RETURN
      END
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
      SUBROUTINE LOCROTZTOX(NTERMS,MPOLE,MROTATE,RD)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
c
c  purpose:
c
c    the rotation matrix used in the subroutine hbrfrc(), and hmkewexp()
c      rotate the east-west expansions from or to the
c      normal direction, depending on rd.
c
c  input :
c    nterms : number of terms in the multipole/local expansion.
c    mpole : the multipole/local expansion.
c    rd : the rotation matrix generated in subroutine
c         hrotgen<- hfstrtn(), it can either be rdplus or rdminus
c
c  output :
c    mrotate : the rotated local expansion coefficients.
c
c  called from :  hmkewexp(), hbrfrc()
c
c  subroutine called : none.
c
c  note : this can be further simplified?
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
c
      IMPLICIT NONE
c
      INTEGER *4 NTERMS
      REAL *8 RD(0:NTERMS,0:NTERMS,-NTERMS:NTERMS)
      COMPLEX *16 MPOLE(0:NTERMS,-NTERMS:NTERMS)
      COMPLEX *16 MROTATE(0:NTERMS,-NTERMS:NTERMS)
c
c-----local variables.
c
      INTEGER *4 L,M,MP
c
      DO L=0,NTERMS
        DO M=-NTERMS,L
          MROTATE(L,M)=MPOLE(L,0)*RD(L,0,M)
          DO MP=1,L
            MROTATE(L,M)=MROTATE(L,M)+MPOLE(L,MP)*RD(L,MP,M)+
     1        MPOLE(L,-MP)*RD(L,MP,-M)
          ENDDO
        ENDDO
      ENDDO
c
      RETURN
      END
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
      SUBROUTINE YHFRMINI(NTERMS,C)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
c
c----- initialization entry point - create factorial scaling factors
c
c     ...... not the most stable way of doing this ......
c
c     the vector fact() is only used here. and c(l,m) will used
c       in subroutine form_mp() and entry brtaev()
c
c     note :
c       the current program is changed a little bit from entry to
c       a subroutine, since the result c is needed at several places.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
c
      IMPLICIT NONE
c
      INTEGER *4 NTERMS
      REAL *8 C(0:NTERMS,0:NTERMS)
c
c-----local variables.
c
      INTEGER *4 L,M
      REAL *8 FACT(0:120)
c
c-----functions called
c
      REAL *8 DBLE
c
      PRINT *, ' initializing'
      FACT(0) = 1.0D0
c
      DO 6000 L=1,2*NTERMS+1
        FACT(L) = FACT(L-1)*DBLE(L)
6000  CONTINUE
c
      DO 6200 L=0,NTERMS
        DO 6100 M = 0,L
          C(L,M) = FACT(L-M)/FACT(L+M)*DBLE(2*L+1)
6100    CONTINUE
6200  CONTINUE
c
      RETURN
      END
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
      SUBROUTINE FRMINI(C,CS,CSINV)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
c
c----- initialization entry point - create factorial scaling factors
c
c     ...... not the most stable way of doing this ......
c
c     the vector fact() is only used here. and c(l,m) will used
c       in subroutine form_mp() and entry brtaev()
c       csinv and cs will be used in brtaev only.
c
c     note :
c       the current program is changed a little bit from entry to
c       a subroutine, since the result c is needed at several places.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
c
      IMPLICIT NONE
c
      REAL *8 C(0:60,0:60),CS(0:60,0:60),CSINV(0:60,0:60)
c
c-----local variables.
c
      INTEGER *4 L,M
      REAL *8 D,FACT(0:120)
c
      D = 1.0D0
      FACT(0) = D
      DO 6000 L=1,120
        D=D*DSQRT(L+0.0D0)
        FACT(L) = D
6000  CONTINUE
c
      CS(0,0) = 1.0D0
      CSINV(0,0) = 1.0D0
      DO 6200 L=1,60
        DO 6100 M = 0,L
          C(L,M) = FACT(L-M)/FACT(L+M)
          CSINV(L,M) = FACT(L-M)*FACT(L+M)
          CS(L,M) = 1.0D0/CSINV(L,M)
6100    CONTINUE
6200  CONTINUE
c
      RETURN
      END
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
      SUBROUTINE YHROTGEN(NTERMS,CARRAY,RDPI2,RDMPI2,RDSQ3,RDMSQ3)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
c
c  purpose:
c    precomputes the rotation matrix for
c     1. mp->mp
c     2. local->local
c     3. east-west expansion
c     4. north-south expansion
c
c  on input :
c    nterms : the number of terms in the multipole expansion.
c
c  on output :
c    rdpi2, rdmpi2 : the rotation matrix for 3 and 4.
c    rdsq3, rdmsq3 : the rotation matrix for 1 and 2.
c
c  workspace :
c    carray : the square root of the binomial numbers.
c      these numbers are only used here in this subroutine.
c
c  subroutine called :
c    bnlcft(), fstrtn(), datan, dacos,
c
c  called from : main()
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
c
      IMPLICIT NONE
c
      INTEGER *4 NTERMS
      REAL *8 CARRAY(0:4*NTERMS,0:4*NTERMS)
      REAL *8 RDPI2(0:NTERMS,0:NTERMS,-NTERMS:NTERMS)
      REAL *8 RDMPI2(0:NTERMS,0:NTERMS,-NTERMS:NTERMS)
      REAL *8 RDSQ3(0:NTERMS,0:NTERMS,-NTERMS:NTERMS)
      REAL *8 RDMSQ3(0:NTERMS,0:NTERMS,-NTERMS:NTERMS)
c
c-----local varialbles.
c
      REAL *8 PI,THETA
c
c-----call initialization routines
c
      PI = 4.0D0*DATAN(1.0D0)
c
      CALL BNLCFT(CARRAY,4*NTERMS)
      THETA = PI/2.0D0
      CALL YHFSTRTN(NTERMS,RDPI2,CARRAY,THETA)
      THETA = -PI/2.0D0
      CALL YHFSTRTN(NTERMS,RDMPI2,CARRAY,THETA)
      THETA = DACOS(DSQRT(3.0D0)/3.0D0)
      CALL YHFSTRTN(NTERMS,RDSQ3,CARRAY,THETA)
      THETA = DACOS(-DSQRT(3.0D0)/3.0D0)
      CALL YHFSTRTN(NTERMS,RDMSQ3,CARRAY,THETA)
c
      RETURN
      END
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
      SUBROUTINE ROTGEN(NTERMS,CARRAY,RDPI2,RDMPI2,RDSQ3,RDMSQ3,DC)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
c
c  purpose:
c    precomputes the rotation matrix for
c     1. mp->mp
c     2. local->local
c     3. east-west expansion
c     4. north-south expansion
c
c  on input :
c    nterms : the number of terms in the multipole expansion.
c
c  on output :
c    rdpi2, rdmpi2 : the rotation matrix for 3 and 4.
c    rdsq3, rdmsq3 : the rotation matrix for 1 and 2.
c
c  workspace :
c    carray : the square root of the binomial numbers.
c      these numbers are only used here in this subroutine.
c
c  subroutine called :
c    bnlcft(), fstrtn(), datan, dacos,
c
c  called from : main()
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
c
      IMPLICIT NONE
c
      INTEGER *4 NTERMS
      REAL *8 CARRAY(0:2*NTERMS,0:2*NTERMS)
      REAL *8 RDPI2(0:NTERMS,0:NTERMS,-NTERMS:NTERMS)
      REAL *8 RDMPI2(0:NTERMS,0:NTERMS,-NTERMS:NTERMS)
      REAL *8 RDSQ3(0:NTERMS,0:NTERMS,-NTERMS:NTERMS)
      REAL *8 RDMSQ3(0:NTERMS,0:NTERMS,-NTERMS:NTERMS)
      REAL *8 DC(0:2*NTERMS,0:2*NTERMS)
c
c-----local varialbles.
c
      REAL *8 PI,THETA
c
c-----call initialization routines
c
      PI = 4.0D0*DATAN(1.0D0)
c
      CALL BNLCFT(DC,2*NTERMS)
c
      THETA = PI/2.0D0
      CALL FSTRTN(NTERMS,RDPI2,DC,THETA)
      THETA = -PI/2.0D0
      CALL FSTRTN(NTERMS,RDMPI2,DC,THETA)
      THETA = DACOS(DSQRT(3.0D0)/3.0D0)
      CALL FSTRTN(NTERMS,RDSQ3,DC,THETA)
      THETA = DACOS(-DSQRT(3.0D0)/3.0D0)
      CALL FSTRTN(NTERMS,RDMSQ3,DC,THETA)
c
      RETURN
      END
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
      SUBROUTINE ADJ_MEM_PTR(MEM_PTR,ASSIGN_PTR,SIZE)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
c
c  purpose:
c    computes the total required memory and the pointer for the
c    memory chunk.
c
c  on input:
c    mem_ptr: the previous pointer to the end of the memory chunk.
c    size: the size of the memory to be added.
c
c  on output:
c    assign_ptr: the pointer to the new memory request.
c    mem_ptr: the updated total memory requirement.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
c
      IMPLICIT NONE
c
      INTEGER*4 MEM_PTR,ASSIGN_PTR,SIZE
c
      ASSIGN_PTR = MEM_PTR
      MEM_PTR = MEM_PTR + SIZE
c
      RETURN
      END
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
      SUBROUTINE ICOPY(N,DX,DY)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
c
c  purpose:
c    copies an integer vector x of length n to another vector y of
c      the same length.
c
c  on input:
c    n: the number of elements to be copied.
c    dx: the vector to be copied.
c
c  on output:
c    dy: the duplicated vector.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
c
      IMPLICIT NONE
c
      INTEGER *4 N
      INTEGER *4 DX(N),DY(N)
      INTEGER I
c
      IF(N.LE.0)RETURN
c
      DO 10 I = 1,N
        DY(I) = DX(I)
10    CONTINUE
c
      RETURN
      END

