      SUBROUTINE PRINI(IP1,IQ1)
      CHARACTER *1 MES(1), AA(1)
      REAL *4 A(1)
      REAL *8 A2(1)
      INTEGER *4 IA(1)
      INTEGER *2 IA2(1)
      SAVE
      IP=IP1
      IQ=IQ1
      RETURN
c
c-----print real numbers.
c
      ENTRY PRIN2(MES,A2,N)
      CALL MESSPR(MES,IP,IQ)
      IF(IP.NE.0 .AND. N.NE.0) WRITE(IP,1400)(A2(J),J=1,N)
      IF(IQ.NE.0 .AND. N.NE.0) WRITE(IQ,1400)(A2(J),J=1,N)
 1400 FORMAT(6(2X,E11.5))
      RETURN
c
c-----print integer vectors.
c
      ENTRY PRINF(MES,IA,N)
      CALL MESSPR(MES,IP,IQ)
      IF(IP.NE.0 .AND. N.NE.0) WRITE(IP,1600)(IA(J),J=1,N)
      IF(IQ.NE.0 .AND. N.NE.0) WRITE(IQ,1600)(IA(J),J=1,N)
 1600 FORMAT(10(1X,I10))
      RETURN
      END
c
      SUBROUTINE MESSPR(MES,IP,IQ)
      CHARACTER *1 MES(1),AST
      DATA AST/'*'/
c
c-----determine the length of the message
c
      I=0
      DO 1400 I=1,10000
        IF(MES(I).EQ.AST) GOTO 1600
        I1=I
 1400 CONTINUE
 1600 CONTINUE
      IF ( (I1.NE.0) .AND. (IP.NE.0) )
     1  WRITE(IP,1800) (MES(I),I=1,I1)
      IF ( (I1.NE.0) .AND. (IQ.NE.0) )
     1  WRITE(IQ,1800) (MES(I),I=1,I1)
 1800 FORMAT(1X,80A1)
c
      RETURN
      END
c

