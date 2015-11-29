cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
      SUBROUTINE MKUPLIST(MYBOX,BOX,CENTER,SIZE,IWORK,
     1  IUALL,NUALL,IXUALL,IYUALL,IU1234,NU1234,IX1234,IY1234)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
c
c  purpose:
c
c     compute iuall and iu1234 interaction lists.
c
c  on input:
c
c     mybox: current box number.
c     box(16): box # mybox information.
c     center(3): center of current box.
c     size: size of current box.
c     iwork: array contains the box and list information.
c
c  on output:
c
c     iuall(*) (integer)
c		boxes at child level receiving up expansion from all
c               eight source box children.
c     iu1234(*) (integer)
c		boxes at child level receiving up expansion from lower
c               set of source box children (i.e. 1,2,3,4).
c     nuall     (integer)
c               number of boxes in iuall list.
c     nu1234    (integer)
c               number of boxes in iu1234 lists.
c     ixuall(*)  (integer)
c               integer x offset of target boxes in iuall.
c     iyuall(*)  (integer)
c               integer y offset of target boxes in iuall.
c     ix1234(*)  (integer)
c               integer x offset of target boxes in iu1234.
c     iy1234(*)  (integer)
c               integer y offset of target boxes in iu1234.
c
c  note: it is still not clear to me whether merging is a good idea.
c        for current implementation, the interaction list is generated
c        each time the code is executed. maybe we can precompute everything
c        in the future, however, it requires more storage.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
c
      IMPLICIT NONE
c
      INTEGER *4 MYBOX,BOX(15),IWORK(1)
      INTEGER *4 IUALL(36),NUALL,IXUALL(36),IYUALL(36)
      INTEGER *4 IU1234(16),NU1234,IX1234(16),IY1234(16)
      REAL *8 CENTER(3),SIZE
c
c-----local variables:
c
      INTEGER *4 IER,DADCOLLS(26),NCOLLS,LUSED
      INTEGER *4 NKIDS,I,J,ICOLL,KID
      INTEGER *4 BOXP(15),BOXC(15),NKIDSP,NKIDSC
      INTEGER *4 IX,IY
      REAL *8 CENTERP(3),CENTERC(3),SIZEP,SIZEC
c
c-----functions called:
c
      INTEGER *4 NINT
c
c-----get current box's colleagues, list #5.
c
      CALL D3MGETLIST(IER,MYBOX,5,DADCOLLS,NCOLLS,IWORK)
c
c-----find the children of the daddy's collegues:
c     note: we are currently not using list 2, but creates
c           the lists directly using parent information.
c
      NUALL=0
      NU1234=0
c
      DO 1600 I = 1,NCOLLS
        ICOLL = DADCOLLS(I)
c
c-------get parent box information.
c
        CALL D3MGETB(IER,ICOLL,BOXP,NKIDSP,CENTERP,SIZEP,IWORK)
c
c-------if childless, skip.
c
        IF (NKIDSP.LE.0) GOTO 1600
c
c-------if in the wrong direction, skip.
c
        IF (CENTERP(3)-CENTER(3) .LE. SIZEP/2.0D0) GOTO 1600
c
c-------right direction, check its children.
c
        DO 1400 J = 6,13
          KID = BOXP(J)
          IF (KID .GT. 0) THEN
c
c-----------get child box information.
c
            CALL D3MGETB(IER,KID,BOXC,NKIDSC,CENTERC,SIZEC,IWORK)
c
c-----------check if current box is inside the region.
c
            IX=NINT( (CENTERC(1)-CENTER(1)+SIZEC/2.0D0)/SIZEC )
            IY=NINT( (CENTERC(2)-CENTER(2)+SIZEC/2.0D0)/SIZEC )
c
            IF (CENTERC(3)-CENTER(3) .GT. 2.0D0*SIZEC) THEN
c
c-------------current child is in upall.
c
              NUALL=NUALL+1
              IUALL(NUALL)=KID
              IXUALL(NUALL)=IX
              IYUALL(NUALL)=IY
            ELSE
c
c-------------check if current box should be added to up1234.
c
              IF (IX.EQ.-2 .OR. IX.EQ.3 .OR. IY.EQ.-2 .OR.
     1          IY.EQ.3) GOTO 1400
              NU1234=NU1234+1
              IU1234(NU1234)=KID
              IX1234(NU1234)=IX
              IY1234(NU1234)=IY
            ENDIF
          ENDIF
1400    CONTINUE
1600  CONTINUE
c
      RETURN
      END
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
      SUBROUTINE MKDNLIST(MYBOX,BOX,CENTER,SIZE,IWORK,IDALL,NDALL,
     1  IXDALL,IYDALL,ID5678,ND5678,IX5678,IY5678)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
c
c  purpose:
c
c     compute idall and id5678 interaction lists.
c
c  on input:
c
c     mybox: current box number.
c     box(16): box # mybox information.
c     center(3): center of current box.
c     size: size of current box.
c     iwork: array contains the box and list information.
c
c  on output:
c
c     idall(*) (integer)
c		boxes at child level receiving down expansion from all
c               eight source box children.
c     id5678(*) (integer)
c		boxes at child level receiving down expansion from
c               upper set of source box children (i.e. 5,6,7,8).
c     ndall    (integer)
c               number of boxes in idnall list.
c     nd5678    (integer)
c               number of boxes in id5678 list.
c     ixdall(*)  (integer)
c               integer x offset of target boxes in idnall.
c     iydall(*)  (integer)
c               integer y offset of target boxes in idnall.
c     ix5678(*)  (integer)
c               integer x offset of target boxes in id5678.
c     iy5678(*)  (integer)
c               integer y offset of target boxes in id5678.
c
c  function called: nint().
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
c
      IMPLICIT NONE
c
      INTEGER *4 MYBOX,BOX(15),IWORK(1)
      INTEGER *4 IDALL(36),NDALL,IXDALL(36),IYDALL(36)
      INTEGER *4 ID5678(16),ND5678,IX5678(16),IY5678(16)
      REAL *8 CENTER(3),SIZE
c
c-----local variables:
c
      INTEGER *4 IER,DADCOLLS(26),NCOLLS,LUSED
      INTEGER *4 NKIDS,I,J,ICOLL,KID
      INTEGER *4 BOXP(15),BOXC(15),NKIDSP,NKIDSC
      INTEGER *4 IX,IY
      REAL *8 CENTERP(3),CENTERC(3),SIZEP,SIZEC
c
c-----functions called.
c
      INTEGER *4 NINT
c
c-----get current box's colleagues, list #5.
c
      CALL D3MGETLIST(IER,MYBOX,5,DADCOLLS,NCOLLS,IWORK)
c
c-----find the children of the daddy's collegues:
c     note: we are currently not using list 2, but creates
c           the lists directly using parent information.
c
      NDALL = 0
      ND5678 = 0
c
      DO 1600 I = 1,NCOLLS
        ICOLL = DADCOLLS(I)
c
c-------get parent box information.
c
        CALL D3MGETB(IER,ICOLL,BOXP,NKIDSP,CENTERP,SIZEP,IWORK)
c
c-------if childless, skip.
c
        IF (NKIDSP.LE.0) GOTO 1600
c
c-------if in the wrong direction, skip.
c
        IF (CENTERP(3)-CENTER(3) .GE. -SIZEP/2.0D0) GOTO 1600
c
c-------right direction, check its children.
c
        DO 1400 J = 6,13
          KID = BOXP(J)
          IF (KID .GT. 0) THEN
c
c-----------get child box information.
c
            CALL D3MGETB(IER,KID,BOXC,NKIDSC,CENTERC,SIZEC,IWORK)
c
c-----------check if current box is inside the region.
c
            IX=NINT( (CENTERC(1)-CENTER(1)+SIZEC/2.0D0)/SIZEC )
            IY=NINT( (CENTERC(2)-CENTER(2)+SIZEC/2.0D0)/SIZEC )
c
            IF (CENTERC(3)-CENTER(3) .LE. -2.0D0*SIZEC) THEN
c
c-------------current child is in udall.
c
              NDALL=NDALL+1
              IDALL(NDALL)=KID
              IXDALL(NDALL)=IX
              IYDALL(NDALL)=IY
            ELSE
c
c-------------check if current box should be added to down5678.
c
              IF (IX.EQ.-2 .OR. IX.EQ.3 .OR. IY.EQ.-2 .OR.
     1          IY.EQ.3) GOTO 1400
              ND5678=ND5678+1
              ID5678(ND5678)=KID
              IX5678(ND5678)=IX
              IY5678(ND5678)=IY
            ENDIF
          ENDIF
1400    CONTINUE
1600  CONTINUE
c
      RETURN
      END
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
      SUBROUTINE MKNOLIST(MYBOX,BOX,CENTER,SIZE,IWORK,
     1  INALL,NNALL,IXNALL,IYNALL,IN1256,NN1256,IX1256,IY1256,
     2  IN12,NN12,IX12,IY12,IN56,NN56,IX56,IY56)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
c
c  purpose:
c
c     compute north interaction lists.
c
c  on input:
c
c     mybox: current box number.
c     box(16): box # mybox information.
c     center(3): center of current box.
c     size: size of current box.
c     iwork: array contains the box and list information.
c
c  on output:
c
c     inall(*) (integer)
c		boxes at child level receiving north expansion from all
c               eight source box children.
c     nnall    (integer)
c               number of boxes in inall list.
c     ixnall(*) (integer)
c               integer x offset of target boxes in inall.
c     iynall(*) (integer)
c               integer y offset of target boxes in inall.
c
c     likewise for all other lists (see algorithm description)...
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
c
      IMPLICIT NONE
c
      INTEGER *4 MYBOX,BOX(15),IWORK(1)
      INTEGER *4 INALL(24),NNALL,IXNALL(24),IYNALL(24)
      INTEGER *4 IN1256(8),NN1256,IX1256(8),IY1256(8)
      INTEGER *4 IN12(4),NN12,IX12(4),IY12(4)
      INTEGER *4 IN56(4),NN56,IX56(4),IY56(4)
      REAL *8 CENTER(3),SIZE
c
c-----local variables.
c
      INTEGER *4 IER,DADCOLLS(26),NCOLLS,LUSED
      INTEGER *4 NKIDS,I,J,ICOLL,KID
      INTEGER *4 BOXP(15),BOXC(15),NKIDSP,NKIDSC
      INTEGER *4 IX,IY
      REAL *8 CENTERP(15),CENTERC(3),SIZEP,SIZEC
c
c-----functions called:
c
      INTEGER *4 NINT
c
c
c-----get current box's colleagues, list #5.
c
      CALL D3MGETLIST(IER,MYBOX,5,DADCOLLS,NCOLLS,IWORK)
c
c-----find the children of the daddy's collegues:
c     note: we are currently not using list 2, but creates
c           the lists directly using parent information.
c
      NNALL = 0
      NN1256 = 0
      NN12 = 0
      NN56 = 0
c
      DO 1600 I = 1,NCOLLS
        ICOLL = DADCOLLS(I)
c
c-------get parent box information.
c
        CALL D3MGETB(IER,ICOLL,BOXP,NKIDSP,CENTERP,SIZEP,IWORK)
c
c-------if childless, skip.
c
        IF (NKIDSP.LE.0) GOTO 1600
c
c-------if in the wrong direction, skip.
c
        IF (CENTERP(2)-CENTER(2) .LE. SIZEP/2.0D0) GOTO 1600
c
c-------right direction, check its children.
c
        DO 1400 J = 6,13
          KID = BOXP(J)
          IF (KID .GT. 0) THEN
c
c-----------get child box information.
c
            CALL D3MGETB(IER,KID,BOXC,NKIDSC,CENTERC,SIZEC,IWORK)
c
c-----------check if current box is inside the region.
c
            IX=NINT( (CENTERC(3)-CENTER(3)+SIZEC/2.0D0)/SIZEC )
            IY=NINT( (CENTERC(1)-CENTER(1)+SIZEC/2.0D0)/SIZEC )
c
            IF (CENTERC(2)-CENTER(2) .GT. 2.0D0*SIZEC) THEN
c
c-------------the upper level.
c
              IF ( IX.NE.-2 .AND. IX.NE.3) THEN
c
c---------------current child is in north-all. (otherwise it
c                 is in up-all and was processed)
c
                NNALL=NNALL+1
                INALL(NNALL)=KID
                IXNALL(NNALL)=IX
                IYNALL(NNALL)=IY
              ENDIF
            ELSE
c
c-------------it is in lower level.
c
              IF ((IX.EQ.0 .OR. IX.EQ.1) .AND. IY.GE.-1 .AND. IY.LE.2)
     1          THEN
c
c---------------current box should be added to north1256.
c
                NN1256=NN1256+1
                IN1256(NN1256)=KID
                IX1256(NN1256)=IX
                IY1256(NN1256)=IY
              ELSEIF (IX.EQ.-1. .AND. (IY.GE.-1 .AND. IY.LE.2)) THEN
c
c---------------current box should be added to north12
c
                NN12=NN12+1
                IN12(NN12)=KID
                IX12(NN12)=IX
                IY12(NN12)=IY
              ELSEIF (IX.EQ.2. .AND. (IY.GE.-1 .AND. IY.LE.2)) THEN
c
c---------------current box should be added to north56
c
                NN56=NN56+1
                IN56(NN56)=KID
                IX56(NN56)=IX
                IY56(NN56)=IY
              ENDIF
            ENDIF
          ENDIF
1400    CONTINUE
1600  CONTINUE
c
      RETURN
      END
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
      SUBROUTINE MKSOLIST(MYBOX,BOX,CENTER,SIZE,IWORK,
     1  ISALL,NSALL,IXSALL,IYSALL,IS3478,NS3478,IX3478,IY3478,
     2  IS34,NS34,IX34,IY34,IS78,NS78,IX78,IY78)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
c
c  purpose:
c
c     compute south interaction lists.
c
c  on input:
c
c     mybox: current box number.
c     box(16): box # mybox information.
c     center(3): center of current box.
c     size: size of current box.
c     iwork: array contains the box and list information.
c
c  on output:
c
c     isall(*) (integer)
c		boxes at child level receiving south expansion from all
c               eight source box children.
c     nsall    (integer)
c               number of boxes in isall list.
c     ixsall(*) (integer)
c               integer x offset of target boxes in isall.
c     iysall(*) (integer)
c               integer y offset of target boxes in isall.
c
c     likewise for all other lists (see algorithm description)...
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
c
      IMPLICIT NONE
c
      INTEGER *4 MYBOX,BOX(15),IWORK(1)
      INTEGER *4 ISALL(24),NSALL,IXSALL(24),IYSALL(24)
      INTEGER *4 IS3478(8),NS3478,IX3478(8),IY3478(8)
      INTEGER *4 IS34(4),NS34,IX34(4),IY34(4)
      INTEGER *4 IS78(4),NS78,IX78(4),IY78(4)
      REAL *8 CENTER(3),SIZE
c
c-----local variables:
c
      INTEGER *4 IER,DADCOLLS(26),NCOLLS,LUSED
      INTEGER *4 NKIDS,I,J,ICOLL,KID
      INTEGER *4 BOXP(15),BOXC(15),NKIDSP,NKIDSC
      INTEGER *4 IX,IY
      REAL *8 CENTERP(3),CENTERC(3),SIZEP,SIZEC
c
c-----functions called:
c
      INTEGER *4 NINT
c
c-----get current box's colleagues, list #5.
c
      CALL D3MGETLIST(IER,MYBOX,5,DADCOLLS,NCOLLS,IWORK)
c
c-----find the children of the daddy's collegues:
c     note: we are currently not using list 2, but creates
c           the lists directly using parent information.
c
      NSALL = 0
      NS3478 = 0
      NS34 = 0
      NS78 = 0
c
      DO 1600 I = 1,NCOLLS
        ICOLL = DADCOLLS(I)
c
c-------get parent box information.
c
        CALL D3MGETB(IER,ICOLL,BOXP,NKIDSP,CENTERP,SIZEP,IWORK)
c
c-------if childless, skip.
c
        IF (NKIDSP.LE.0) GOTO 1600
c
c-------if in the wrong direction, skip.
c
        IF (CENTERP(2)-CENTER(2) .GE. -SIZEP/2.0D0) GOTO 1600
c
c-------right direction, check its children.
c
        DO 1400 J = 6,13
          KID = BOXP(J)
          IF (KID .GT. 0) THEN
c
c-----------get child box information.
c
            CALL D3MGETB(IER,KID,BOXC,NKIDSC,CENTERC,SIZEC,IWORK)
c
c-----------check if current box is inside the region.
c
            IX=NINT( (CENTERC(3)-CENTER(3)+SIZEC/2.0D0)/SIZEC )
            IY=NINT( (CENTERC(1)-CENTER(1)+SIZEC/2.0D0)/SIZEC )
c
            IF (CENTERC(2)-CENTER(2) .LE. -2.0D0*SIZEC) THEN
c
c-------------the lowest level.
c
              IF ( IX.NE.-2 .AND. IX.NE.3) THEN
c
c---------------current child is in south-all. (otherwise it
c                 is in down-all and was processed)
c
                NSALL=NSALL+1
                ISALL(NSALL)=KID
                IXSALL(NSALL)=IX
                IYSALL(NSALL)=IY
              ENDIF
            ELSE
c
c-------------it is in lower level.
c
              IF ((IX.EQ.0 .OR. IX.EQ.1) .AND. IY.GE.-1 .AND. IY.LE.2)
     1          THEN
c
c---------------current box should be added to south3478.
c
                NS3478=NS3478+1
                IS3478(NS3478)=KID
                IX3478(NS3478)=IX
                IY3478(NS3478)=IY
              ELSEIF (IX.EQ.-1. .AND. (IY.GE.-1 .AND. IY.LE.2)) THEN
c
c---------------current box should be added to south34
c
                NS34=NS34+1
                IS34(NS34)=KID
                IX34(NS34)=IX
                IY34(NS34)=IY
              ELSEIF (IX.EQ.2. .AND. (IY.GE.-1 .AND. IY.LE.2)) THEN
c
c---------------current box should be added to south78
c
                NS78=NS78+1
                IS78(NS78)=KID
                IX78(NS78)=IX
                IY78(NS78)=IY
              ENDIF
            ENDIF
          ENDIF
1400    CONTINUE
1600  CONTINUE
c
      RETURN
      END
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
      SUBROUTINE MKEALIST(MYBOX,BOX,CENTER,SIZE,IWORK,
     1  IEALL,NEALL,IXEALL,IYEALL,IE1357,NE1357,IX1357,IY1357,
     2  IE13,NE13,IX13,IY13,IE57,NE57,IX57,IY57,IE1,NE1,IX1,IY1,
     3  IE3,NE3,IX3,IY3,IE5,NE5,IX5,IY5,IE7,NE7,IX7,IY7)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
c
c  purpose:
c
c     compute east interaction lists.
c
c  on input:
c
c     mybox: current box number.
c     box(16): box # mybox information.
c     center(3): center of current box.
c     size: size of current box.
c     iwork: array contains the box and list information.
c
c  on output:
c
c     ieall(*) (integer)
c		boxes at child level receiving east expansion from all
c               eight source box children.
c     neall    (integer)
c               number of boxes in ieall list.
c     ixeall(*) (integer)
c               integer x offset of target boxes in ieall list.
c     iyeall(*) (integer)
c               integer y offset of target boxes in ieall list.
c
c     likewise for all other lists (see algorithm description)...
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
c
      IMPLICIT NONE
c
      INTEGER *4 MYBOX,BOX(15),IWORK(1)
      INTEGER *4 IEALL(16),NEALL,IXEALL(16),IYEALL(16)
      INTEGER *4 IE1357(4),NE1357,IX1357(4),IY1357(4)
      INTEGER *4 IE13(2),NE13,IX13(2),IY13(2)
      INTEGER *4 IE57(2),NE57,IX57(2),IY57(2)
      INTEGER *4 IE1(3),NE1,IX1(3),IY1(3)
      INTEGER *4 IE3(3),NE3,IX3(3),IY3(3)
      INTEGER *4 IE5(3),NE5,IX5(3),IY5(3)
      INTEGER *4 IE7(3),NE7,IX7(3),IY7(3)
      REAL *8 CENTER(3),SIZE
c
c-----local variables.
c
      INTEGER *4 IER,DADCOLLS(26),NCOLLS,LUSED
      INTEGER *4 NKIDS,I,J,ICOLL,KID
      INTEGER *4 BOXP(15),BOXC(15),NKIDSP,NKIDSC
      INTEGER *4 IX,IY
      REAL *8 CENTERP(3),CENTERC(3),SIZEP,SIZEC
c
c-----functions called:
c
      INTEGER *4 NINT
c
c
c-----get current box's colleagues, list #5.
c
      CALL D3MGETLIST(IER,MYBOX,5,DADCOLLS,NCOLLS,IWORK)
c
c-----find the children of the daddy's collegues:
c     note: we are currently not using list 2, but creates
c           the lists directly using parent information.
c
      NEALL = 0
      NE1357 = 0
      NE13 = 0
      NE57 = 0
      NE1 = 0
      NE3 = 0
      NE5 = 0
      NE7 = 0
c
      DO 1600 I = 1,NCOLLS
        ICOLL = DADCOLLS(I)
c
c-------get parent box information.
c
        CALL D3MGETB(IER,ICOLL,BOXP,NKIDSP,CENTERP,SIZEP,IWORK)
c
c-------if childless, skip.
c
        IF (NKIDSP.LE.0) GOTO 1600
c
c-------if in the wrong direction, skip.
c
        IF (CENTERP(1)-CENTER(1) .LE. SIZEP/2.0D0) GOTO 1600
c
c-------right direction, check its children.
c
        DO 1400 J = 6,13
          KID = BOXP(J)
          IF (KID .GT. 0) THEN
c
c-----------get child box information.
c
            CALL D3MGETB(IER,KID,BOXC,NKIDSC,CENTERC,SIZEC,IWORK)
c
c-----------check if current box is inside the region.
c
            IX=NINT( (CENTERC(3)-CENTER(3)+SIZEC/2.0D0)/SIZEC )
            IY=NINT( (CENTERC(2)-CENTER(2)+SIZEC/2.0D0)/SIZEC )
c
            IF (CENTERC(1)-CENTER(1) .GT. 2.0D0*SIZEC) THEN
c
c-------------the upper level.
c
              IF ( IX.GE.-1 .AND. IX.LE.2 .AND. IY.GE.-1 .AND. IY.LE.2)
     1          THEN
c
c---------------current child is in east-all. (otherwise it
c                 is in ud-all or ns-all and was processed)
c
                NEALL=NEALL+1
                IEALL(NEALL)=KID
                IXEALL(NEALL)=-IX
                IYEALL(NEALL)=IY
              ENDIF
            ELSE
c
c-------------it is in lower level.
c
              IF ((IX.EQ.0 .OR. IX.EQ.1) .AND. (IY.EQ.0 .OR. IY.EQ.1))
     1          THEN
c
c---------------current box should be added to east1357.
c
                NE1357=NE1357+1
                IE1357(NE1357)=KID
                IX1357(NE1357)=-IX
                IY1357(NE1357)=IY
              ELSEIF (IX.EQ.-1 .AND. (IY.EQ.0 .OR. IY.EQ.1)) THEN
c
c---------------current box should be added to east13
c
                NE13=NE13+1
                IE13(NE13)=KID
                IX13(NE13)=-IX
                IY13(NE13)=IY
              ELSEIF (IX.EQ.2 .AND. (IY.EQ.0 .OR. IY.EQ.1)) THEN
c
c---------------current box should be added to east57
c
                NE57=NE57+1
                IE57(NE57)=KID
                IX57(NE57)=-IX
                IY57(NE57)=IY
              ELSEIF (IY.EQ.-1) THEN
                IF (IX.GE.-1 .AND. IX.LE.1) THEN
c
c-----------------current box should be added to east1
c
                  NE1=NE1+1
                  IE1(NE1)=KID
                  IX1(NE1)=-IX
                  IY1(NE1)=IY
                ENDIF
c
                IF (IX.GE.0 .AND. IX.LE.2) THEN
c
c-----------------current box should be added to east5
c
                  NE5=NE5+1
                  IE5(NE5)=KID
                  IX5(NE5)=-IX
                  IY5(NE5)=IY
                ENDIF
              ELSEIF (IY.EQ.2) THEN
                IF (IX.GE.-1 .AND. IX.LE.1) THEN
c
c-----------------current box should be added to east3
c
                  NE3=NE3+1
                  IE3(NE3)=KID
                  IX3(NE3)=-IX
                  IY3(NE3)=IY
                ENDIF
                IF (IX.GE.0 .AND. IX.LE.2) THEN
c
c-----------------current box should be added to east7
c
                  NE7=NE7+1
                  IE7(NE7)=KID
                  IX7(NE7)=-IX
                  IY7(NE7)=IY
                ENDIF
              ENDIF
            ENDIF
          ENDIF
1400    CONTINUE
1600  CONTINUE
c
      RETURN
      END
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
      SUBROUTINE MKWELIST(MYBOX,BOX,CENTER,SIZE,IWORK,
     1  IWALL,NWALL,IXWALL,IYWALL,IW2468,NW2468,IX2468,IY2468,IW24,
     2  NW24,IX24,IY24,IW68,NW68,IX68,IY68,IW2,NW2,IX2,IY2,
     3  IW4,NW4,IX4,IY4,IW6,NW6,IX6,IY6,IW8,NW8,IX8,IY8)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
c
c  purpose:
c
c     compute west interaction lists.
c
c  on input:
c
c     mybox: current box number.
c     box(16): box # mybox information.
c     center(3): center of current box.
c     size: size of current box.
c     iwork: array contains the box and list information.
c
c  on output:
c
c     iwall(*) (integer)
c		boxes at child level receiving east expansion from all
c               eight source box children.
c     nwall    (integer)
c               number of boxes in iwall list.
c     ixwall(*) (integer)
c               integer x offset of target boxes in iwall list.
c     iywall(*) (integer)
c               integer y offset of target boxes in iwall list.
c
c     likewise for all other lists (see algorithm description)...
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      IMPLICIT NONE
c
      INTEGER *4 MYBOX,BOX(15),IWORK(1)
      INTEGER *4 IWALL(16),NWALL,IXWALL(16),IYWALL(16)
      INTEGER *4 IW2468(4),NW2468,IX2468(4),IY2468(4)
      INTEGER *4 IW24(2),NW24,IX24(2),IY24(2)
      INTEGER *4 IW68(2),NW68,IX68(2),IY68(2)
      INTEGER *4 IW2(3),NW2,IX2(3),IY2(3)
      INTEGER *4 IW4(3),NW4,IX4(3),IY4(3)
      INTEGER *4 IW6(3),NW6,IX6(3),IY6(3)
      INTEGER *4 IW8(3),NW8,IX8(3),IY8(3)
      REAL *8 CENTER(3),SIZE
c
c-----local variables.
c
      INTEGER *4 IER,DADCOLLS(26),NCOLLS,LUSED
      INTEGER *4 NKIDS,I,J,ICOLL,KID
      INTEGER *4 BOXP(15),BOXC(15),NKIDSP,NKIDSC
      INTEGER *4 IX,IY
      REAL *8 CENTERP(3),CENTERC(3),SIZEP,SIZEC
c
c-----functions called:
c
      INTEGER *4 NINT
c
c-----get current box's colleagues, list #5.
c
      CALL D3MGETLIST(IER,MYBOX,5,DADCOLLS,NCOLLS,IWORK)
c
c-----find the children of the daddy's collegues:
c     note: we are currently not using list 2, but creates
c           the lists directly using parent information.
c
      NWALL = 0
      NW2468 = 0
      NW24 = 0
      NW68 = 0
      NW2 = 0
      NW4 = 0
      NW6 = 0
      NW8 = 0
c
      DO 1600 I = 1,NCOLLS
        ICOLL = DADCOLLS(I)
c
c-------get parent box information.
c
        CALL D3MGETB(IER,ICOLL,BOXP,NKIDSP,CENTERP,SIZEP,IWORK)
c
c-------if childless, skip.
c
        IF (NKIDSP.LE.0) GOTO 1600
c
c-------if in the wrong direction, skip.
c
        IF (CENTERP(1)-CENTER(1) .GT. -SIZEP/2.0D0) GOTO 1600
c
c-------right direction, check its children.
c
        DO 1400 J = 6,13
          KID = BOXP(J)
          IF (KID .GT. 0) THEN
c
c-----------get child box information.
c
            CALL D3MGETB(IER,KID,BOXC,NKIDSC,CENTERC,SIZEC,IWORK)
c
c-----------check if current box is inside the region.
c
            IX=NINT( (CENTERC(3)-CENTER(3)+SIZEC/2.0D0)/SIZEC )
            IY=NINT( (CENTERC(2)-CENTER(2)+SIZEC/2.0D0)/SIZEC )
c
            IF (CENTERC(1)-CENTER(1) .LT. -2.0D0*SIZEC) THEN
c
c-------------the lowest level.
c
              IF ( IX.GE.-1 .AND. IX.LE.2 .AND. IY.GE.-1 .AND. IY.LE.2)
     1          THEN
c
c---------------current child is in east-all. (otherwise it
c                 is in ud-all or ns-all and was processed)
c
                NWALL=NWALL+1
                IWALL(NWALL)=KID
                IXWALL(NWALL)=-IX
                IYWALL(NWALL)=IY
              ENDIF
            ELSE
c
c-------------it is in lower level.
c
              IF ((IX.EQ.0 .OR. IX.EQ.1) .AND. (IY.EQ.0 .OR. IY.EQ.1))
     1          THEN
c
c---------------current box should be added to west2468.
c
                NW2468=NW2468+1
                IW2468(NW2468)=KID
                IX2468(NW2468)=-IX
                IY2468(NW2468)=IY
              ELSEIF (IX.EQ.-1 .AND. (IY.EQ.0 .OR. IY.EQ.1)) THEN
c
c---------------current box should be added to west24
c
                NW24=NW24+1
                IW24(NW24)=KID
                IX24(NW24)=-IX
                IY24(NW24)=IY
              ELSEIF (IX.EQ.2 .AND. (IY.EQ.0 .OR. IY.EQ.1)) THEN
c
c---------------current box should be added to west68
c
                NW68=NW68+1
                IW68(NW68)=KID
                IX68(NW68)=-IX
                IY68(NW68)=IY
              ELSEIF (IY.EQ.-1) THEN
                IF (IX.GE.-1 .AND. IX.LE.1) THEN
c
c-----------------current box should be added to west2
c
                  NW2=NW2+1
                  IW2(NW2)=KID
                  IX2(NW2)=-IX
                  IY2(NW2)=IY
                ENDIF
c
                IF (IX.GE.0 .AND. IX.LE.2) THEN
c
c-----------------current box should be added to west6
c
                  NW6=NW6+1
                  IW6(NW6)=KID
                  IX6(NW6)=-IX
                  IY6(NW6)=IY
                ENDIF
              ELSEIF (IY.EQ.2) THEN
                IF ( IX.GE.-1 .AND. IX.LE.1) THEN
c
c-----------------current box should be added to w4
c
                  NW4=NW4+1
                  IW4(NW4)=KID
                  IX4(NW4)=-IX
                  IY4(NW4)=IY
                ENDIF
c
                IF (IX.GE.0 .AND. IX.LE.2) THEN
c
c-----------------current box should be added to west8
c
                  NW8=NW8+1
                  IW8(NW8)=KID
                  IX8(NW8)=-IX
                  IY8(NW8)=IY
                ENDIF
              ENDIF
            ENDIF
          ENDIF
1400    CONTINUE
1600  CONTINUE
c
      RETURN
      END
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
c
c  this is the 3-d adaptive box data structure code
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
c
c  important notes:
c
c    1) in three-dimensional space, chidren are ordered as follows:
c
c                3 4  |  7 8
c                1 2  |  5 6
c
c    2) the list of chidren may contain holes
c
c       i.e. if 6,7,8,9,10,11,12,13 - the list of children of the
c         box ibox (eight of them, and the child is identified by
c         its address in the array boxes), then the 6,7,8,9-th
c         elements may be zero and 10,11,12,13-th elements may be
c         non-zero.
c
c    3) corners array currently stores the sides of the boxes only
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
        SUBROUTINE D3MSTRCR(IER,Z,N,NBOX,NBOXES,IZ,LADDR,NLEV,
     1    CENTER,SIZE,IWORK,LW,LUSED777,EPSCLOSE,NINIRE)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
c
c purpose:
c
c   this subroutine constructs the logical structure for the
c     fully adaptive fmm in three dimensions and stores it in the
c     array w in the form of a link-list. including, a quad-tree
c     of boxes formed, an array of has_kid indicaters, an arrary
c     of centers for all boxes, and all (5 among them 4 will be used
c     by fmm) lists.
c
c   the user can access the information about various boxes and
c     lists in it by calling the entries d3mgetb, d3mgetl, d3mlinfo
c     of this subroutine (see).
c
c on input:
c
c   z(3,n): the user-specified points in the space
c   n: the number of elements in array z
c   nbox: the maximum number of points allowed in a childless box
c   lw: the amount of memory in the array w (in integer*4 elements)
c   epsclose: the smallest size a box can be (the length of sides).
c   ninire: the ratio for storing real*8 and integer *4. this is
c           required as only iwork is used for storage.
c           this can be avoided later if a separate workspace is introduced.
c           note: for ifort compiler:
c              real *8 = 8 bytes. integer *4 = 4 bytes. so ninire=2.
c
c on output:
c
c   ier: error return code
c        ier=0   means successful execution
c        ier=16 means that the subroutine attempted to construct more
c          than 199 levels of subdivision; the code will pause
c          and one should stop to check, for any physically
c          meaningful distribution of points, it should not go
c          higher than that.
c        ier=32  means that the amount lw of space in array w
c          is insufficient
c   nboxes: the total number of boxes created
c   iz(n): the integer array addressing the particles in
c          all boxes. explanation: for a box ibox, the particles living
c          in it are:
c          (z(1,j),z(2,j),z(3,j)),(z(1,j+1),z(2,j+1),z(3,j+1)),
c          (z(1,j+2),z(2,j+2),z(3,j+3)), . . .
c          (z(1,j+nj-1),z(2,j+nj-1),z(3,j+nj-1)),
c          (z(1,j+nj),z(2,j+nj),z(3,j+nj)),
c          with j=boxes(14,ibox), and nj=boxes(15,ibox)-1
c   laddr(2,200): an integer array (2,nlev), describing the
c          numbers of boxes on various levels of subdivision, so that
c          the first box on level (i-1) has sequence number laddr(1,i),
c          and there are laddr(2,i) boxes on level i-1
c   nlev: the maximum level number on which any boxes have
c         been created. the maximim number possible is 200.
c         it is recommended that the array laddr above be dimensioned
c         at least (2,200), in case the user underestimates
c         the number of levels required.
c   center(3): the center of the box on the level 0, containing
c         the whole simulation
c   size: the side of the box on the level 0
c   iwork(lw): the array containing all tables describing boxes,
c         lists, etc. it is a link-list (for the most part), and can
c         only be accessed via the entries d3mgetb, d3mgetl, d3mlinfo,
c         of this subroutine (see below). the first lused777 integer*4
c         locations of this array should not be altered between the
c         call to this entry and subsequent calls to the entries
c         d3mgetb, d3mgetl, d3mlinfo, of this  subroutine
c   lused777: the amount of space in the array w (in integer*4
c            words) that is occupied by various tables on exit from this
c           subroutine. this space should not be altered between the
c           call to this entry and subsequent calls to entries d3mgetb,
c           d3mgetl, d3mlinfo, of this  subroutine (see below).
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
c
        IMPLICIT NONE
        SAVE
c
        INTEGER *4 IER,N,NBOX,NBOXES,NLEV,LW,LUSED777
        INTEGER *4 NINIRE
        INTEGER *4 IZ(N),IWORK(LW),LADDR(2,200)
        REAL *8 Z(3,N),CENTER(3),SIZE,EPSCLOSE
c
c-------local variables.
c
        INTEGER *4 I,IIWORK,LIWORK,IBOXES,LBOXES,MAXBOXES
        INTEGER *4 NN,IKIDS,LKIDS,ICENTERS,LCENTERS,ICORNERS,LCORNERS
        INTEGER *4 IWLISTS,LWLISTS,LUSED
        INTEGER *4 IIBOXES,IICENTERS,IICORNERS,IIWLISTS,IIKIDS
        DATA IIBOXES/1/,IICENTERS/2/,IICORNERS/3/,IIWLISTS/4/,IIKIDS/5/
c
c-------first construct the quad-tree structure for given set of points:
c
        IF( (LW .LT. 20*N) .OR. (LW .LT. 5000) ) THEN
c
c---------the array is simply too small, don't even try the solution
c
          IER=32
          RETURN
        ENDIF
c
c-------initialize the array iwork.
c
        DO I=1,LW
          IWORK(I)=0
        ENDDO
c
        IIWORK = 1
        LIWORK = N+4
c
        IBOXES = IIWORK+LIWORK
        LBOXES = LW-N-5
        MAXBOXES = LBOXES/15-1
c
c-------generate the octree structure.
c
        CALL D3MALLB(IER,Z,N,NBOX,IWORK(IBOXES),MAXBOXES,
     1    NBOXES,IZ,LADDR,NLEV,CENTER,SIZE,IWORK(IIWORK),EPSCLOSE)
c
c-------if the memory is insufficient - bomb
c
        IF (IER .EQ. 4) THEN
          CALL PRINF('fatal: not enough memory in d3mallb, ier=*',IER,1)
          IER = 32
          RETURN
        ENDIF
c
c-------compress iwork, so that the first 15*nboxes contain the information
c         in boxes derived in d3mallb(). allocating extra space for the
c         following pointers: iboxes, icenters, icorners, iwlists, ikids.
c
c-------1. shift the boxes information. iwork(6:6+15*nboxes)
c
        NN = NBOXES*15
        DO 1200 I = 1,NN
          IWORK(I + 5) = IWORK(IBOXES+I-1)
1200    CONTINUE
c
c-------2. iwork(iiboxes=1)=6. the starting position of box info in iwork.
c       3. iwork(iikids=5) = ikids. the starting position of kid information.
c          construct the array has_kid for all boxes: store it in the second
c            block of iw(starts at ikids, of length lkids), to indicate if
c            certain box is childless or not.
c       4. iwork(iicenters=2) = icenters.
c          construct the centers for all boxes in the octree. the
c            third block of iw (starts at icenters, of size lcenters) stores
c            the center coordinates for all boxes (note: we are storing
c            real *8 numbers in integer*4 storage space, not good, but for
c            simplicity)
c       5. iwork(iicorners=3) = icorners.
c       6. iwork(iiwlists=4) = iwlists.
c
        IBOXES = 1 + 5
        LBOXES = NBOXES*15+16
        IWORK(IIBOXES)=IBOXES
c
        IKIDS = IBOXES+LBOXES
        LKIDS = NBOXES+1
        IWORK(IIKIDS)=IKIDS
c
        ICENTERS = IKIDS+LKIDS
        LCENTERS = (NBOXES*3+2)*NINIRE
        IWORK(IICENTERS)=ICENTERS
c
        ICORNERS = ICENTERS+LCENTERS
        LCORNERS = (NBOXES+2)*NINIRE
        IWORK(IICORNERS)=ICORNERS
c
        IWLISTS = ICORNERS+LCORNERS
        LWLISTS = LW-IWLISTS-6
        IWORK(IIWLISTS)=IWLISTS
c
c-------create an array to indicate if the box has kids.
c
        CALL D3MHASKID(IWORK(IWORK(IIBOXES)),
     1     NBOXES,IWORK(IWORK(IIKIDS)))
c
c-------generate the centers and the side length of all boxes.
c
        CALL D3MCENTC(CENTER,SIZE,IWORK(IWORK(IIBOXES)),NBOXES,
     1     IWORK(IWORK(IICENTERS)),IWORK(IWORK(IICORNERS)))
c
c-------now, construct all lists for all boxes
c
        CALL D3MCREALISTS(IER,IWORK(IWORK(IIBOXES)),NBOXES,SIZE,
     1    IWORK(IWORK(IICENTERS)),
     2    IWORK(IWORK(IIKIDS)),IWORK(IWORK(IIWLISTS)),LWLISTS,LUSED)
c
        LUSED777 = LUSED+IWLISTS
        IF (IER .EQ. 32)
     1  CALL PRINF('fatal: not enough memory in d3mcrealists, ier=*',
     2     IER,1)
c
        RETURN
        END
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
        SUBROUTINE D3MCREALISTS(IER,BOXES,NBOXES,SIZE,CENTERS,HASKID,
     1     W,LW,LUSED)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
c
c purpose:
c
c   this subroutine constructs all lists for all boxes
c     and stores them in the storage area w in the form
c     of a link list. the resulting data can be accessed
c     by calls to various entries of the subroutine d3mlinkinit (see).
c
c on input:
c
c   boxes: an integer array dimensioned (15,nboxes), as created by
c           the subroutine d3mallb (see).
c   nboxes: the total number of boxes created
c   size: the side length of whole box
c   centers: the center of all boxes
c   haskid: an array of flag indicate if box has kid or not
c   lw: the total amount of storage in array w(in integer*4 words)
c
c output:
c
c   ier: the error return code.
c     ier=0 means successful execution
c     ier=32 means that the amount lw of space in array w is
c            insufficient. it is a fatal error.
c   w: storage area containing all lists for all boxes in
c      the form of link-lists, accessible by the subroutine
c      d3mlinkretr (see).
c   lused: the amount of space in the array w (in integer*4 words)
c          that is occupied by various tables on exit from this
c          subroutine. this space should not be altered between the
c          call to this subroutine and subsequent calls to the
c          entries d3mlinkretr, etc. of the subroutine d3mlinkinit (see).
c
c note on the list conventions:
c
c    list 1 of the box ibox - the list of all boxes with which the
c           box ibox interacts directly, including the boxes on the
c           same level as ibox, boxes on the finer levels, and boxes
c           on the coarser levels. obviously, list 1 is empty for any
c           box that is not childless. all boxes in list 1 must also
c           be childless.
c
c           note: direct interactions, two boxes are touching.
c
c           note: see cheng's adaptive 3d fmm paper, list 1.
c
c    list 2 of the box ibox - the list of all boxes with which the
c           box ibox interacts in the regular multipole fashion, i.e.
c           boxes on the same level as ibox that are separated from it
c           but whose daddies are not separated from the daddy of ibox.
c
c           note: this is the interaction list. the interactions are
c             computed using mpole->exp->local.
c
c    list 3 of the box ibox - for a childless ibox, the list of all
c           boxes on the levels finer than that of ibox, which are
c           separated from ibox, and whose daddys are not separated
c           from ibox. for a box with children, list 3 is empty.
c
c           note: see cheng's adaptive 3d fmm paper, list 3.
c
c           this list is computed using
c
c    list 4 is dual to list 3, i.e. jbox is on the list 4 of ibox if
c           and only if ibox is on the list 3 of jbox. all boxes in list
c            4 are of higher level and are childless.
c
c    list 5 is the colleages of box ibox - adjacent boxes on the same
c            level. also called neighbors.
c
c    list 6 of the box ibox - for a childless ibox, the list of all
c           boxes on the levels finer than that of ibox, which touch
c           the box ibox, and which have children. for a box jbox
c           with children, the list of all childless boxes on the
c           coarser levels, which touch the box jbox. obviously, this
c           is a generalization of list 1 for the boxes with children.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
c
        IMPLICIT NONE
c
        INTEGER *4 IER,NBOXES,LUSED,LW
        INTEGER *4 BOXES(15,1),W(1),HASKID(1)
        REAL *8 SIZE
        REAL *8 CENTERS(3,1)
c
        INTEGER *4 I,J,LUSED2,KID,IFINTER
        INTEGER *4 NTYPES,NTYPES4,ITYPE,ITYPE3,ITYPE4,ITYPE5,ITYPE2
        INTEGER *4 IBOX,JBOX,IDAD,NKIDS,NLIST,JER,NLIST1,ICOLL,NCOLLS
        INTEGER *4 COLLKIDS(5000),DADCOLLS(200),LIST5(2000),STACK(6000)
c
c-------initialize the storage-retrieval routine for all boxes:
c         also, set up correct lused2 to avoid the problem when
c         there is only one box and no lists.
c
        IER = 0
        NTYPES = 6
        LUSED = 0
        LUSED2 = 32+NTYPES*NBOXES+10
c
        CALL D3MLINKINIT(IER,NBOXES,NTYPES,W,LW)
c
        IF( IER .NE. 0 ) RETURN
c
c-------first for all boxes, construct - lists 5: colleagues of ibox,
c         and list 2: direct interaction using multipole expansion
c         note that list 2 becomes necessary for adaptive structures.
c         ===========================================================
c
        ITYPE5 = 5
        ITYPE2 = 2
c
        DO 2000 IBOX = 2,NBOXES
c
c---------find daddy of current box.
c
          IDAD = BOXES(5,IBOX)
c
c---------retrieve daddy's collegues, including daddy himself:
c           note that the box ibox itself is not stored in the colleague
c           list.
c
          DADCOLLS(1) = IDAD
c
          CALL D3MLINKRETR(IER,ITYPE5,IDAD,DADCOLLS(2),NCOLLS,W,LUSED)
c
          NCOLLS = NCOLLS+1
c
c---------find the children of the daddy's collegues:
c
          NKIDS = 0
          DO 1600 I = 1,NCOLLS
            ICOLL = DADCOLLS(I)
            IF (HASKID(ICOLL) .LE. 0) GOTO 1600
            DO 1400 J = 6,13
              KID = BOXES(J,ICOLL)
              IF (KID .GT. 0) THEN
                IF (KID .EQ. IBOX) GOTO 1400
                NKIDS = NKIDS+1
                COLLKIDS(NKIDS) = KID
              ENDIF
1400        CONTINUE
1600      CONTINUE
c
c---------sort the kids of the daddy's collegues into the
c           lists 2, 5 of the box ibox
c
          NLIST1 = 1
          DO 1800 I = 1,NKIDS
c
c-----------check if this kid is touching the box ibox, if touching, this kid
c             belongs to list 5 of box ibox, if not, it belongs to list 2
c             of box ibox:
c
            KID = COLLKIDS(I)
            CALL D3MINTERS(SIZE,CENTERS(1,IBOX),BOXES(1,IBOX),
     1        CENTERS(1,KID),BOXES(1,KID),IFINTER)
            IF (IFINTER .EQ. 1)
     1        CALL D3MLINKSTOR(IER,ITYPE5,IBOX,KID,NLIST1,W,LUSED)
            IF (IFINTER .EQ. 0)
     1        CALL D3MLINKSTOR(IER,ITYPE2,IBOX,KID,NLIST1,W,LUSED)
1800      CONTINUE
c
c---------possible error exiting, not enough memory space:
c
          IF (IER .EQ. 32) THEN
            RETURN
          ENDIF
c
2000    CONTINUE
c
c-------second, construct lists 1, 3:
c       ============================
c
        DO 3000 I = 1,NBOXES
c
c---------if this box has kids - its lists 1, 3 are empty, so skip it,
c           otherwise get its list 5, from which to construct lists 1,3:
c
          IF (HASKID(I) .GT. 0) GOTO 3000
c
          CALL D3MLINKRETR(JER,ITYPE5,I,LIST5,NLIST,W,LUSED)
c
c---------if list 5 is empty, skip, otherwise call the list 1, 3 maker:
c
          IF(JER .EQ. 4) GOTO 3000
          DO 3200 J = 1,NLIST
            JBOX = LIST5(J)
            CALL D3MLISTS31(IER,I,JBOX,BOXES,NBOXES,SIZE,CENTERS,
     1         HASKID,W,STACK,LUSED)
c
c-----------possible error exiting, not enough memory space:
c
            IF (IER .EQ. 32) THEN
              CALL PRINF('fatal: memory exhausted in building
     1          LISTS 1,3, IER=*',ier,1)
              RETURN
            ENDIF
3200      CONTINUE
3000    CONTINUE
c
c-------copy all elements of lists 1, 2, 3, 5, 6 while skipping list 4:
c       =========================================================
c
        NTYPES4 = 6
c
        CALL D3MLINKINIT(IER,NBOXES,NTYPES4,W(LUSED+1),LW-(LUSED+5))
        IF( IER .NE. 0 ) RETURN
        DO 3600 IBOX = 1,NBOXES
          DO 2400 ITYPE = 1,6
            CALL D3MLINKRETR(JER,ITYPE,IBOX,LIST5,NLIST,W,LUSED)
            IF (JER .EQ. 4) GOTO 2400
            CALL D3MLINKSTOR
     1          (IER,ITYPE,IBOX,LIST5,NLIST,W(LUSED+1),LUSED2)
            IF (IER .EQ. 32) THEN
               CALL PRINF('fatal: memory exhausted in realigning lists
     1                1,2,3, IER=*',ier,1)
               RETURN
            ENDIF
 2400     CONTINUE
 3600   CONTINUE
c
c-------compress array w:
c
        DO I = 1,LUSED2
          W(I) = W(LUSED+I)
        ENDDO
        LUSED = LUSED2
c
c-------finally, construct the lists 4 for all boxes that need them:
c       ===========================================================
c
        ITYPE3 = 3
        ITYPE4 = 4
        NLIST1 = 1
c
        DO 4000 IBOX = 1,NBOXES
c
          CALL D3MLINKRETR(JER,ITYPE3,IBOX,LIST5,NLIST,W,LUSED)
c
          IF (JER .EQ. 4) GOTO 4000
          DO J = 1,NLIST
            CALL D3MLINKSTOR(IER,ITYPE4,LIST5(J),IBOX,NLIST1,W,LUSED2)
            IF (IER .EQ. 32) THEN
              CALL PRINF('fatal: memory exhausted in building
     1          LIST 4, IER=*',ier,1)
              RETURN
            ENDIF
          ENDDO
4000    CONTINUE
        LUSED = LUSED2
c
        RETURN
        END
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
        SUBROUTINE D3MLISTS31(IER,IBOX,JBOX0,BOXES,NBOXES,SIZE,CENTERS,
     1    HASKID,W,STACK,LUSED)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
c
c purpose:
c
c   this subroutine constructs all elements of lists 1 and 3
c   resulting from the subdivision of one element of list 5
c   of the box ibox. all these elements of lists 1, 3 are
c   stored in the link-lists by the subroutine linstro (see)
c
c input:
c
c   ibox: the box whose lists are being constructed
c   jbox0: the element of list 5 of the box ibox that is being
c          subdivided
c   boxes: the array boxes as created by d3mallb
c   nboxes: the number of boxes in array boxes
c   size: the side length of whole box
c   centers: the array of centers of all the boxes
c   haskid: the array of indicaters of how many kid each box has
c   w: the storage area formatted by the subroutine d3mlinkinit (see)
c      to be used to store the elements of lists 1, 3 constructed
c      by this subroutine. obviously, by this time, it contains
c      plenty of other lists.
c
c output:
c
c   ier: the error return code.
c     =0 means successful execution
c     =32 means that the amount lw of space in array w is
c         insufficient. it is a fatal error.
c   w: the augmented storage area, containing all the boxes just
c      created, in addition to whatever had been stored previously
c   lused: the total length of array w (in integer words) used
c          on exit from this subroutine
c
c working arrays:
c          stack - must be at least 600 integer locations long
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        IMPLICIT NONE
c
        INTEGER *4 IER,IBOX,JBOX0,NBOXES,LUSED
        INTEGER *4 BOXES(15,1),HASKID(1)
        INTEGER *4 STACK(3,1)
        REAL *8 SIZE
        REAL *8 W(1),CENTERS(3,1)
c
        INTEGER *4 KK,NSONS,NNSONS,JBOX,ISTACK,IFINTER
        INTEGER *4 ITYPE1,ITYPE2,ITYPE3,ITYPE4,ITYPE5,ITYPE6,NLIST1
        DATA ITYPE1/1/,ITYPE2/2/,ITYPE3/3/,ITYPE4/4/,ITYPE5/5/,NLIST1/1/
c
c-------initialization: jbox0 is put in as the first stack member
c
        JBOX = JBOX0
        ISTACK = 1
        STACK(1,1) = 1
        STACK(2,1) = JBOX
c
        NSONS = HASKID(JBOX)
        STACK(3,1) = NSONS
c
c-------process the stack until it's empty:
c
        DO 5000 WHILE (ISTACK .GT. 0)
c
c---------check to see if this box is separated from ibox, if yes -
c           store it in list 3, continue to previous stack member;
c           if not - further checking based on if it has kid or not:
c
          CALL D3MINTERS(SIZE,CENTERS(1,IBOX),BOXES(1,IBOX),
     1      CENTERS(1,JBOX),BOXES(1,JBOX),IFINTER)
          IF (IFINTER .EQ. 0) THEN
            CALL D3MLINKSTOR(IER,ITYPE3,IBOX,JBOX,NLIST1,W,LUSED)
            IF(IER .EQ. 32) RETURN
c
            ISTACK = ISTACK-1
            STACK(3,ISTACK) = STACK(3,ISTACK)-1
            JBOX = STACK(2,ISTACK)
            GOTO 5000
          ENDIF
c
c---------this box is touching box ibox: check to see if it is childless,
c           if yes - enter it in list 1, go on to previous stack member,
c           if not - check to see if it still has unprocessed kids:
c
          IF (HASKID(JBOX) .LE. 0) THEN
c
c-----------enters jbox in the list1 of ibox; if jbox is on the finer
c             level than ibox - also enter ibox in the list 1 of jbox
c
            CALL D3MLINKSTOR(IER,ITYPE1,IBOX,JBOX,NLIST1,W,LUSED)
c
            IF (IER .EQ. 32) RETURN
            IF (BOXES(1,JBOX) .GT. BOXES(1,IBOX)) THEN
              CALL D3MLINKSTOR(IER,ITYPE1,JBOX,IBOX,NLIST1,W,LUSED)
              IF (IER .EQ. 32) RETURN
            ENDIF
c
c-----------if we have processed the whole box jbox0, get out
c             of the subroutine
c
            IF (JBOX .EQ. JBOX0) RETURN
c
            ISTACK = ISTACK-1
            STACK(3,ISTACK) = STACK(3,ISTACK)-1
            JBOX = STACK(2,ISTACK)
            GOTO 5000
          ENDIF
c
c---------this box is touching box ibox, and has kids. check to see if
c           the number of unprocessed sons of this box is zero,
c           if yes - continue to previous stack member,
c           if not - contruct the stack member for the son, and process it
c
          NSONS = STACK(3,ISTACK)
          IF (NSONS .LT. 1) THEN
c
            IF (JBOX .EQ. JBOX0) RETURN
c
c
c-----------this box is not separated from ibox, and has children,
c             store this box in list 6
c
            ITYPE6=6
c
            CALL D3MLINKSTOR(IER,ITYPE6,IBOX,JBOX,NLIST1,W,LUSED)
            CALL D3MLINKSTOR(IER,ITYPE6,JBOX,IBOX,NLIST1,W,LUSED)
c
            ISTACK = ISTACK-1
            STACK(3,ISTACK) = STACK(3,ISTACK)-1
            JBOX = STACK(2,ISTACK)
            GOTO 5000
          ENDIF
c
c---------this box is touching box ibox and has unprocessed kids,
c           construct a stack element for next unprocessed son, and
c           process that son instead:
c
          NNSONS = 0
          DO KK = 6, 13
            IF (BOXES(KK,JBOX) .GT. 0) NNSONS = NNSONS + 1
            IF (NNSONS .EQ. NSONS) GOTO 555
          ENDDO
 555      JBOX = BOXES(KK,JBOX)
c
          ISTACK = ISTACK+1
          NSONS = HASKID(JBOX)
          STACK(1,ISTACK) = ISTACK
          STACK(2,ISTACK) = JBOX
          STACK(3,ISTACK) = NSONS
c
 5000   CONTINUE
c
        RETURN
        END
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
        SUBROUTINE D3MINTERS(SIZE,C1,LEVEL1,C2,LEVEL2,IFINTER)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
c
c purpose:
c
c   this subroutine determines if two boxes intersect (or touch).
c
c on input:
c
c   size: the side length of the whole box
c   c1: the center of the first box
c   level1: level of the first box
c   c2: the center of the second box
c   level2: level of the second box
c
c on output:
c
c   ifinter: the indicator.
c     ifinter=1 means that the boxes intersect
c     ifinter=0 means that the boxes do not intersect
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
        IMPLICIT NONE
        INTEGER *4 IFINTER,LEVEL1,LEVEL2
        REAL *8 C1(3),C2(3),SIZE
c
        REAL *8 SIDE1,SIDE1H,SIDE2,SIDE2H,EPS
        REAL *8 XMIN1,YMIN1,ZMIN1,XMAX1,YMAX1,ZMAX1
        REAL *8 XMIN2,YMIN2,ZMIN2,XMAX2,YMAX2,ZMAX2
c
c-------find the side length for both boxes:
c
        SIDE1 = SIZE/2**LEVEL1
        SIDE1H = SIDE1/2
        SIDE2 = SIZE/2**LEVEL2
        SIDE2H = SIDE2/2
c
c-------find the maximum and minimum coordinates for both boxes:
c
        XMIN1=C1(1) - SIDE1H
        YMIN1=C1(2) - SIDE1H
        ZMIN1=C1(3) - SIDE1H
c
        XMAX1=C1(1) + SIDE1H
        YMAX1=C1(2) + SIDE1H
        ZMAX1=C1(3) + SIDE1H
c
        XMIN2=C2(1) - SIDE2H
        YMIN2=C2(2) - SIDE2H
        ZMIN2=C2(3) - SIDE2H
c
        XMAX2=C2(1) + SIDE2H
        YMAX2=C2(2) + SIDE2H
        ZMAX2=C2(3) + SIDE2H
c
c-------decide if the boxes intersect:
c
        EPS=XMAX1-XMIN1
        IF (EPS .GT. YMAX1-YMIN1) EPS=YMAX1-YMIN1
        IF (EPS .GT. ZMAX1-ZMIN1) EPS=ZMAX1-ZMIN1
c
        IF (EPS .GT. XMAX2-XMIN2) EPS=XMAX2-XMIN2
        IF (EPS .GT. YMAX2-YMIN2) EPS=YMAX2-YMIN2
        IF (EPS .GT. ZMAX2-ZMIN2) EPS=ZMAX2-ZMIN2
c
        EPS=EPS/10000
c
        IFINTER=1
c
c-------the code defaults the two boxes as touching, but if any of
c         the following happens, it means they are seperated:
c
        IF (XMIN1 .GT. XMAX2+EPS) IFINTER=0
        IF (XMIN2 .GT. XMAX1+EPS) IFINTER=0
c
        IF (YMIN1 .GT. YMAX2+EPS) IFINTER=0
        IF (YMIN2 .GT. YMAX1+EPS) IFINTER=0
c
        IF (ZMIN1 .GT. ZMAX2+EPS) IFINTER=0
        IF (ZMIN2 .GT. ZMAX1+EPS) IFINTER=0
c
        RETURN
        END
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
        SUBROUTINE D3MALLB(IER,Z,NPTS,NBOX,BOXES,MAXBXS,
     1    NBOXES,IZ,LADDR,NLEV,CENTER0,SIZE,IWORK,EPSCLOSE)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
c
c purpose:
c
c   this subroutine constructs a quad-tree corresponding to the
c     user-specified collection of points in space.
c     note that the interaction lists and other lists are not
c     generated in this subroutine.
c
c on input:
c   z(3,npts): the set of given points.
c   npts: the number of elements in z
c   nbox: the maximum number of points permitted in a box on
c         the finest level. in other words, a box will be further
c         subdivided if it contains more than nbox points.
c   maxbxs: the maximum total number of boxes the subroutine
c        is permitted to create.
c
c on output:
c   ier: the error return mode:
c        ier=0 means successful execution
c        ier=4 means that the subroutine attempted to create more
c              than maxbxs boxes
c        ier=16 means that the subroutine attempted to construct more
c               than 199 levels of subdivision.
c
c   boxes: main output, an integer array dimensioned (15,nboxes).
c        each 15-element column describes one box, as follows:
c
c       1 - the level of subdivision on which this box
c             was constructed;
c       2,3,4 - the coordinates of this box among  all
c             boxes on this level, for x, y, and z directions.
c       5 - the daddy of this box, identified by its address
c             in array boxes
c       6,7,8,9,10,11,12,13 - the  list of children of this box
c             (eight of them, and the child is identified by its address
c             in the array boxes)
c       14 - the location in the array iz of the particles
c             living in this box
c       15 - the number of particles living in this box
c
c   nboxes: the total number of boxes created
c
c   iz: the integer array addressing the particles in all boxes.
c       explanation: for a box ibox, the particles living in it are:
c         (z(1,j),z(2,j),z(3,j)),(z(1,j+1),z(2,j+1),z(3,j+1)),
c         (z(1,j+2),z(2,j+2),z(3,j+3)), . . .
c         (z(1,j+nj-1),z(2,j+nj-1),z(3,j+nj-1)),
c         (z(1,j+nj),z(2,j+nj),z(3,j+nj)),
c         with j=boxes(14,ibox), and nj=boxes(15,ibox)-1
c
c   laddr: an integer array dimensioned (2,numlev), containing
c          the map of array boxes, as follows:
c
c       laddr(1,i) is the starting location in array boxes of the
c         information pertaining to level=i-1
c       laddr(2,i) is the number of boxes created on the level i-1
c
c   nlev: the maximum level number on which any boxes have
c         been created. the maximim number possible is 200.
c         it is recommended that the array laddr above be
c         dimensioned at least (2,200), in case the user underestimates
c         the number of levels required.
c
c   center0: the center of the box on the level 0, containing
c         the whole simulation
c
c   size: the side of the box on the level 0
c
c working arrays:
c
c   iwork - must be at least n+2 integer*4 elements long.
c
c note:
c   1. unlike previous version, the particles are not restricted to
c      a unit box.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
c
        IMPLICIT NONE
c
        INTEGER *4 IER,NPTS,NBOX,MAXBXS,NBOXES,NLEV
        INTEGER *4 BOXES(15,1),IZ(1),LADDR(2,1),IWORK(1)
        REAL *8 EPSCLOSE,SIZE
        REAL *8 Z(3,1),CENTER0(3)
c
c-------local variables.
c
        INTEGER *4 I,LEVEL,NLEVSON,LLL
        INTEGER *4 IS(8),NS(8),IISONS(8),JJSONS(8),KKSONS(8)
        INTEGER *4 IDADSON,MAXSON,MAXLEV
        INTEGER *4 ISON,IDAD,NUMPDAD,II,JJ,KK,IIZ,NZ
        INTEGER *4 IDAD0,IDAD1
        REAL *8 CENTER(3),XMIN,XMAX,YMIN,YMAX,ZMIN,ZMAX
        REAL *8 SIZEY,SIZEZ,XLEV
c
        DATA IISONS/1,2,1,2,1,2,1,2/,JJSONS/1,1,2,2,1,1,2,2/,
     1      KKSONS/1,1,1,1,2,2,2,2/
c
c-------execution starts:
c
        IER = 0
        IDADSON = 5
c
c-------first find the smallest cubic box containing all points:
c
        XMIN = 1.0D50
        XMAX = -XMIN
        YMIN = 1.0D50
        YMAX = -YMIN
        ZMIN = 1.0D50
        ZMAX = -ZMIN
c
        DO 1020 I=1, NPTS
          IF(Z(1,I) .LT. XMIN) XMIN=Z(1,I)
          IF(Z(1,I) .GT. XMAX) XMAX=Z(1,I)
          IF(Z(2,I) .LT. YMIN) YMIN=Z(2,I)
          IF(Z(2,I) .GT. YMAX) YMAX=Z(2,I)
          IF(Z(3,I) .LT. ZMIN) ZMIN=Z(3,I)
          IF(Z(3,I) .GT. ZMAX) ZMAX=Z(3,I)
 1020   CONTINUE
c
c-------get the biggest side length:
c
        SIZE = XMAX-XMIN
        SIZEY = YMAX-YMIN
        SIZEZ = ZMAX-ZMIN
        IF(SIZEY .GT. SIZE) SIZE=SIZEY
        IF(SIZEZ .GT. SIZE) SIZE=SIZEZ
c
c-------get the center location for the whole box:
c
        CENTER0(1) = (XMIN+XMAX)/2
        CENTER0(2) = (YMIN+YMAX)/2
        CENTER0(3) = (ZMIN+ZMAX)/2
c
c-------first box is the whole box itself, initialize its 15 elements:
c       1: level; 2-4: coordinates; 5: daddy; 6-13: children; 14: location in iz.
c       15: # or particles in box.
c
        BOXES(1,1) = 0
        BOXES(2,1) = 1
        BOXES(3,1) = 1
        BOXES(4,1) = 1
        BOXES(5,1) = 0
        DO 1040 I = 6, 13
          BOXES(I,1) = 0
1040    CONTINUE
        BOXES(14,1) = 1
        BOXES(15,1) = NPTS
c
c-------for level 0, laddr(*,1) is obvious:
c         laddr(1,i) : starting location in array boxes
c         laddr(2,i) : number of boxes created on the level i-1
c
        LADDR(1,1) = 1
        LADDR(2,1) = 1
c
c-------first let the order of points be as given:
c
        DO 1060 I = 1,NPTS
          IZ(I) = I
1060    CONTINUE
c
c-------this is the main loop of the subroutine:
c         recursively (one level after another) subdivide all
c         boxes till none are left with more than nbox particles
c
        MAXSON = MAXBXS
c
        MAXLEV = 200
c
c-------now, use the size of the box and the user specification
c         for epsclose in order to find out the maximum number of
c         levels allowed here.
c
        IF(EPSCLOSE .NE. 0D0)THEN
          XLEV = DLOG(SIZE/EPSCLOSE) / DLOG(2D0)
          MAXLEV = INT(XLEV) - 1
          MAXLEV = MIN(INT(XLEV)-1,200)
        ELSE
          MAXLEV = 200
        ENDIF
c
        ISON = 1
        NLEV = 0
c
c-------now loop over all levels.
c
        DO 3000 LEVEL = 0, MAXLEV
c
c---------the starting position in the array.
c
          LADDR(1,LEVEL+2) = LADDR(1,LEVEL+1)+LADDR(2,LEVEL+1)
          NLEVSON = 0
c
c---------all previous boxes will be daddy boxes.
c
          IDAD0 = LADDR(1,LEVEL+1)
          IDAD1 = IDAD0+LADDR(2,LEVEL+1)-1
c
c---------for each box in dad's level, do subdivision:
c
          DO 2000 IDAD = IDAD0,IDAD1
c
c-----------number of points in this box:
c
            NUMPDAD = BOXES(15,IDAD)
c
c-----------if number of points in this box is less than nbox, no need to
c             further divide it:
c
            IF (NUMPDAD .LE. NBOX) GOTO 2000
c
c-----------numbering of the box idad at this level:
c
            II = BOXES(2,IDAD)
            JJ = BOXES(3,IDAD)
            KK = BOXES(4,IDAD)
c
c-----------get the center of this box:
c
            CALL D3MCENTF(CENTER0,SIZE,LEVEL,II,JJ,KK,CENTER)
c
c-----------reorder the points in this box according to its eight children:
c
            IIZ = BOXES(14,IDAD)
            NZ  = BOXES(15,IDAD)
            CALL D3MSEPA1(CENTER,Z,IZ(IIZ),NZ,IWORK,IS,NS)
c
c-----------for the eight children just created by the routine d3msepa1,
c             check to see if they were empty, if not, put it in the array
c             boxes appropriately:
c
            DO 1000 I = 1, 8
c
c-------------if this child is empty, it shall be forgotten:
c
              IF (NS(I) .EQ. 0) GOTO 1000
c
              NLEVSON = NLEVSON+1
              ISON = ISON+1
              NLEV = LEVEL+1
c
c-------------if i had created too many boxes, bomb out:
c
              IF (ISON .GT. MAXSON) THEN
                IER = 4
                RETURN
              ENDIF
c
c-------------store in array boxes all information about this son:
c
              BOXES(1,ISON) = LEVEL+1
              BOXES(2,ISON) = (II-1)*2+IISONS(I)
              BOXES(3,ISON) = (JJ-1)*2+JJSONS(I)
              BOXES(4,ISON) = (KK-1)*2+KKSONS(I)
              BOXES(5,ISON) = IDAD
c
              DO LLL = 6, 13
                BOXES(LLL,ISON) = 0
              ENDDO
c
              BOXES(14,ISON) = IS(I)+IIZ-1
              BOXES(15,ISON) = NS(I)
c
c-------------rewrite the dad-son information to dad's place:
c
              BOXES(IDADSON+I,IDAD) = ISON
c
1000        CONTINUE
c
2000      CONTINUE
          LADDR(2,LEVEL+2) = NLEVSON
          IF (NLEVSON .EQ. 0) GOTO 4000
3000    CONTINUE
c
c-------if the previous loop ends natually and coming to here, it means
c         that we have attempted to divide deep into maxlev = 200, bomb!
c
4000    CONTINUE
c
c-------total number of boxes:
c
        NBOXES = ISON
c
        RETURN
        END
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
        SUBROUTINE D3MCENTF(CENTER0,SIZE,LEVEL,I,J,K,CENTER)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
c
c purpose:
c
c   this subroutine finds the center of the box
c   number (i,j,k) on the level level. note that the
c   box on level 0 is assumed to have the center
c   center0, and the side size.
c
c on input:
c
c   center0(3): center of the whole box.
c   size: side length of the whole box.
c   level: current level.
c   i,j,k: the box numbering at current level.
c
c on output:
c
c   center(3): center of box (i,j,k) at level level.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
c
        IMPLICIT NONE
        INTEGER *4 LEVEL,I,J,K
        REAL *8 CENTER(1),CENTER0(1)
        REAL *8 SIZE
c
        INTEGER *4 LEVEL0
        REAL *8 SIDE,SIDE2,X0,Y0,Z0
        DATA LEVEL0/-1/
c
c-------execution starts:
c
        SIDE = SIZE/2.0D0**LEVEL
        SIDE2 = SIDE/2
        X0 = CENTER0(1)-SIZE/2
        Y0 = CENTER0(2)-SIZE/2
        Z0 = CENTER0(3)-SIZE/2
        CENTER(1) = X0+(I-1)*SIDE+SIDE2
        CENTER(2) = Y0+(J-1)*SIDE+SIDE2
        CENTER(3) = Z0+(K-1)*SIDE+SIDE2
c
        RETURN
        END
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
        SUBROUTINE D3MSEPA1(CENT,Z,IZ,N,IWORK,IS,NS)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
c
c purpose:
c
c   this subroutine subdivides a box into 8 child boxes and
c     reorders the particles in a box, so that each of the 8
c     children occupies a contigious chunk of array iz.
c
c   note that we are using a strange numbering convention
c     for the children:
c
c             z
c             ^
c             |         y
c             |        /
c             |      7/    8
c             |      /
c             |    5/    6
c             |----/-----------
c             |   /
c             |  /   3     4
c             | /
c             |/   1     2
c             |-------------------------> x
c
c on input:
c
c   cent: the center of the box to be subdivided
c   z: the list of all points in the box to be subdivided
c   iz: the integer array specifying the transposition already
c       applied to the points z, before the subdivision of
c       the box into children
c   n: the total number of points in array z
c
c on output:
c
c   iz: the integer array specifying the transposition already
c       applied to the points z, after the subdivision of
c       the box into children
c   is: an integer array of length 8 containing the locations
c       of the sons in array iz
c   ns: an integer array of length 8 containig the numbers of
c       elements in the sons
c
c working arrays:
c   iwork - must be n integer elements long
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
c
        IMPLICIT NONE
c
        INTEGER *4 N
        INTEGER *4 IZ(1),IWORK(1),IS(1),NS(1)
        REAL *8 CENT(1),Z(3,1)
c
c-------local variables.
c
        INTEGER *4 N1,N2,N3,N4,N5,N6,N7,N8,N12,N34,N56,N78,N1234,N5678
        INTEGER *4 IDIR
        REAL *8 THRESH
c
c-------execution starts:
c
c-------n#: the number of particles in box #.
c
        N1 = 0
        N2 = 0
        N3 = 0
        N4 = 0
        N5 = 0
        N6 = 0
        N7 = 0
        N8 = 0
        N12 = 0
        N34 = 0
        N56 = 0
        N78 = 0
        N1234 = 0
        N5678 = 0
c
c-------divide in the z-direction:
c
        IDIR=3
        THRESH=CENT(3)
        CALL D3MSEPA0(Z,IZ,N,IDIR,THRESH,IWORK,N1234)
        N5678=N-N1234
c
c-------divide in the y-direction:
c
        IDIR=2
        THRESH=CENT(2)
        IF(N1234 .NE. 0)
     1    CALL D3MSEPA0(Z,IZ,N1234,IDIR,THRESH,IWORK,N12)
        N34=N1234-N12
c
        IF(N5678 .NE. 0)
     1    CALL D3MSEPA0(Z,IZ(N1234+1),N5678,IDIR,THRESH,IWORK,N56)
        N78=N5678-N56
c
c-------divide in the x-direction:
c
        IDIR=1
        THRESH=CENT(1)
        IF(N12 .NE. 0)
     1    CALL D3MSEPA0(Z,IZ,N12,IDIR,THRESH,IWORK,N1)
        N2=N12-N1
c
        IF(N34 .NE. 0)
     1    CALL D3MSEPA0(Z,IZ(N12+1),N34,IDIR,THRESH,IWORK,N3)
        N4=N34-N3
c
        IF(N56 .NE. 0)
     1    CALL D3MSEPA0(Z,IZ(N1234+1),N56,IDIR,THRESH,IWORK,N5)
        N6=N56-N5
c
        IF(N78 .NE. 0)
     1    CALL D3MSEPA0(Z,IZ(N1234+N56+1),N78,IDIR,THRESH,IWORK,N7)
        N8=N78-N7
c
c-------write to is(8) and ns(8):
c
        IS(1)=1
        NS(1)=N1
c
        IS(2)=IS(1)+NS(1)
        NS(2)=N2
c
        IS(3)=IS(2)+NS(2)
        NS(3)=N3
c
        IS(4)=IS(3)+NS(3)
        NS(4)=N4
c
        IS(5)=IS(4)+NS(4)
        NS(5)=N5
c
        IS(6)=IS(5)+NS(5)
        NS(6)=N6
c
        IS(7)=IS(6)+NS(6)
        NS(7)=N7
c
        IS(8)=IS(7)+NS(7)
        NS(8)=N8
c
        RETURN
        END
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
        SUBROUTINE D3MSEPA0(Z,IZ,N,IDIR,THRESH,IWORK,N1)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
c
c purpose:
c
c       subdivide the points in this box according to the directional
c       parameter idir = 1 (x), 2 (y), 3 (z).
c
c input:
c
c       z(n): array of points in this box.
c       iz(*): old ordering of points.
c       idir: which direction to divide.
c       thresh: middle point of this box in this direction.
c
c output:
c
c       iz(*): rearranged ordering of points.
c       n1: number of points in the lower half of the box.
c
c working arrays:
c
c       iwork(*)
c
c note: this subroutine explains the reason for reordering the 8 children
c       boxes.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
c
        IMPLICIT NONE
c
        INTEGER *4 N,IDIR,N1
        INTEGER *4 IZ(1),IWORK(1)
        REAL *8 THRESH,Z(3,1)
c
c-------local parameters.
c
        INTEGER *4 I1,I2,I,J
c
c-------execution starts:
c
        I1 = 0
        I2 = 0
        DO 1400 I = 1,N
          J = IZ(I)
          IF (Z(IDIR,J) .LE. THRESH) THEN
            I1 = I1+1
            IZ(I1) = J
          ELSE
            I2 = I2+1
            IWORK(I2) = J
          ENDIF
1400    CONTINUE
c
c-------copy over from temperary work space the upper half.
c
        DO 1600 I = 1,I2
          IZ(I1+I) = IWORK(I)
1600    CONTINUE
c
        N1 = I1
c
        RETURN
        END
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
        SUBROUTINE D3MHASKID(BOXES,NBOXES,HASKID)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
c
c purpose:
c
c   this subroutine produces an integer array indicates
c   whether a box has kid or not.
c
c input:
c
c   boxes(15,nboxes): see d3mallb. box information.
c   nboxes: number of boxes.
c
c output:
c
c   haskid(nboxes) - an array of flag indicates if a box has a kid
c     =0 childless
c     >0 has "haskid" kids
c
c note: in fact, haskid stores the number of kids for each box
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
        IMPLICIT NONE
c
        INTEGER *4 NBOXES
        INTEGER *4 BOXES(15,1), HASKID(1)
c
        INTEGER *4 I,J
c
c-------initialize to childless:
c
        DO 1200 I = 1, NBOXES
          HASKID(I) = 0
1200    CONTINUE
c
c-------checking one by one:
c
        DO 1600 I = 1, NBOXES
          DO 1400 J = 6, 13
            IF (BOXES(J,I) .GT. 0) THEN
              HASKID(I) = HASKID(I) + 1
            ENDIF
1400      CONTINUE
1600    CONTINUE
c
        RETURN
        END
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
        SUBROUTINE D3MCENTC(CENTER0,SIZE,BOXES,NBOXES,CENTERS,SIDES)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
c
c purpose:
c
c   this subroutine produces arrays of centers for all boxes
c     in the quad-tree structure.
c
c on input:
c
c   center0: the center of the box on the level 0, containing
c            the whole simulation
c   size: the side of the box on the level 0
c   boxes: an integer array dimensioned (15,nboxes), as produced
c          by the subroutine d3mallb (see)
c   nboxes: the total number of boxes created
c
c on output:
c
c   centers: the centers of all boxes in the array boxes
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
c
        IMPLICIT NONE
c
        INTEGER BOXES(15,1),NBOXES
        REAL *8 SIZE
        REAL *8 CENTERS(3,1),SIDES(1),CENTER0(3)
c
        INTEGER *4 I,LEVEL,II,JJ,KK
        REAL *8 X00,Y00,Z00,SIDE,SIDE2
c
c-------execution starts:
c
        X00 = CENTER0(1)-SIZE/2
        Y00 = CENTER0(2)-SIZE/2
        Z00 = CENTER0(3)-SIZE/2

        DO 1200 I = 1,NBOXES
          LEVEL = BOXES(1,I)
          SIDE = SIZE/2**LEVEL
          SIDE2 = SIDE/2
c
          II = BOXES(2,I)
          JJ = BOXES(3,I)
          KK = BOXES(4,I)
c
          CENTERS(1,I) = X00+(II-1)*SIDE+SIDE2
          CENTERS(2,I) = Y00+(JJ-1)*SIDE+SIDE2
          CENTERS(3,I) = Z00+(KK-1)*SIDE+SIDE2
          SIDES(I)=SIDE
1200    CONTINUE
c
        RETURN
        END
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
        SUBROUTINE D3MLINKINIT(IER,NBOXES,NTYPES,W,LW)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
c
c purpose:
c
c   this is the initialization entry point for the link-list
c   storage-retrieval facility. it formats the array w, to
c   be used by the entries d3mlinkstor, d3mlinkretr, linkrem,
c   d3mlinkinfo below.
c
c on input:
c
c   nboxes: the number of boxes total
c   ntypes: the number of types of lists that will be stored
c           currently set to 6.
c   lw: the total amount of space in the area w to be used for
c       storage (in integer locations)
c
c on output:
c
c   ier: error return code;
c     ier=0 means successful execution
c     ier=1024 means that the amount of space in array w is grossly
c              insufficient
c   w: the formatted area for storage
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
c
        IMPLICIT NONE
c
        INTEGER *4 IER,NBOXES,NTYPES,LW
        INTEGER *4 W(1)
c
c-------local variables.
c
        INTEGER *4 IILISTAD,IILISTS,INUMELE,INBOXES,INTYPES,ILW
        INTEGER *4 ILTOT,INUMS(20)
        INTEGER *4 ILISTADD,NLISTADD,ILISTS,NUMELE,LTOT,I
c
        DATA IILISTAD/1/,IILISTS/2/,INUMELE/3/,INBOXES/4/,INTYPES/5/,
     1     ILW/6/,ILTOT/7/, INUMS/11,12,13,14,15,16,17,18,19,20,21,
     2      22,23,24,25,26,27,28,29,30/
c
c-------allocate memory for the storage facility
c
        IER = 0
c
        ILISTADD = 32
        NLISTADD = NBOXES*NTYPES+10
c
        ILISTS = ILISTADD+NLISTADD
        NUMELE = 0
        LTOT = ILISTS+NUMELE*2+10
c
c-------if storage is not enough, bomb and exit:
c
        IF (LTOT+100 .GE. LW) THEN
          IER = 1024
          RETURN
        ENDIF
c
        DO 1200 I = 1,20
          W(INUMS(I)) = 0
1200    CONTINUE
c
c-------where the information is stored.
c       iilistad:1, iilists:2, inumele:3, inboxes:4, intypes:5
c       ilw:6, iltot:7
c
        W(IILISTAD) = ILISTADD
        W(IILISTS) = ILISTS
        W(INUMELE) = NUMELE
        W(INBOXES) = NBOXES
        W(INTYPES) = NTYPES
        W(ILW) = LW
        W(ILTOT) = LTOT
c
        CALL D3MLINKINI0(W(ILISTADD),NBOXES,NTYPES)
c
        RETURN
        END
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
        SUBROUTINE D3MLINKSTOR(IER,ITYPE,IBOX,LIST,NLIST,W,LUSED)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
c
c purpose:
c
c   this entry stores dynamically a chunk of infomation for lists
c     (ibox,itype) into the linked lists in the storage array w.
c
c on input:
c   itype: the type of the elements being stored
c   ibox: the box to which these elements corresponds
c   list: the list of positive integer elements to be stored
c   nlist: the number of elements in the array list
c   w: the storage area used by this subroutine; must be first
c      formatted by the entry d3mlinkinit of this subroutine(see above)
c
c on output:
c   ier: error return code;
c     =0 means successful execution
c     =32 means that the storage area w would be exceeded
c         by this storage request
c   w: the storage area used by this subroutine
c   lused: the number of integer elements used in array
c          w after this call.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
c
        IMPLICIT NONE
c
        INTEGER *4 IER,ITYPE,IBOX,NLIST,LUSED
        INTEGER *4 W(1),LIST(1)
c
        INTEGER *4 IILISTAD,IILISTS,INUMELE,INBOXES,INTYPES,ILW,ILTOT
        INTEGER *4 INUMS(20)
c
        DATA IILISTAD/1/,IILISTS/2/,INUMELE/3/,INBOXES/4/,INTYPES/5/,
     1     ILW/6/,ILTOT/7/, INUMS/11,12,13,14,15,16,17,18,19,20,21,
     2      22,23,24,25,26,27,28,29,30/
c
c-------execution starts:
c
        IER=0
c
c-------if memory not enough, bomb:
c
        IF (W(IILISTS)+W(INUMELE)*2+NLIST*2 .GE. W(ILW) ) THEN
          IER=32
          RETURN
        ENDIF
c
c-------otherwise, store the user-specified list in array w
c
        CALL D3MLINKSTO0(ITYPE,IBOX,LIST,NLIST,W(W(IILISTAD)),
     1      W(INBOXES),W(W(IILISTS)),W(INUMELE),W(INUMS(ITYPE)) )
c
c-------update the amount of storege used
c
        LUSED=W(IILISTS)+W(INUMELE)*2+10
c
        RETURN
        END
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
        SUBROUTINE D3MLINKRETR(IER,ITYPE,IBOX,LIST,NLIST,W,LUSED)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
c
c purpose:
c
c   this entry retrieves from the storage area  w
c     a list of type (ibox,itype): meaning the interaction lists type
c     'itype' for box 'ibox'.
c
c on input:
c   itype: the type of the elements to be retrieved
c   ibox: the box to which these elements correspond
c   w: the storage area from which the information is to be
c      retrieved
c
c on output:
c   ier - error return code;
c     =0 means successful execution, nothing to report
c     =4 means that no elements are present of the type
c        itype and the box ibox
c   list: the list of positive integer elements retrieved
c   nlist: the number of elements in the array list
c   lused: the number of integer elements used in array
c          w after this call.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
        IMPLICIT NONE
c
        INTEGER *4 IER,ITYPE,IBOX,NLIST,LUSED
        INTEGER *4 W(1),LIST(1)
c
        INTEGER *4 IILISTAD,IILISTS,INUMELE,INBOXES,INTYPES,ILW,ILTOT
        INTEGER *4 INUMS(20)
        DATA IILISTAD/1/,IILISTS/2/,INUMELE/3/,INBOXES/4/,INTYPES/5/,
     1     ILW/6/,ILTOT/7/, INUMS/11,12,13,14,15,16,17,18,19,20,21,
     2      22,23,24,25,26,27,28,29,30/
c
c-------execution starts:
c
        CALL D3MLINKRET0(IER,ITYPE,IBOX,W(W(IILISTAD)),
     1      W(W(IILISTS)),LIST,W(INBOXES),NLIST)
c
        LUSED=W(IILISTS)+W(INUMELE)*2+10
c
        RETURN
        END
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        SUBROUTINE D3MLINKSTO0(ITYPE,IBOX,LIST,NLIST,LISTADDR,
     1      NBOXES,LISTS,NUMELE,NUMTYPE)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c purpose:
c
c   this entry stores dynamically a list of box numbers of list
c     type 'itype' for box 'ibox'. this list is inputed in 'list',
c     there are 'nlist' members (however, in list building code,
c     we are only adding one at a time, effectively nlist is always
c     1). this list is to be stored in global storage space 'lists',
c     which is consisted of (ntypes = 5 x nboxes) linked lists. the
c     appropriate link handle for each type and each box is stored in
c     'listaddr'.
c
c on input:
c
c   itype: the type of lists being stored
c   ibox: the box to which these elements corresponds
c   list: the list to be stored
c   nlist: the number of elements in the array list
c   nboxes: the total number of boxes indexing the elements
c           in array lists
c   listaddr: the addressing array for the main storage array
c             lists; it is assumed that it has been formatted by a call
c             to the entry d3mlinkini0 of this subroutine (see below).
c   lists: the main storage area used by this subroutine
c   numele: the number of elements of type itype for box ibox
c           already stored in array lists on entry to this subroutine
c   numtype: the number of elements of type itype for all boxes
c            stored in array lists on entry to this subroutine
c
c on output:
c
c   listaddr: the addressing array for the main storage array
c             lists; it is assumed that it has been formatted by a call
c             to the entry d3mlinkini0 of this subroutine (see below).
c   lists: the main storage area used by this subroutine
c   numele: the number of elements stored in array lists on exit
c           from this subroutine
c   numtype: the total number of elements in all lists of the
c            type itype after this call
c
c explanation:
c   lists(1,i): the pointer points to the preceding element in
c              the list with (box,type) combination as the user-supplied
c              (ibox,itype). lists(1,i) < 0 means that this is the head of
c              the link list of its type,
c   lists(2,i): the box number on the list.
c   listaddr(ibox,itype): the handle for the link list with (ibox,
c                         itype), i.e. the tail address for the link list.
c                         listaddr(ibox,itype) < 0 means that there are no
c                         elements in the list of this type for this box.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
        IMPLICIT NONE
c
        INTEGER *4 ILAST,IBOX,ITYPE,NUMELE,NUMTYPE,NLIST,NBOXES
        INTEGER *4 LISTADDR(NBOXES,1),LISTS(2,1),LIST(1)
c
        INTEGER *4 I
c
c-------execution starts:
c
c-------get the handle for the link list (ibox, itype):
c
        ILAST = LISTADDR(IBOX,ITYPE)
c
c-------add new member to the link list:
c
        DO 1200 I = 1,NLIST
          NUMELE = NUMELE+1
          NUMTYPE = NUMTYPE+1
          LISTS(1,NUMELE) = ILAST
          LISTS(2,NUMELE) = LIST(I)
          ILAST = NUMELE
1200    CONTINUE
c
c-------update the new handle:
c
        LISTADDR(IBOX,ITYPE) = ILAST
c
        RETURN
        END
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
        SUBROUTINE D3MLINKRET0(IER,ITYPE,IBOX,LISTADDR,LISTS,LIST,
     1      NBOXES,NLIST)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
c
c purpose:
c
c   this entry retrieves from the main storage array lists
c   the list of type itype for box ibox.
c
c on input:
c   itype: the type of the elements being to be retrieved
c   ibox: the box to which these element corresponds
c   listaddr: the array of handles for the array lists;
c   nboxes: the total number of boxes indexing the elements
c           in array lists
c   lists: the main storage area used by this subroutine
c
c on output:
c   ier: error return code;
c     ier=0 means successful execution, nothing to report
c     ier=4 means that no elements are present of the type
c           itype and the box ibox
c   list: the list of (ibox,itype) being retrived
c   nlist: the number of elements in the array list
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
        IMPLICIT NONE
c
        INTEGER *4 IER,ITYPE,IBOX,NBOXES,NLIST
        INTEGER *4 LISTADDR(NBOXES,1),LISTS(2,1),LIST(1)
c
        INTEGER *4 I,J,ILAST
c
c-------execution starts:
c
        IER=0
c
c-------get the handle:
c
        ILAST=LISTADDR(IBOX,ITYPE)
c
c-------if handle indicates no such lists, so be it:
c
        IF (ILAST .LE. 0) THEN
          NLIST=0
          IER=4
          RETURN
        ENDIF
c
c-------now retrieve one by one until the whole list has been retrieved:
c         (note: the if statement inside the loop is unlikely to happen,
c          i left it in there just for insurance)
c
        NLIST=0
        DO 2400 WHILE (ILAST .GT .0)
          IF (LISTS(2,ILAST) .GT. 0) THEN
            NLIST=NLIST+1
            LIST(NLIST)=LISTS(2,ILAST)
          ENDIF
          ILAST=LISTS(1,ILAST)
2400    CONTINUE
c
c-------if i did not get anything from above, quit empty handed:
c
        IF (NLIST .LE. 0) THEN
          IER=4
          RETURN
        ENDIF
c
c-------the retrieved array needs to fliped to get the old ordering:
c
        IF (NLIST .EQ. 1) RETURN
        DO 2700 I=1,NLIST/2
          J=LIST(I)
          LIST(I)=LIST(NLIST-I+1)
          LIST(NLIST-I+1)=J
2700    CONTINUE
c
        RETURN
        END
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
        SUBROUTINE D3MLINKINI0(LISTADDR,NBOXES,NTYPES)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
c
c purpose:
c   this subroutine initializes the array listaddr to be used
c     later by d3mlinksto0.
c
c on input:
c   nboxes: the number of boxes
c   ntypes: the number of types of lists that will be stored
c                  currently set to 6.
c
c on output:
c   listaddr: the array of handles for all link lists in
c          the main storage array lists; it is initialized to -1 here
c           so it does not points to anywhere.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
c
        IMPLICIT NONE
c
        INTEGER *4 NBOXES,NTYPES
        INTEGER *4 LISTADDR(NBOXES,1)
c
        INTEGER *4 K,I
c
c-------execution starts:
c
        DO 1400 K = 1,NTYPES
          DO 1200 I = 1,NBOXES
            LISTADDR(I,K) = -1
1200      CONTINUE
1400    CONTINUE
c
        RETURN
        END
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
        SUBROUTINE D3MGETB(IER,IBOX,BOX,NKIDS,CENTER,SIDE,IWORK)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
c
c purpose:
c   this entry returns to the user the characteristics of
c   user-specified box ibox.
c
c on input:
c   ibox: the box number for which the information is desired
c   w: storage area as created by the entry d3mstrcr (see above)
c
c on output:
c   ier: the error return code.
c      ier=0 means successful execution
c      ier=4 means that ibox is either greater than the number of
c            boxes in the structure or less than 1.
c   box: an integer array dimensioned box(15). its elements
c        describe the box number ibox, as explained in d3mallb.
c   nkids: number of kids box ibox has
c   center: the center of the box number ibox
c   side: the box length.
c
c note: return to the user all information about the box ibox
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
        IMPLICIT NONE
c
        INTEGER *4 IER,IBOX,NKIDS,BOX(15),IWORK(1)
        REAL *8 SIDE,CENTER(3)
c
        INTEGER *4 I,IIBOXES,IICENTERS,IICORNERS,IIWLISTS,IIKIDS,IBOX0
        INTEGER *4 IKID0
        DATA IIBOXES/1/,IICENTERS/2/,IICORNERS/3/,IIWLISTS/4/,IIKIDS/5/
        SAVE
c
        IER = 0
        IBOX0 = IWORK(IIBOXES)+(IBOX-1)*15-1
        DO 4400 I = 1,15
           BOX(I) = IWORK(IBOX0+I)
4400    CONTINUE
c
c-------return to the user the number of kids of the box ibox
c
        IKID0 = IWORK(IIKIDS)+IBOX-1
        NKIDS = IWORK(IKID0)
c
c-------return to the user the center and number of kids of the box ibox
c
        CALL D3MGETCENTR(IWORK(IWORK(IICENTERS)),IBOX,CENTER)
        CALL D3MGETSIDE(IWORK(IWORK(IICORNERS)),IBOX,SIDE)
c
        RETURN
        END
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
        SUBROUTINE D3MGETCENTR(CENTERS,IBOX,CENTER)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
c
c purpose: this subroutine returns the center of box ibox.
c
c on input:
c
c   ibox: the box number
c   centers: centers of all boxes.
c
c on output:
c
c   center: the center of box ibox.
c
c note: this subroutine can only be called by d3mgetb().
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
        IMPLICIT NONE
c
        INTEGER *4 IBOX
        REAL *8 CENTERS(3,1),CENTER(3)
c
        CENTER(1)=CENTERS(1,IBOX)
        CENTER(2)=CENTERS(2,IBOX)
        CENTER(3)=CENTERS(3,IBOX)
c
        RETURN
        END
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
        SUBROUTINE D3MGETSIDE(SIDES,IBOX,SIDE)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
c
c purpose:
c
c    this subroutine returns the side length of box ibox.
c
c on input:
c
c    ibox: the box number.
c    sides: sides of all boxes.
c
c on output:
c
c    side: side of box ibox.
c
c note: this subroutine can only be called by d3mgetb().
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
        IMPLICIT NONE
c
        INTEGER *4 IBOX
        REAL *8 SIDE,SIDES(1)
c
        SIDE=SIDES(IBOX)
c
        RETURN
        END
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
        SUBROUTINE D3MGETLIST(IER,IBOX,ITYPE,LIST,NLIST,IWORK)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccd
c
c  purpose:
c
c    return to the user the specified list of boxes for box 'ibox'.
c
c  on input:
c
c    ibox: the box number for which the information is desired
c    itype: the type of the desired list for the box ibox
c    w: storage area as created by the entry d3mstrcr (see above)
c
c  on output:
c
c    ier: the error return code.
c       ier=0 means successful execution
c       ier=4 means that the list itype for the box 'ibox  is empty
c    list: the list  itype  for the box  ibox
c    nlist: the number of elements in array  list
c
c    return to the user the list number itype for the box ibox
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
        IMPLICIT NONE
c
        INTEGER *4 IER,IBOX,ITYPE,LIST(1),NLIST,IWORK(1)
c
c-------local variables.
c
        INTEGER IIWLISTS,LUSED
        DATA IIWLISTS/4/
        SAVE
c
c-------get the box information from the correct storage location.
c
        CALL D3MLINKRETR(IER,ITYPE,IBOX,LIST,NLIST,
     $     IWORK(IWORK(IIWLISTS)),LUSED)
c
        RETURN
        END
c
