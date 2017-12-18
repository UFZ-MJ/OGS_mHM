C     ------   SUBROUTINE NFINE  ---------------------------
C
C     PURPOSE : FINE TRIANGULATION BY INSERTING NEW POINT
C
      SUBROUTINE NFINE(NODE,PX,PY,PD,NER,IFIX,NELM,MTJ,JAC,IDM,
     &               ISTACK,KTE,KTJ,KCM,KBD,
     &               NEX,IBEX,IBNO,NBK,IBREAK,NBREAK,IONL,DPP,IERR)
      IMPLICIT NONE
      INTEGER KTE,KTJ,KCM,KBD
      INTEGER NODE,IFIX(KTJ),NELM
      INTEGER MTJ(KTE,3),JAC(KTE,3),IDM(KTE),ISTACK(KTJ)
      INTEGER NEX,IBEX(KBD),IBNO(KBD,KTJ),NER(KTJ+3)
      INTEGER NBK,IBREAK(KBD),NBREAK(KBD,KTJ)
      INTEGER IONL(KTJ),IERR
      REAL*8 PX(KTJ+3),PY(KTJ+3),PD(KTJ+3),DPP
C     ==functions
      REAL*8 DLENGTH
C     ==temp
      INTEGER I,J,K
      INTEGER JADRES(KTJ+3,KBD),IAN(KTJ+3),LIST(KTJ)
      INTEGER NED(KTJ),IEDN(KTJ,KBD)
      INTEGER NP,JJ,IBK,NB,JB,IFN,NEL,N1,N2,IP1,IP2,IP3,IDIF
      INTEGER KEY,IELM,JELM,IED,JED
      REAL*8 DX,DY,DD,PD1,PD2,PD3,APD,DL1,DL2,DL3,TP,TL,DIF
      REAL*8 COM1,COM2,COM3,XP,YP


C     INITIALIZE
      DO 400 I=1,KTJ+3
        IAN(I)=0
        DO 400 J=1,KBD
          JADRES(I,J)=0
  400 CONTINUE
      NP=0
      JJ=0
      IBK=0
      DO I=1,KTJ
        NED(I)=0
        DO J=1,KBD
          IEDN(I,J)=0
        ENDDO
      ENDDO

C     MAKE LIST
      DO 410 I=1,NBK
        IBK=NP+1
        NB=I
        NP=IBREAK(NB)+IBK
        DO 420 J=IBK,NP-1
          JB=J+1-IBK
          JJ=JJ+1
          LIST(JJ)=NBREAK(NB,JB)
  420   CONTINUE
        JJ=JJ+1
C== Modified 07/3/1
        LIST(JJ)=0
C        LIST(JJ)=1
C==
  410 CONTINUE

C     MAKE JADRES
      DO 430 I=1,NP
        IF (LIST(I).GT.0) THEN
          IAN(LIST(I))=IAN(LIST(I))+1
          JADRES(LIST(I),IAN(LIST(I)))=I+1
        ENDIF
  430 CONTINUE
C==
C      DO 440 I=1,KBD
C        JADRES(1,I)=0
C  440 CONTINUE
C==

C     ADD NODE ON CLOSED BOUNDARY
      IFN=1
      DO 200 I=1,NEX
*        WRITE(*,*) '[NFINE] CLOSED BOUNDARY No = ',I
*        CALL OUTEX(NODE,NELM,MTJ,JAC,IDM,PX,PY,KTJ)
        DO 210 J=1,IBEX(I)
*          WRITE(*,*) J,'/',IBEX(I)

          N1=IBNO(I,J)
          IF(J.EQ.IBEX(I))THEN
            N2=IBNO(I,1)
          ELSE
            N2=IBNO(I,J+1)
          END IF

*          WRITE(*,*) '[NFINE] N1=',N1,' N2=',N2
*          IF (I.EQ.3) THEN
*            CALL OUTEX(NEX,NODE,NELM,MTJ,JAC,IDM,PX,PY,KTJ)
*          ENDIF

          DO K=1,NED(N1)
            IF(IEDN(N1,K).EQ.N2)THEN
              GOTO 210
            ENDIF
          ENDDO

          IF(NED(N1).EQ.KBD .OR. NED(N2).EQ.KBD)THEN
            IERR=301
            RETURN
          ENDIF
          NED(N1)=NED(N1)+1
          IEDN(N1,NED(N1))=N2
          NED(N2)=NED(N2)+1
          IEDN(N2,NED(N2))=N1

C== comment out
C          IED(N1)=IED(N1)+1
C          IED(N2)=IED(N2)+1
C          IF(N1.EQ.IBNO(I,1).AND.IED(N1).GT.1.AND.IED(N1).EQ.IED(N2))
C     &     THEN
C            GO TO 210
C          ENDIF
C          IF(N2.EQ.IBNO(I,1).AND.IED(N1).GT.2.AND.IED(N1).EQ.IED(N2))
C     &     THEN
C            GO TO 210
C          ENDIF
C          IF(IED(N1).GT.2.AND.IED(N1).EQ.IED(N2)+1)THEN
C            GO TO 210
C          ENDIF
C*          IF(N2.NE.IBNO(I,1).AND.IED(N1).NE.1.AND.IED(N2).NE.1)THEN
C*            GO TO 210
C*          ENDIF
C==
          DX=PX(N2)-PX(N1)
          DY=PY(N2)-PY(N1)
          DD=DSQRT(DX*DX+DY*DY)
          PD1=PD(N1)
          PD2=PD(N2)
          APD=(PD1+PD2)*0.5D0
          NEL=AINT((DD+0.001D0)/APD)
          IF(NEL.LT.1) GO TO 210
          IF(NEL.GT.2) GO TO 220
          NEL=2

C         ADDITION OF NEW POINT ON BOUNDARY
*          WRITE(*,*) '[NFINE] ADDONE'
          CALL ADDONE(NODE,PX,PY,NER,JADRES,DX,DY,N1,PD1,PD2,APD,
     &                PD,IAN,IFIX,NELM,MTJ,JAC,IDM,ISTACK,
     &                KTE,KTJ,KCM,KBD,NBK,IBREAK,NBREAK,NEL,IFN,IERR)
          IF (IERR.NE.0) RETURN
          KEY=1
          GO TO 210

  220     CONTINUE
*          WRITE(*,*) '[NFINE] ADDMOR',NEL
*          WRITE(*,*) '[NFINE] N1=',PX(N1),PY(N1),' N2=',PX(N2),PY(N2)
          CALL ADDMOR(NODE,PX,PY,NER,DX,DY,DD,N1,N2,PD1,PD2,APD,DPP,
     &                  PD,JADRES,IAN,IFIX,NELM,MTJ,JAC,IDM,ISTACK,
     &                  KTE,KTJ,KCM,KBD,NBK,IBREAK,NBREAK,NEL,IFN,IERR)
          IF (IERR.NE.0) RETURN
          KEY=1
  210   CONTINUE
  200 CONTINUE

C     POINT ADDITION TO BREAKLINE
*      WRITE(*,*) '[NFINE] POINT ADDITION TO OPEN BOUNDARIES'
      DO 300 I=1,NBK
*        WRITE(*,*) '[NFINE] OPEN BOUNDARY No = ',I
        DO 310 J=1,IBREAK(I)-1
*          WRITE(*,*) '[NFINE] LINE = ',NBREAK(I,J),NBREAK(I,J+1)
*          CALL OUTEX(NODE,NELM,MTJ,JAC,IDM,PX,PY,KTJ)
          N1=NBREAK(I,J)
          N2=NBREAK(I,J+1)
          PD1=PD(N1)
          PD2=PD(N2)
          DX=PX(N2)-PX(N1)
          DY=PY(N2)-PY(N1)
          DD=DSQRT(DX*DX+DY*DY)
          APD=(PD1+PD2)*0.5D0
          IF(DD.LT.APD*1.5D0) GO TO 310
          NEL=AINT(DD/APD)

*          WRITE(*,*) '[NFINE] PD1=',PD1,' PD2=',PD2
*          WRITE(*,*) '[NFINE] DD=',DD,' APD=',APD,' NEL=',NEL

          IF(NEL.LT.1) GO TO 310
          IF(NEL.GT.2) GO TO 320
          NEL=2
C         ADDITION OF NEW POINT ON BREAKLINE
*          WRITE(*,*) '[NFINE] CALL ADDONE'
          CALL ADDONE(NODE,PX,PY,NER,JADRES,DX,DY,N1,PD1,PD2,APD,
     &                  PD,IAN,IFIX,NELM,MTJ,JAC,IDM,ISTACK,
     &                  KTE,KTJ,KCM,KBD,NBK,IBREAK,NBREAK,NEL,IFN,IERR)
          IF (IERR.NE.0) RETURN
          GO TO 310

  320     CONTINUE
*          WRITE(*,*) '[NFINE] CALL ADDMOR'
          CALL ADDMOR(NODE,PX,PY,NER,DX,DY,DD,N1,N2,PD1,PD2,APD,DPP,
     &                  PD,JADRES,IAN,IFIX,NELM,MTJ,JAC,IDM,ISTACK,
     &                  KTE,KTJ,KCM,KBD,NBK,IBREAK,NBREAK,NEL,IFN,IERR)
          IF (IERR.NE.0) RETURN
  310   CONTINUE
  300 CONTINUE

*      WRITE(*,*) '[NFINE] JADRES='
*      DO I=1,NODE
*        WRITE(*,*) I,(JADRES(I,J),J=1,IAN(I))
*      ENDDO
*      PAUSE

*      WRITE(*,*) '[NFINE] ADD NODE INSIDE'
*      CALL OUTEX(NODE,NELM,MTJ,JAC,IDM,PX,PY,KTJ)

   20 CONTINUE
      KEY=0
      DO 10 I=1,NELM
*        WRITE(*,*) '[NFINE] ELEMENT No = ',I,' Node-Nelm',NODE,NELM
*        CALL OUTEX(NODE,NELM,MTJ,JAC,IDM,PX,PY,KTJ)
        N1=I
        IF(IDM(N1).EQ.0)GO TO 10
        IP1=MTJ(N1,1)
        IP2=MTJ(N1,2)
        IP3=MTJ(N1,3)
        PD1=PD(IP1)
        PD2=PD(IP2)
        PD3=PD(IP3)
        TP=PD1+PD2+PD3
        DL1=DLENGTH(IP1,IP2,PX,PY,KTJ)
        DL2=DLENGTH(IP2,IP3,PX,PY,KTJ)
        DL3=DLENGTH(IP3,IP1,PX,PY,KTJ)
        TL=DL1+DL2+DL3
        IF(TL.LT.TP*1.5D0)GO TO 10
*        WRITE(*,*) 'IP=',IP1,IP2,IP3
*        WRITE(*,*) 'PD=',PD1,PD2,PD3
*        WRITE(*,*) 'DL=',DL1,DL2,DL3
*        IF (PD1*PD2*PD3.EQ.0) THEN
*          WRITE(*,*) 'PD IS ZERO'
*          PAUSE
*        ENDIF
        DIF=0.D0
        COM1=DL1/(PD1+PD2)
        COM2=DL2/(PD2+PD3)
        COM3=DL3/(PD3+PD1)
*        WRITE(*,*) 'COM=',COM1,COM2,COM3
        IF(COM1.GT.DIF)THEN
          DIF=COM1
          IDIF=1
        END IF
        IF(COM2.GT.DIF)THEN
          DIF=COM2
          IDIF=2
        END IF
        IF(COM3.GT.DIF)THEN
          DIF=COM3
          IDIF=3
        END IF

        IF(IDIF.EQ.1)THEN
          N1=IP1
          N2=IP2
          DX=PX(N2)-PX(N1)
          DY=PY(N2)-PY(N1)
          DD=DSQRT(DX*DX+DY*DY)
          APD=(PD1+PD2)*0.5D0
        ELSE IF(IDIF.EQ.2)THEN
          N1=IP2
          N2=IP3
          DX=PX(N2)-PX(N1)
          DY=PY(N2)-PY(N1)
          DD=DSQRT(DX*DX+DY*DY)
          APD=(PD2+PD3)*0.5D0
        ELSE
          N1=IP3
          N2=IP1
          DX=PX(N2)-PX(N1)
          DY=PY(N2)-PY(N1)
          DD=DSQRT(DX*DX+DY*DY)
          APD=(PD3+PD1)*0.5D0
        END IF
        NODE=NODE+1
        IF(NODE.GT.KTJ)THEN
          IERR=302
          RETURN
        ENDIF
        PX(NODE)=PX(N1)+(DX*PD(N1))/(PD(N1)+PD(N2))
        PY(NODE)=PY(N1)+(DY*PD(N1))/(PD(N1)+PD(N2))
        PD(NODE)=APD
        XP=PX(NODE)
        YP=PY(NODE)
*        WRITE(*,*) '[NFINE] CALL WLOCATE'
        CALL WLOCATE(XP,YP,PX,PY,MTJ,JAC,NELM,KTE,KTJ,IELM,JELM,IERR)
        IF(IERR.NE.0)RETURN
        IF(IDM(IELM).NE.IDM(JELM)) IFIX(NODE)=1
*        WRITE(*,*) '[NFINE] IELM=',IELM,' JELM=',JELM
*        WRITE(*,*) '[NFINE] CALL EDGE'
        CALL EDGE(IELM,JELM,JAC,KTE,IED,IERR)
        IF(IERR.NE.0) RETURN
*        WRITE(*,*) '[NFINE] CALL EDGE'
        CALL EDGE(JELM,IELM,JAC,KTE,JED,IERR)
        IF(IERR.NE.0) RETURN
*        WRITE(*,*) '[NFINE] CALL REMESH'
        CALL REMESH(IELM,IED,JELM,JED,NODE,PX,PY,NER,JADRES,
     &                  IAN,IFIX,NELM,MTJ,JAC,IDM,ISTACK,
     &                  KTE,KTJ,KCM,KBD,NBK,IBREAK,NBREAK,0,IERR)
        IF(IERR.NE.0)RETURN
        KEY=1
   10 CONTINUE

      IF(KEY.EQ.1) GO TO 20

      RETURN
      END
C     -------   SUBROUTINE ADDONE   ------------------------
C
C     PURPOSE :
C       ï”è„Ç…ì_Ç1Ç¬í«â¡Ç∑ÇÈ
C     ARGUMENTS : 
C
      SUBROUTINE ADDONE(NODE,PX,PY,NER,JADRES,DX,DY,N1,PD1,PD2,APD,
     &                  PD,IAN,IFIX,NELM,MTJ,JAC,IDM,ISTACK,
     &                  KTE,KTJ,KCM,KBD,NBK,IBREAK,NBREAK,NEL,IFN,IERR)
      IMPLICIT NONE
      INTEGER KTE,KTJ,KCM,KBD
      INTEGER NODE,JADRES(KTJ+3,KBD),N1
      INTEGER IAN(KTJ+3),IFIX(KTJ),NELM,MTJ(KTE,3),JAC(KTE,3),IDM(KTE)
      INTEGER ISTACK(KTJ),NBK,IBREAK(KBD),NBREAK(KBD,KTJ),NEL,IFN,IERR
      INTEGER NER(KTJ+3)
      REAL*8 PX(KTJ+3),PY(KTJ+3),PD(KTJ+3),DX,DY,PD1,PD2,APD
C     ==temp
      INTEGER IELM,JELM,IED,JED
      REAL*8 XP,YP

      NODE=NODE+1
      IF(NODE.GT.KTJ)THEN
        IERR=311
        RETURN
      ENDIF
      PX(NODE)=PX(N1)+(DX*PD1)/(PD1+PD2)
      PY(NODE)=PY(N1)+(DY*PD1)/(PD1+PD2)
      PD(NODE)=APD
      XP=PX(NODE)
      YP=PY(NODE)
      CALL WLOCATE(XP,YP,PX,PY,MTJ,JAC,NELM,KTE,KTJ,IELM,JELM,IERR)
      IF(IERR.NE.0)RETURN
*      WRITE(*,*) '[ADDMOR] IELM=',IELM,' JELM=',JELM
      CALL EDGE(IELM,JELM,JAC,KTE,IED,IERR)
      IF(IERR.NE.0)RETURN
      CALL EDGE(JELM,IELM,JAC,KTE,JED,IERR)
      IF(IERR.NE.0)RETURN
*      WRITE(*,*) '[ADDMOR] IED=',IED,' JED=',JED
      CALL REMESH(IELM,IED,JELM,JED,NODE,PX,PY,NER,JADRES,
     &                  IAN,IFIX,NELM,MTJ,JAC,IDM,ISTACK,
     &                  KTE,KTJ,KCM,KBD,NBK,IBREAK,NBREAK,0,IERR)
      IF(IERR.NE.0)RETURN
      IF(IFN.EQ.1) IFIX(NODE)=1
      RETURN
      END
C     -------   SUBROUTINE ADDMOR   ------------------------
C
C     PURPOSE :
C       ï”è„Ç…ï°êîÇÃì_Çí«â¡Ç∑ÇÈ
C     ARGUMENTS : 
C
      SUBROUTINE ADDMOR(NODE,PX,PY,NER,DX,DY,DD,N1,N2,PD1,PD2,APD
     &                  ,DPP,PD,JADRES,IAN,IFIX,NELM,MTJ,JAC,IDM,ISTACK,
     &                  KTE,KTJ,KCM,KBD,NBK,IBREAK,NBREAK,NEL,IFN,IERR)
      IMPLICIT NONE
      INTEGER KTE,KTJ,KCM,KBD
      INTEGER NODE,N1,N2,JADRES(KTJ+3,KBD)
      INTEGER IAN(KTJ+3),IFIX(KTJ),NELM,MTJ(KTE,3),JAC(KTE,3),IDM(KTE)
      INTEGER ISTACK(KTJ),NBK,IBREAK(KBD),NBREAK(KBD,KTJ),NEL,IFN,IERR
      INTEGER NER(KTJ+3)
      REAL*8 PX(KTJ+3),PY(KTJ+3),PD(KTJ+3),DX,DY,DD,PD1,PD2,APD,DPP
C     ==temp
      INTEGER K
      INTEGER IELM,JELM,IED,JED
      REAL*8 XP,YP,TL(NEL),DIFF
      
      DO K=1,NEL
        TL(K)=PD1+(PD2-PD1)/(NEL-1)*(K-1)
      END DO

      DO 10 K=1,NEL-1
        IF(NODE.EQ.KTJ)THEN
          IERR=321
          RETURN
        ENDIF

        XP=PX(N1)+TL(K)*DX/DD
        YP=PY(N1)+TL(K)*DY/DD
        
        DIFF=XP*XP+PX(N2)*PX(N2)-2.D0*XP*PX(N2)
        DIFF=DIFF+YP*YP+PY(N2)*PY(N2)-2.D0*YP*PY(N2)
        DIFF=DSQRT(DIFF)
        IF(DIFF.LT.PD2*0.75D0)THEN
*          WRITE(*,*) '[ADDMOR] DIFF=',DIFF,' PD1=',PD1,' PD2=',PD2
          DIFF=PX(N2)*PX(N2)+PX(NODE)*PX(NODE)
     &        -2.D0*PX(N2)*PX(NODE)
          DIFF=DIFF+PY(N2)*PY(N2)+PY(NODE)*PY(NODE)
     &        -2.D0*PY(N2)*PY(NODE)
          DIFF=DSQRT(DIFF)
*          WRITE(*,*) '[ADDMOR] DIFF=',DIFF
          IF(DIFF.GT.PD2*1.5D0)THEN
            XP=(PX(N2)+PX(NODE))*0.5D0
            YP=(PY(N2)+PY(NODE))*0.5D0
          ELSE
            GOTO 20
          ENDIF
        ENDIF

        NODE=NODE+1
        PX(NODE)=XP
        PY(NODE)=YP
        PD(NODE)=(TL(K)+TL(K+1))*0.5D0

*        WRITE(*,*) '[ADDMOR] NODE=',NODE,' K=',K
*        WRITE(*,*) '[ADDMOR] XP=',XP,' YP=',YP
*        WRITE(*,*) '[ADDMOR] CALL WLOCATE'
        CALL WLOCATE(XP,YP,PX,PY,MTJ,JAC,NELM,KTE,KTJ,IELM,JELM,IERR)
        IF(IERR.NE.0)RETURN
*        WRITE(*,*) '[ADDMOR] IELM=',IELM,' JELM=',JELM
*        WRITE(*,*) '[ADDMOR] CALL EDGE'
        CALL EDGE(IELM,JELM,JAC,KTE,IED,IERR)
        IF(IERR.NE.0)RETURN
*        WRITE(*,*) '[ADDMOR] CALL EDGE'
        CALL EDGE(JELM,IELM,JAC,KTE,JED,IERR)
        IF(IERR.NE.0)RETURN
*        WRITE(*,*) '[ADDMOR] CALL REMESH'
        CALL REMESH(IELM,IED,JELM,JED,NODE,PX,PY,NER,JADRES,
     &                  IAN,IFIX,NELM,MTJ,JAC,IDM,ISTACK,
     &                  KTE,KTJ,KCM,KBD,NBK,IBREAK,NBREAK,0,IERR)
        IF(IERR.NE.0)RETURN
        N1=NODE
        IF(IFN.EQ.1) IFIX(NODE)=1
   10 CONTINUE

   20 CONTINUE

      RETURN
      END
