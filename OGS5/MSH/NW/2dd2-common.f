C     -------   SUBROUTINE EDGE   --------------------------
C
C     PURPOSE : FIND EDGE IN TRIANGLE L WHICH IS
C               ADJACENT TO TRIANGLE K
C
      SUBROUTINE EDGE(L,K,JAC,KTE,IEDGE,IERR)
      IMPLICIT REAL*8(A-H,O-Z) 
      DIMENSION JAC(KTE,3)
C
      DO 10 I=1,3
        IF(JAC(L,I).EQ.K)THEN
          IEDGE=I
          RETURN
        END IF
   10 CONTINUE
      IERR = 1
C      WRITE(*,'('' ***ERROR IN SUBROUTINE EDGE***'')')
C      WRITE(*,'('' ***ELEMENTS NOT ADJACENT***'')')
C      PAUSE
C      STOP
C
      RETURN
      END
C     ------   SUBROUTINE INCR   ---------------------------
C
C     PURPOSE : PLACE TRIANGLE L ON ADJACENT LIST OF POINT N
C                AND INCREMENT LIST SIZE
C     LAST MODIFIED : 22 JUN 1990
C
      SUBROUTINE INCR(N,L,NER,KTJ,KCM,IERR)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION NER(KTJ+3)

      NER(N)=L

      RETURN
      END
C     ------   SUBROUTINE DECR   ---------------------------
C
C     PURPOSE : REMOVE TRIANGLE L FROM ADJACENT LIST OF
C               POINT N
C                AND DECREMENT LIST SIZE
C     LAST MODIFIED : 22 JUN 1990
C
      SUBROUTINE DECR(N,L,MTJ,JAC,NER,KTJ,KTE,KCM,IERR)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION NER(KTJ+3),MTJ(KTE,3),JAC(KTE,3)

      IF(NER(N).EQ.L)THEN
        NER(N)=0
      ENDIF

      DO I=1,3
        JELM=JAC(L,I)
        IF(JELM.NE.0)THEN
          DO J=1,3
            IF(MTJ(JELM,J).EQ.N)THEN
              NER(N)=JELM
            ENDIF
          ENDDO
        ENDIF
      ENDDO

      RETURN
      END
CC     -------   FUNCTION THETA   ---------------------------
CC
CC     PURPOSE : COMPUTATION OF ANGLE
CC     LAST MODIFIED : 11 JUN 1990
CC
C      FUNCTION THETA(X0,Y0,X1,Y1,X2,Y2)
C      IMPLICIT REAL*8(A-H,O-Z)
C      PARAMETER(PI=3.141592653589793D0,ERROR=1.0D-15)
C
C      XA=X1-X0
C      YA=Y1-Y0
C      XB=X2-X0
C      YB=Y2-Y0
C      PRDIN=XA*XB+YA*YB
C      PRDEX=XA*YB-XB*YA
C      IF(DABS(PRDIN).LT.ERROR)THEN
C        THETA=PI*0.5D0
C      ELSE
C        THETA=DATAN(PRDEX/PRDIN)
C        IF(THETA.LT.0)THEN
C        THETA=THETA+PI
C        END IF
C        IF(THETA.GT.PI)THEN
C        THETA=THETA-PI
C        END IF
C      END IF
C      RETURN
C      END
C     -------   FUNCTION AREA   ----------------------------
C
C     PURPOSE : COMPUTATION OF AREA
C     LAST MODIFIED : 11 JUN 1990
C
      FUNCTION AREA(IELM,MTJ,PX,PY,KTJ,KTE)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION MTJ(KTE,3),PX(KTJ+3),PY(KTJ+3)
C
      J1=MTJ(IELM,1)
      J2=MTJ(IELM,2)
      J3=MTJ(IELM,3)
      AREA=0.5D0*(PX(J1)*PY(J2)+PX(J2)*PY(J3)+PX(J3)*PY(J1)
     &           -PX(J1)*PY(J3)-PX(J2)*PY(J1)-PX(J3)*PY(J2))
      RETURN
      END
C     ------   FUNCTION DLENGTH   ---------------------------
C
      FUNCTION DLENGTH(IQ,IP,PX,PY,KTJ)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION PX(KTJ+3),PY(KTJ+3)
C
      X1=PX(IQ)
      Y1=PY(IQ)
      X2=PX(IP)
      Y2=PY(IP)
      DLENGTH=DSQRT(X1*X1+X2*X2+Y1*Y1+Y2*Y2-2.D0*(X1*X2+Y1*Y2))
      RETURN
      END
C     -------   FUNCTION IPUSH   ---------------------------
C     PURPOSE : PLACE ITEM ON LIFO STACK AND
C               INCREMENT STACK SIZE
C
      FUNCTION IPUSH(ITEM,MAXSTK,ITOP,ISTACK,NNN,IERR)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION ISTACK(NNN)
C
      IF(ITOP.GT.MAXSTK)THEN
        IERR = 5
        RETURN
C        WRITE(*,'('' ***ERROR IN FUNCTION IPUSH***'')')
C        WRITE(*,'('' ***STACK OVERFLOW***'')')
C        PAUSE
C        STOP
      END IF
      IPUSH=ITEM
C
      RETURN
      END
CC     ------   SUBROUTINE SEDGE   --------------------------
CC
CC     PURPOSE : FIND TRIANGLE AND JAC WHICH HAVE POINT IP , IQ
CC
C      SUBROUTINE SEDGE(IP,IQ,KTJ,KTE,KCM,MTJ,JAC,NELM,JNB,NEI,
C     &                 IELM,JELM,IED)
C        IMPLICIT REAL*8(A-H,O-Z)
C        DIMENSION MTJ(KTE,3),JAC(KTE,3)
C        DIMENSION JNB(KTJ+3),NEI(KTJ+3,KCM)
CC
C      IELM=0
C      JELM=0
C      DO 10 I=1,JNB(IP)
C        N1=NEI(IP,I)
C        NP=0
C        NQ=0
C        DO 20 J=1,3
C          IF(MTJ(N1,J).EQ.IP)NP=J
C          IF(MTJ(N1,J).EQ.IQ)NQ=J
C   20   CONTINUE
C        IF(NP.EQ.0)GO TO 10
C        IF(NQ.EQ.0)GO TO 10
C        NP1=MOD(NQ,3)+1
C        IF(NP1.EQ.NP)THEN
C          IELM=N1
C          IED=NQ
C          IF(JELM.NE.0)RETURN
C        ELSE
C          JELM=N1
C          IF(IELM.NE.0)RETURN
C        END IF
C   10 CONTINUE
CC
C      RETURN
C        END
C     -------   FUNCTION IVERT   ---------------------------
C
C     PURPOSE : FIND VERTEX IN TRIANGLE L
C               WHICH INCLUDES POINT K
C     LAST MODIFIED :  2 JUL 1990
C
      FUNCTION IVERT(L,K,MTJ,KTE,IERR)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION MTJ(KTE,3)

      DO 10 I=1,3
        IF(MTJ(L,I).EQ.K)THEN
          IVERT=I
          RETURN
        END IF
   10 CONTINUE

      IERR=6
      RETURN
C      WRITE(*,'('' ***ERROR IN FUNCTION IVERT***'')')
C      WRITE(*,'('' ***VERTICES NOT INCLUDES***'')')
C      PAUSE
C      STOP
C
      RETURN
      END
CC     ------   SUBROUTINE LAPLAS   -------------------------
CC
CC     PURPOSE :  LAPLACIAN METHOD
CC     MODIFIED : 31 MAY 1990
CC
C      SUBROUTINE LAPLAS(NODE,NELM,IFIX,MTJ,PX,PY,JNB,NEI,
C     &                  KTE,KTJ,KCM)
C      IMPLICIT REAL*8(A-H,O-Z)
C      DIMENSION IFIX(KTJ),MTJ(KTE,3),PX(KTJ+3),PY(KTJ+3)
C      DIMENSION NEI(KTJ+3,KCM),JNB(KTJ+3)
CC
CC     ITERA : ITERATION TIMES
CC
C      ITERA=5
C      DO 10 IT=1,ITERA
C        DO 20 I=1,NODE
C        IF(IFIX(I).EQ.0)THEN
C          GX=0.D0
C          GY=0.D0
C          AR=0.D0
C          DO 30 J=1,JNB(I)
C            IELM=NEI(I,J)
C            J1=MTJ(IELM,1)
C            J2=MTJ(IELM,2)
C            J3=MTJ(IELM,3)
C            S=AREA(IELM,MTJ,PX,PY,KTJ,KTE)
C            XC=(PX(J1)+PX(J2)+PX(J3))/3.D0
C            YC=(PY(J1)+PY(J2)+PY(J3))/3.D0
C            AR=AR+S
C            GX=GX+S*XC
C            GY=GY+S*YC
C   30     CONTINUE
C          CGRAX=GX/AR
C          CGRAY=GY/AR
C          PX(I)=CGRAX
C          PY(I)=CGRAY
C        END IF
C   20   CONTINUE
C   10 CONTINUE
C      RETURN
C      END

C     -------   SUBROUTINE LAPLAS2  ------------------------
      SUBROUTINE LAPLAS2(NODE,NELM,IFIX,MTJ,JAC,PX,PY,KTE,KTJ,KCM)
      IMPLICIT NONE
      INTEGER KTE,KTJ,KCM
      INTEGER NODE,NELM,IFIX(KTJ),MTJ(KTE,3),JAC(KTE,3)
      REAL*8 PX(KTJ+3),PY(KTJ+3)
C     == functions
      REAL*8 AREA
C     == temp
      INTEGER I,J,K,L,M
      INTEGER ITRCNT,IDONE(NODE),INODE,NEIB(KCM),NBCNT,IELM,JELM
      INTEGER J1,J2,J3
      REAL*8 GX,GY,AR,S,XC,YC,CGRAX,CGRAY

      CALL ARRAYSET(NEIB,KCM,0)

      ITRCNT = 5

      DO 10 I=1,ITRCNT
        CALL ARRAYSET(IDONE,NODE,0)
        DO 20 J=1,NELM
          DO 30 K=1,3
            INODE=MTJ(J,K)
            IF(INODE.GT.KTJ) GOTO 30
            IF(IFIX(INODE).NE.0) GOTO 30
            IF(IDONE(INODE).NE.0) GOTO 30
            IDONE(INODE)=1

C           construct NEIB
            NBCNT=1
            NEIB(NBCNT)=J
            IELM=J
   40       CONTINUE
            DO 50 L=1,3
              JELM=JAC(IELM,L)
              DO M=1,NBCNT
                IF(NEIB(M).EQ.JELM) GOTO 50
              ENDDO
              DO M=1,3
                IF(MTJ(JELM,M).EQ.INODE)THEN
                  NBCNT=NBCNT+1
                  NEIB(NBCNT)=JELM
                  IELM=JELM
                  GOTO 40
                ENDIF
              ENDDO
   50       CONTINUE

C           move coordinates
            GX=0.D0
            GY=0.D0
            AR=0.D0
            DO 60 L=1,NBCNT
              IELM=NEIB(L)
              J1=MTJ(IELM,1)
              J2=MTJ(IELM,2)
              J3=MTJ(IELM,3)
              S=AREA(IELM,MTJ,PX,PY,KTJ,KTE)
              XC=(PX(J1)+PX(J2)+PX(J3))/3.D0
              YC=(PY(J1)+PY(J2)+PY(J3))/3.D0
              AR=AR+S
              GX=GX+S*XC
              GY=GY+S*YC
   60       CONTINUE
            CGRAX=GX/AR
            CGRAY=GY/AR
            PX(INODE)=CGRAX
            PY(INODE)=CGRAY

   30     CONTINUE
   20   CONTINUE
   10 CONTINUE

      RETURN
      END

C     -------   SUBROUTINE DISTOR   ------------------------
C
C     PURPOSE : COMPUTATION OF DISTORTION RATE
C     LAST MODIFIED : 11 JUN 1990
C
C      SUBROUTINE DISTOR(N,MTJ,PX,PY,DRATE,JEDG,KTE,KTJ)
C      IMPLICIT REAL*8(A-H,O-Z)
C      DIMENSION MTJ(KTE,3),PX(KTJ),PY(KTJ),DST(3)
CC
C      PI=3.1415926535
C      XA=PX(MTJ(N,1))
C      YA=PY(MTJ(N,1))
C      XB=PX(MTJ(N,2))
C      YB=PY(MTJ(N,2))
C      XC=PX(MTJ(N,3))
C      YC=PY(MTJ(N,3))
C      ANG1=THETA(XC,YC,XA,YA,XB,YB)
C      ANG2=THETA(XA,YA,XB,YB,XC,YC)
C      ANG3=PI-ANG1-ANG2
C      DST(1)=ANG1-PI/3.D0
C      DST(2)=ANG2-PI/3.D0
C      DST(3)=ANG3-PI/3.D0
C      DRATE=DABS(DST(1))+DABS(DST(2))+DABS(DST(3))
C      JEDG=1
C      IF(DST(2).GT.DST(JEDG))JEDG=2
C      IF(DST(3).GT.DST(JEDG))JEDG=3
C      RETURN
C      END
C
C     -------   SUBROUTINE SWAP   --------------------------
C
C     PURPOSE : CHECK IF POINT LIES INSIDE THE CIRCUMCIRCLE
C
      SUBROUTINE SWAP(X1,Y1,X2,Y2,X3,Y3,XP,YP,ISWAP)
      IMPLICIT REAL*8(A-H,O-Z)
C
      X13=X1-X3
      Y13=Y1-Y3
      X23=X2-X3
      Y23=Y2-Y3
      X1P=X1-XP
      Y1P=Y1-YP
      X2P=X2-XP
      Y2P=Y2-YP
      COSA=X13*X23+Y13*Y23
      COSB=X2P*X1P+Y1P*Y2P
      IF((COSA.GE.0.D0).AND.(COSB.GE.0.D0))THEN
        ISWAP=0
      ELSEIF((COSA.LT.0.D0).AND.(COSB.LT.0.D0))THEN
        ISWAP=1
      ELSE
        SINA=X13*Y23-X23*Y13
        SINB=X2P*Y1P-X1P*Y2P
        IF((SINA*COSB+SINB*COSA).LT.0.D0)THEN
          ISWAP=1
        ELSE
         ISWAP=0
        END IF
      ENDIF

      RETURN
      END
C     ------   SUBROUTINE MEA   --------------------------
C     DISを計算する.
C     MODIFIED : 23 JAN 2006
C       DIS - the distance between the line(A-B) and P
C       XX  - 傾き
C       YY  - 切片
      SUBROUTINE MEA(DIS,IA,IB,XP,YP,XX,YY,PX,PY,KTJ)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION PX(KTJ+3),PY(KTJ+3)

      NEA=0
C     点A,Bを通る直線の形状を取得する
      CALL EQUATION(IA,IB,PX,PY,A,B,IDEN,IHOR,KTJ)
C     垂直線の場合
      IF(IDEN.EQ.1)THEN
        DI=DABS(XP-A)
        IF((YP.LT.PY(IA).AND.YP.GT.PY(IB)).OR.
     &     (YP.LT.PY(IB).AND.YP.GT.PY(IA)))THEN
          DIS=DI*DI
        ELSE
          DIS=1.0D0
        END IF
        XX=A
        YY=YP
C     水平線の場合
      ELSE IF(IHOR.EQ.1)THEN
        DI=DABS(YP-B)
        IF((XP.LT.PX(IA).AND.XP.GT.PX(IB)).OR.
     &     (XP.LT.PX(IB).AND.XP.GT.PX(IA)))THEN
          DIS=DI*DI
        ELSE
          DIS=1.0D0
        END IF
        XX=XP
        YY=B
C     その他
      ELSE
        AA=-1.D0/A
        IF(XP.NE.0.0D0)THEN
          BB=YP-(AA*XP)
        ELSE
          BB=0.0D0
        END IF
        XX=(BB-B)/(A-AA)
        YY=AA*XX+BB
        IF((XX.LT.PX(IA).AND.XX.GT.PX(IB)).OR.
     &     (XX.LT.PX(IB).AND.XX.GT.PX(IA)))THEN
          DIS=XX*XX+XP*XP+YY*YY+YP*YP-2*(XX*XP+YY*YP)
        ELSE
          DIS=1.0D0
        END IF
      END IF
      RETURN
      END
C     ------   SUBROUTINE EQUATION   ---------------------
C
C     PURPOSE : 
C       指定した2点を通る直線の形状を取得する
C     ARGUMENTS : 
C       MP    - 点M
C       NP    - 点N
C       PX    - X座標データ
C       PY    - Y座標データ
C       A     - OUT：傾き
C       B     - OUT：切片
C       IDEN  - OUT：線MNが垂直か
C       IHOR  - OUT：線MNが水平か
C       KTJ   - 座標データ数
C
      SUBROUTINE EQUATION(MP,NP,PX,PY,A,B,IDEN,IHOR,KTJ)
      IMPLICIT NONE
C     ==arguments==
      INTEGER KTJ
      INTEGER MP,NP,IDEN,IHOR
      REAL*8 PX(KTJ+3),PY(KTJ+3),A,B
C     ==temp==
      REAL*8 X1,Y1,X2,Y2,DIFX

      IDEN=0
      IHOR=0
      X1=PX(MP)
      Y1=PY(MP)
      X2=PX(NP)
      Y2=PY(NP)
      DIFX=DABS(X1-X2)
      IF(DIFX.GT.1.0D-5)THEN
        A=(Y1-Y2)/(X1-X2)
        IF(DABS(A).LT.1.0D-5)THEN
          IHOR=1
        END IF
        B=Y1-A*X1
      ELSE
        IDEN=1
        A=(X1+X2)/2
        B=0.0D0
      END IF
      RETURN
      END
C     ------   SUBROUTINE WLOCATE   -------------------------
C
C     PURPOSE : LOCATE TRIANGLE WHICH ENCLOSES POINT (XP,YP)
C       点を辺上にもつ2つの要素を探す
C
      SUBROUTINE WLOCATE(XP,YP,PX,PY,MTJ,JAC,NELM,
     &                 KTE,KTJ,ITRI1,ITRI2,IERR)
      IMPLICIT NONE
      INTEGER KTE,KTJ
      INTEGER MTJ(KTE,3),JAC(KTE,3),NELM,ITRI1,ITRI2,IERR
      REAL*8 PX(KTJ+3),PY(KTJ+3),XP,YP
C     == temp
      INTEGER I,J
      INTEGER IA,IB,IC,MINTRI,KEY
      REAL*8 DIFF,DIFMIN,DMAXX,DMINX,DMAXY,DMINY,DAB,DAP,DBP

      ITRI1=0
      ITRI2=0

      CALL LOCATE(XP,YP,PX,PY,MTJ,JAC,NELM,KTE,KTJ,ITRI1,IERR)
      IF(IERR.NE.0) RETURN

      DIFMIN=1.0D10
      MINTRI=0

      DO 30 I=1,3
        IC=JAC(ITRI1,I)
        KEY=0
        DO 40 J=1,3
          IA=MTJ(IC,J)
          IB=MTJ(IC,MOD(J,3)+1)
*          WRITE(*,*) '[WLOCATE] IA2=',IA,' IB2=',IB
*          WRITE(*,*) '[WLOCATE] A=',PX(IA),PY(IA)
*          WRITE(*,*) '[WLOCATE] B=',PX(IB),PY(IB)
*          WRITE(*,*) '[WLOCATE] P=',XP,YP
C         座標値のMAX-MINで判断
          DMAXX=DMAX1(PX(IA),PX(IB))
          DMINX=DMIN1(PX(IA),PX(IB))
          DMAXY=DMAX1(PY(IA),PY(IB))
          DMINY=DMIN1(PY(IA),PY(IB))
          IF(XP.GT.DMAXX+1.0D-10 .OR. XP.LT.DMINX-1.0D-10
     &    .OR. YP.GT.DMAXY+1.0D-10 .OR. YP.LT.DMINY-1.0D-10)
     &    THEN
            GOTO 40
          ENDIF

C         三辺の内外判断
          IF((PY(IA)-YP)*(PX(IB)-XP).LE.
     &                         (PX(IA)-XP)*(PY(IB)-YP)+0.0001D0)THEN
            KEY=KEY+1
            IF(KEY.EQ.3)THEN
              ITRI2=IC
              RETURN
            END IF
          END IF

C         距離で判断
          DAB=PX(IA)*PX(IA)+PX(IB)*PX(IB)-2.D0*PX(IA)*PX(IB)
          DAB=DAB+PY(IA)*PY(IA)+PY(IB)*PY(IB)-2.D0*PY(IA)*PY(IB)
C          DAB=PX(IA)*PX(IA)-2.D0*PX(IA)*PX(IB)+PX(IB)*PX(IB)
C          DAB=DAB+PY(IA)*PY(IA)-2.D0*PY(IA)*PY(IB)+PY(IB)*PY(IB)
C          DAB=(PX(IA)-PX(IB))*(PX(IA)-PX(IB))
C          DAB=DAB+(PY(IA)-PY(IB))*(PY(IA)-PY(IB))
          DAP=PX(IA)*PX(IA)+XP*XP+PY(IA)*PY(IA)+YP*YP
          DAP=DAP-2.D0*(PX(IA)*XP+PY(IA)*YP)
C          DAP=(PX(IA)-XP)*(PX(IA)-XP)
C          DAP=DAP+(PY(IA)-YP)*(PY(IA)-YP)
          DBP=PX(IB)*PX(IB)+XP*XP+PY(IB)*PY(IB)+YP*YP
          DBP=DBP-2.D0*(PX(IB)*XP+PY(IB)*YP)
C          DBP=PX(IB)*PX(IB)-2.D0*PX(IB)*XP+XP*XP
C          DBP=DBP+PY(IB)*PY(IB)-2.D0*PY(IB)*YP+YP*YP
C          DBP=(PX(IB)-XP)*(PX(IB)-XP)
C          DBP=DBP+(PY(IB)-YP)*(PY(IB)-YP)
          DIFF=DABS(DSQRT(DAB)-(DSQRT(DAP)+DSQRT(DBP)))
*          WRITE(*,*) '[WLOCATE] DF=',DIFF
          
          IF(DIFF.LT.1.0D-10)THEN
              ITRI2=IC
              RETURN
          ENDIF

          IF(DIFF.LT.DIFMIN)THEN
            DIFMIN=DIFF
            MINTRI=IC
          ENDIF
C==
   40   CONTINUE
   30 CONTINUE

C     応急処置
      IF(MINTRI.NE.0 .AND. DIFMIN.LT.1.0D-5)THEN
        ITRI2=IC
      ENDIF

      IF(ITRI2.EQ.0) THEN
        IERR=7
        RETURN
      ENDIF

      RETURN
      END
C     ------   SUBROUTINE LOCATE   -------------------------
C
C     PURPOSE : LOCATE TRIANGLE WHICH ENCLOSES POINT (XP,YP)
C
      SUBROUTINE LOCATE(XP,YP,PX,PY,MTJ,JAC,NELM,
     &                 KTE,KTJ,ITRI,IERR)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION PX(KTJ+3),PY(KTJ+3),MTJ(KTE,3),JAC(KTE,3)

      ITRI=NELM
   10 CONTINUE
*      WRITE(*,*) '[LOCATE] ITRI=',ITRI
      DO 20 I=1,3
        IA=MTJ(ITRI,I)
        IB=MTJ(ITRI,MOD(I,3)+1)
        DIFF=(PY(IA)-YP)*(PX(IB)-XP)-(PX(IA)-XP)*(PY(IB)-YP)
*        WRITE(*,*) '[LOCATE] A=',PX(IA),PY(IA)
*        WRITE(*,*) '[LOCATE] B=',PX(IB),PY(IB)
*        WRITE(*,*) '[LOCATE] P=',XP,YP
*        WRITE(*,*) '[LOCATE] DIFF=',DIFF
        IF(DIFF.GT.1.D-10)THEN
          ITRI=JAC(ITRI,I)
          IF(ITRI.EQ.0)THEN
            IERR=1141
            RETURN
          ENDIF
          GO TO 10
        END IF
   20 CONTINUE

      RETURN
      END
C     -------   SUBROUTINE REMESH   ------------------------
C
C     PURPOSE : SUBDIVIDE GIVEN DOMAIN BY NEW DATA POINT
C               2つの三角形の対角線上に点を追加して、4つの三角形に
C               分割する
C
C     ARGUMENTS:
C       IELM    ELEMENT A ID
C       IED     ELEMENT A EDGE INDEX
C       JELM    ELEMENT B ID
C       JED     ELEMENT B EDGE INDEX
C       JADRES  
C       IAN     SIZE OF EACH JADRES ARRAY
C
      SUBROUTINE REMESH(IELM,IED,JELM,JED,NODE,PX,PY,NER,
     &                 JADRES,IAN,IFIX,NELM,MTJ,JAC,IDM,ISTACK,KTE,KTJ,
     &                 KCM,KBD,NBK,IBREAK,NBREAK,IUPNEI,IERR)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION PX(KTJ+3),PY(KTJ+3)
      DIMENSION JADRES(KTJ+3,KBD),IAN(KTJ+3),ISTACK(KTJ)
      DIMENSION JAC(KTE,3),IDM(KTE),MTJ(KTE,3),NER(KTJ+3)
      DIMENSION IBREAK(KBD),IFIX(KTJ),NBREAK(KBD,KTJ)

      ITOP=0
      MAXSTK=NODE
      IP=NODE
      XP=PX(IP)
      YP=PY(IP)

      IBJ=MOD(IED,3)+1
      JBJ=MOD(JED,3)+1
C     GET ELEMENTS AROUND A,B
      IA=JAC(IELM,IBJ)
      IB=JAC(IELM,MOD(IBJ,3)+1)
      IC=JAC(JELM,JBJ)
      ID=JAC(JELM,MOD(JBJ,3)+1)
C     GET NODE AROUND EDGE
      IV1=MTJ(IELM,IBJ)
      IV2=MTJ(IELM,MOD(IBJ,3)+1)
      IV3=MTJ(JELM,JBJ)
      IV4=MTJ(JELM,MOD(JBJ,3)+1)

*      WRITE(*,*) '[REMESH] IP=',IP,' IV1=',IV1,' IV2=',IV2,' IV3=',IV3,
*     &           ' IV4=',IV4
C     DIVIDE 2 ELEMENTS TO 4 ELEMENTS WITH IP
*      WRITE(*,*) '[REMESH] DIVIDE 2 ELEMENTS TO 4 ELEMENTS WITH IP'
      MTJ(IELM,1)=IP
      MTJ(IELM,2)=IV1
      MTJ(IELM,3)=IV2
      JAC(IELM,1)=NELM+2
      JAC(IELM,2)=IA
      JAC(IELM,3)=NELM+1

      NELM=NELM+1
      IDM(NELM)=IDM(IELM)
      MTJ(NELM,1)=IP
      MTJ(NELM,2)=IV2
      MTJ(NELM,3)=IV3
      JAC(NELM,1)=IELM
      JAC(NELM,2)=IB
      JAC(NELM,3)=JELM

      MTJ(JELM,1)=IP
      MTJ(JELM,2)=IV3
      MTJ(JELM,3)=IV4
      JAC(JELM,1)=NELM
      JAC(JELM,2)=IC
      JAC(JELM,3)=NELM+1

      NELM=NELM+1
      IDM(NELM)=IDM(JELM)
      MTJ(NELM,1)=IP
      MTJ(NELM,2)=IV4
      MTJ(NELM,3)=IV1
      JAC(NELM,1)=JELM
      JAC(NELM,2)=ID
      JAC(NELM,3)=IELM

C     RECONSTRCUT NEIBOR ELEMENTS INFO AROUND NODE
      IF(IUPNEI.NE.0)THEN
*        WRITE(*,*) '[REMESH] RECONSTRCUT NEIBOR ELEMENTS INFO'
        IF(IV1.LE.KTJ)THEN
          IF(NER(IV1).EQ.JELM)THEN
            NER(IV1)=NELM
          ENDIF
        END IF
        CALL INCR(IV2,NELM-1,NER,KTJ,KCM,IERR)
        IF (IERR.NE.0) RETURN
        IF(IV3.LE.KTJ)THEN
          IF(NER(IV3).EQ.IELM)THEN
            NER(IV3)=NELM-1
          ENDIF
        END IF
        CALL INCR(IV4,NELM,NER,KTJ,KCM,IERR)
        IF (IERR.NE.0) RETURN
        NER(IP)=IELM
      ENDIF

C     対角辺の2点が境界線上で連続しているか
      IDF=10
      DO 30 I=1,IAN(IV1)
        II=I
        DO 30 J=1,IAN(IV3)
          JJ=J
          KIDF=IABS(JADRES(IV1,I)-JADRES(IV3,J))
          IF(KIDF.LT.IDF)THEN
            IDF=KIDF
            IF(IDF.EQ.1) GO TO 50
          END IF
   30 CONTINUE

   50 CONTINUE
      MS=JADRES(IV1,1)*JADRES(IV3,1)
      IF(IDF.EQ.1 .AND. MS.NE.0) THEN
C       新規点は境界線上である
        IFIX(NODE)=1
C==
        !TODO これで正しい？
        IMAX=MAX0(JADRES(IV1,II),JADRES(IV3,JJ))
C        IF(JADRES(IV1,I).GT.JADRES(IV3,J)) THEN
C          IMAX=JADRES(IV1,I)
C        ELSE
C          IMAX=JADRES(IV3,J)
C        END IF
C==
        JADRES(NODE,1)=IMAX
        IAN(NODE)=1
C       JADRESを再構築
        IKEY=JADRES(NODE,1)
        DO 20 I=1,NODE-1
          DO 60 J=1,IAN(I)
            IF(JADRES(I,J).GE.IKEY) THEN
              JADRES(I,J)=JADRES(I,J)+1
            END IF
   60     CONTINUE
   20   CONTINUE
      END IF

C     SWAP対象を探し、スタックに蓄積する
C     ・同じ領域の三角形である
C     ・境界線を壊さない
      IF(IA.NE.0)THEN
        CALL EDGE(IA,IELM,JAC,KTE,IEDGE,IERR)
        IF (IERR.NE.0) RETURN
        IBO1=MTJ(IA,IEDGE)
        IBO2=MTJ(IA,MOD(IEDGE,3)+1)
        MS=JADRES(IBO1,1)*JADRES(IBO2,1)
        IDF=10
        DO 100 I=1,IAN(IBO1)
          DO 100 J=1,IAN(IBO2)
            KIDF=IABS(JADRES(IBO1,I)-JADRES(IBO2,J))
            IF(KIDF.LT.IDF)IDF=KIDF
  100   CONTINUE
*        WRITE(*,*) '[REMESH] IBO1=',IBO1,' IBO2=',IBO2
*        WRITE(*,*) '[REMESH] JADRES(IBO1,1)=',JADRES(IBO1,1),
*     &             ' JADRES(IBO2,1)=',JADRES(IBO2,1)
*        WRITE(*,*) '[REMESH] IDF=',IDF
        IF(IDM(IA).EQ.IDM(IELM).AND.((MS.EQ.0).OR.(IDF.NE.1)))THEN
          ITOP=ITOP+1
          ISTACK(ITOP)=IPUSH(IELM,MAXSTK,ITOP,ISTACK,KTJ,IERR)
          IF (IERR.NE.0) RETURN
        END IF
      END IF
      IF(IB.NE.0)THEN
        CALL EDGE(IB,IELM,JAC,KTE,IEDGE,IERR)
        IF (IERR.NE.0) RETURN
        JAC(IB,IEDGE)=NELM-1
        CALL EDGE(IB,NELM-1,JAC,KTE,IEDGE,IERR)
        IF (IERR.NE.0) RETURN
        IBO1=MTJ(IB,IEDGE)
        IBO2=MTJ(IB,MOD(IEDGE,3)+1)
        MS=JADRES(IBO1,1)*JADRES(IBO2,1)
        IDF=10
        DO 110 I=1,IAN(IBO1)
          DO 110 J=1,IAN(IBO2)
            KIDF=IABS(JADRES(IBO1,I)-JADRES(IBO2,J))
            IF(KIDF.LT.IDF)IDF=KIDF
  110   CONTINUE
*        WRITE(*,*) '[REMESH] IBO1=',IBO1,' IBO2=',IBO2
*        WRITE(*,*) '[REMESH] JADRES(IBO1,1)=',JADRES(IBO1,1),
*     &             ' JADRES(IBO2,1)=',JADRES(IBO2,1)
*       WRITE(*,*)'[REMESH] IDF=',IDF,'IDM1=',IDM(B),' IDM2=',IDM(NELM-1)
        IF(IDM(IB).EQ.IDM(NELM-1).AND.((MS.EQ.0).OR.(IDF.NE.1)))THEN
          ITOP=ITOP+1
          ISTACK(ITOP)=IPUSH(NELM-1,MAXSTK,ITOP,ISTACK,KTJ,IERR)
          IF (IERR.NE.0) RETURN
        END IF
      END IF
      IF(IC.NE.0)THEN
        CALL EDGE(IC,JELM,JAC,KTE,IEDGE,IERR)
        IF (IERR.NE.0) RETURN
        IBO1=MTJ(IC,IEDGE)
        IBO2=MTJ(IC,MOD(IEDGE,3)+1)
        MS=JADRES(IBO1,1)*JADRES(IBO2,1)
        IDF=10
        DO 120 I=1,IAN(IBO1)
          DO 120 J=1,IAN(IBO2)
            KIDF=IABS(JADRES(IBO1,I)-JADRES(IBO2,J))
            IF(KIDF.LT.IDF)IDF=KIDF
  120   CONTINUE
*        WRITE(*,*) '[REMESH] IBO1=',IBO1,' IBO2=',IBO2
*        WRITE(*,*) '[REMESH] JADRES(IBO1,1)=',JADRES(IBO1,1),
*     &             ' JADRES(IBO2,1)=',JADRES(IBO2,1)
*        WRITE(*,*) '[REMESH] IDF=',IDF
        IF(IDM(IC).EQ.IDM(JELM).AND.((MS.EQ.0).OR.(IDF.NE.1)))THEN
          ITOP=ITOP+1
          ISTACK(ITOP)=IPUSH(JELM,MAXSTK,ITOP,ISTACK,KTJ,IERR)
          IF (IERR.NE.0) RETURN
        END IF
      END IF
      IF(ID.NE.0)THEN
        CALL EDGE(ID,JELM,JAC,KTE,IEDGE,IERR)
        IF (IERR.NE.0) RETURN
        JAC(ID,IEDGE)=NELM
        CALL EDGE(ID,NELM,JAC,KTE,IEDGE,IERR)
        IF (IERR.NE.0) RETURN
        IBO1=MTJ(ID,IEDGE)
        IBO2=MTJ(ID,MOD(IEDGE,3)+1)
        MS=JADRES(IBO1,1)*JADRES(IBO2,1)
        IDF=10
        DO 130 I=1,IAN(IBO1)
          DO 130 J=1,IAN(IBO2)
            KIDF=IABS(JADRES(IBO1,I)-JADRES(IBO2,J))
            IF(KIDF.LT.IDF)IDF=KIDF
  130   CONTINUE
*        WRITE(*,*) '[REMESH] IBO1=',IBO1,' IBO2=',IBO2
*        WRITE(*,*) '[REMESH] JADRES(IBO1,1)=',JADRES(IBO1,1),
*     &             ' JADRES(IBO2,1)=',JADRES(IBO2,1)
*        WRITE(*,*) '[REMESH] IDF=',IDF
        IF(IDM(ID).EQ.IDM(NELM).AND.((MS.EQ.0).OR.(IDF.NE.1)))THEN
          ITOP=ITOP+1
          ISTACK(ITOP)=IPUSH(NELM,MAXSTK,ITOP,ISTACK,KTJ,IERR)
          IF (IERR.NE.0) RETURN
        END IF
      END IF

C     スタックに蓄積した要素に対して、SWAP処理をおこなう
   10 CONTINUE
      IF(ITOP.GT.0)THEN
        IL=ISTACK(ITOP)
        ITOP=ITOP-1
        IR=JAC(IL,2)

        CALL EDGE(IR,IL,JAC,KTE,IERL,IERR)
        IF (IERR.NE.0) RETURN
        IERA=MOD(IERL,3)+1
        IERB=MOD(IERA,3)+1
        IV1=MTJ(IR,IERL)
        IV2=MTJ(IR,IERA)
        IV3=MTJ(IR,IERB)
        CALL SWAP(PX(IV1),PY(IV1),PX(IV2),PY(IV2),PX(IV3),
     &            PY(IV3),XP,YP,ISWAP)
        IF(ISWAP.EQ.1)THEN
*          WRITE(*,*) '[REMESH] CALL SWAP-',IV1,IV2,IV3

          IA=JAC(IR,IERA)
          IB=JAC(IR,IERB)
          IC=JAC(IL,3)
C          WRITE(*,*) '[REMESH] IA=',IA,' IB=',IB,' IC=',IC
C          WRITE(*,*) '[REMESH] IP=',IP,' IL=',IL

          MTJ(IL,3)=IV3
          JAC(IL,2)=IA
          JAC(IL,3)=IR

          MTJ(IR,1)=IP
          MTJ(IR,2)=IV3
          MTJ(IR,3)=IV1
          JAC(IR,1)=IL
          JAC(IR,2)=IB
          JAC(IR,3)=IC

          IF(IUPNEI.NE.0)THEN
            CALL DECR(IV1,IL,MTJ,JAC,NER,KTJ,KTE,KCM,IERR)
            IF (IERR.NE.0) RETURN
            CALL DECR(IV2,IR,MTJ,JAC,NER,KTJ,KTE,KCM,IERR)
            IF (IERR.NE.0) RETURN
            CALL INCR(IP,IR,NER,KTJ,KCM,IERR)
            IF (IERR.NE.0) RETURN
            CALL INCR(IV3,IL,NER,KTJ,KCM,IERR)
            IF (IERR.NE.0) RETURN
          ENDIF

          IF(IA.NE.0)THEN
            CALL EDGE(IA,IR,JAC,KTE,IEDGE,IERR)
            IF (IERR.NE.0) RETURN
            JAC(IA,IEDGE)=IL
            CALL EDGE(IA,IL,JAC,KTE,IEDGE,IERR)
            IF (IERR.NE.0) RETURN
            IBO1=MTJ(IA,IEDGE)
            IBO2=MTJ(IA,MOD(IEDGE,3)+1)
            MS=JADRES(IBO1,1)*JADRES(IBO2,1)
            IDF=10
            DO 140 I=1,IAN(IBO1)
              DO 140 J=1,IAN(IBO2)
                KIDF=IABS(JADRES(IBO1,I)-JADRES(IBO2,J))
                IF(KIDF.LT.IDF)IDF=KIDF
  140       CONTINUE
            IF(IDM(IA).EQ.IDM(IL).AND.((MS.EQ.0).OR.(IDF.NE.1)))THEN
              ITOP=ITOP+1
              ISTACK(ITOP)=IPUSH(IL,MAXSTK,ITOP,ISTACK,KTJ,IERR)
              IF (IERR.NE.0) RETURN
            END IF
          END IF
          IF(IB.NE.0)THEN
            CALL EDGE(IB,IR,JAC,KTE,IEDGE,IERR)
            IF (IERR.NE.0) RETURN
            IBO1=MTJ(IB,IEDGE)
            IBO2=MTJ(IB,MOD(IEDGE,3)+1)
            MS=JADRES(IBO1,1)*JADRES(IBO2,1)
            IDF=10
            DO 150 I=1,IAN(IBO1)
              DO 150 J=1,IAN(IBO2)
                KIDF=IABS(JADRES(IBO1,I)-JADRES(IBO2,J))
                IF(KIDF.LT.IDF)IDF=KIDF
  150       CONTINUE
            IF(IDM(IB).EQ.IDM(IR).AND.((MS.EQ.0).OR.(IDF.NE.1)))THEN
              ITOP=ITOP+1
              ISTACK(ITOP)=IPUSH(IR,MAXSTK,ITOP,ISTACK,KTJ,IERR)
              IF (IERR.NE.0) RETURN
            END IF
          END IF
          IF(IC.NE.0)THEN
            CALL EDGE(IC,IL,JAC,KTE,IEDGE,IERR)
            IF (IERR.NE.0) RETURN
            JAC(IC,IEDGE)=IR
          END IF
        END IF
        GO TO 10
      END IF

      RETURN
      END

C     ------   SUBROUTINE ARRAYSET   ---------------------
C     配列に指定された値をセットする
C
      SUBROUTINE ARRAYSET(ARRAY, LENGTH, VALUE)

      IMPLICIT NONE
      INTEGER ARRAY(*),LENGTH,I,VALUE

      DO 10 I=1,LENGTH
        ARRAY(I)=VALUE
   10 CONTINUE

      RETURN
      END

C     ------   SUBROUTINE ARRAYSETD  ---------------------
C     配列に指定された値をセットする
C
      SUBROUTINE ARRAYSETD(ARRAY, LENGTH, VALUE)

      IMPLICIT NONE
      INTEGER LENGTH,I
      REAL*8 ARRAY(*),VALUE

      DO 10 I=1,LENGTH
        ARRAY(I)=VALUE
   10 CONTINUE

      RETURN
      END

C     ------   
C
      SUBROUTINE CHKBND(IERR,NELM,NEX,IBEX,IBNO,NIDM,NBK,IBREAK,NBREAK,
     &  MTJ,JAC,IDM,NER,PX,PY
     &  ,KBD,KTE,KTJ,KCM)

      IMPLICIT NONE
      INTEGER KBD,KTE,KTJ,KCM
      INTEGER NELM,NEX,IBEX(KBD),IBNO(KBD,KTJ),NIDM(KBD),MTJ(KTE,3)
      INTEGER JAC(KTE,3),IDM(KTE),IBREAK(KBD)
      INTEGER NBK,NBREAK(KBD,KTJ),IERR,NER(KTJ+3)
      REAL*8 PX(KTJ+3),PY(KTJ+3)
C     ==tmp
      INTEGER I,J,K,L
      INTEGER INODE,JNODE,IELM,NEI(KCM),NEICNT,CNT,IA,IB,IC
      REAL*8 XA,YA,XB,YB
      LOGICAL RET

C     閉境界のチェック
      DO I=1,NEX
*        WRITE(*,*) '[CHKBND] CLOSED BOUNDARY = ',I
        DO J=1,IBEX(I)
          INODE=IBNO(I,J)
          JNODE=IBNO(I,MOD(J,IBEX(I))+1)
          RET = .FALSE.
          CALL NEIGH(IERR,KTJ,KTE,KCM,MTJ,JAC,NER,INODE,NEI,NEICNT)
          DO K=1,NEICNT
            IELM=NEI(K)
            DO L=1,3
              IF(MTJ(IELM,L).EQ.JNODE) THEN
                RET = .TRUE.
              ENDIF
            ENDDO
          ENDDO
          IF(RET.EQV..FALSE.)THEN
*            WRITE(*,*) '[CHKBND] INODE = ',INODE,' JNODE = ',JNODE
            IERR = 7
            RETURN
          ENDIF
        ENDDO
      ENDDO

C     開境界のチェック
      DO I=1,NBK
*        WRITE(*,*) '[CHKBND] OPEN BOUNDARY = ',I
        DO J=1,IBREAK(I)-1
          INODE=NBREAK(I,J)
          JNODE=NBREAK(I,J+1)
          RET = .FALSE.
          CALL NEIGH(IERR,KTJ,KTE,KCM,MTJ,JAC,NER,INODE,NEI,NEICNT)
          DO K=1,NEICNT
            IELM=NEI(K)
            DO L=1,3
              IF(MTJ(IELM,L).EQ.JNODE) THEN
                RET = .TRUE.
              ENDIF
            ENDDO
          ENDDO
          IF(RET.EQV..FALSE.)THEN
*            WRITE(*,*) '[CHKBND] INODE = ',INODE,' JNODE = ',JNODE
            IERR = 8
            RETURN
          ENDIF
        ENDDO
      ENDDO

C     
*      WRITE(*,*) '[CHKBND] CHECK TOPOLOGY '
      CNT=0
      DO 300 I=1,NELM
        DO J=1,NELM
          CNT=0
          IF(I.NE.J)THEN
            DO K=1,3
              IF(MTJ(J,K).EQ.MTJ(I,1)) CNT=CNT+1
              IF(MTJ(J,K).EQ.MTJ(I,2)) CNT=CNT+1
            ENDDO
            IF(CNT.EQ.2) GOTO 310
          ENDIF
        ENDDO
  310   CONTINUE
        IF(CNT.NE.2 .AND. JAC(I,1).NE.0) THEN
          IERR=9
          RETURN
        ENDIF
        DO J=1,NELM
          CNT=0
          IF(I.NE.J)THEN
            DO K=1,3
              IF(MTJ(J,K).EQ.MTJ(I,1)) CNT=CNT+1
              IF(MTJ(J,K).EQ.MTJ(I,3)) CNT=CNT+1
            ENDDO
            IF(CNT.EQ.2) GOTO 320
          ENDIF
        ENDDO
  320   CONTINUE
        IF(CNT.NE.2 .AND. JAC(I,3).NE.0) THEN
          IERR=9
          RETURN
        ENDIF
        DO J=1,NELM
          CNT=0
          IF(I.NE.J)THEN
            DO K=1,3
              IF(MTJ(J,K).EQ.MTJ(I,2)) CNT=CNT+1
              IF(MTJ(J,K).EQ.MTJ(I,3)) CNT=CNT+1
            ENDDO
            IF(CNT.EQ.2) GOTO 330
          ENDIF
        ENDDO
  330   CONTINUE
        IF(CNT.NE.2 .AND. JAC(I,2).NE.0) THEN
          IERR=9
          RETURN
        ENDIF
  300 CONTINUE

C     CHECK NOT-CLOCKWISE TRIANGLE
*      WRITE(*,*) '[CHKBND] CHECK NOT-CLOCKWISE TRIANGLE '
      DO I=1,NELM
        IA=MTJ(I,1)
        IB=MTJ(I,2)
        IC=MTJ(I,3)
        XA=PX(IB)-PX(IA)
        YA=PY(IB)-PY(IA)
        XB=PX(IC)-PX(IA)
        YB=PY(IC)-PY(IA)
        IF(XA*YB-XB*YA.GT.1.0D-15)THEN
        ELSE
          IERR=10
          RETURN
        ENDIF
      ENDDO

      RETURN
      END

C     ------   
C
      SUBROUTINE NEIGH(IERR,KTJ,KTE,KCM,MTJ,JAC,NER,IP,NEI,NEICNT)
      IMPLICIT NONE
      INTEGER KTJ,KTE,KCM
      INTEGER MTJ(KTE,3),JAC(KTE,3),NER(KTJ+3),IP,NEI(KCM),NEICNT,IERR
C     == temp
      INTEGER I,J
      INTEGER IELM,JELM

      CALL ARRAYSET(NEI,KCM,0)

      IELM=NER(IP)
      IF(IELM.EQ.0)THEN
        IERR=9
        RETURN
      ENDIF
      NEICNT=1
      NEI(1)=IELM

*      WRITE(*,*) '[NEIGH] IP=',IP,' IELM=',IELM

   10 CONTINUE
        DO 30 I=1,3
          JELM=JAC(IELM,I)
*          WRITE(*,*) '[NEIGH] IELM=',IELM,' I=',I,' JELM=',JELM
          IF(JELM.EQ.0) GOTO 30
          DO J=1,NEICNT
            IF(NEI(J).EQ.JELM) GOTO 30
          ENDDO

          DO J=1,3
            IF(MTJ(JELM,J).EQ.IP)THEN
              NEICNT=NEICNT+1
              NEI(NEICNT)=JELM
              IELM=JELM
              GOTO 10
            ENDIF
          ENDDO
   30   CONTINUE

*      WRITE(*,*) '[NEIGH] NEI=',(NEI(I),I=1,NEICNT)

      RETURN
      END
C     -------   FUNCTION THETA   ---------------------------
C
C     PURPOSE : COMPUTATION OF ANGLE
C     LAST MODIFIED : 11 JUN 1990
C
      FUNCTION THETA(X0,Y0,X1,Y1,X2,Y2)
      IMPLICIT REAL*8(A-H,O-Z)

      PI=3.141592653589793D0
      ERROR=1.0D-15
      XA=X1-X0
      YA=Y1-Y0
      XB=X2-X0
      YB=Y2-Y0
      PRDIN=XA*XB+YA*YB
      DA=DSQRT(XA*XA+YA*YA)
      DB=DSQRT(XB*XB+YB*YB)

      S=PRDIN/(DA*DB)
      THETA=DACOS(S)
*      WRITE(*,*) '[THETA]',THETA

      RETURN
      END

      FUNCTION THETA2(X0,Y0,X1,Y1,X2,Y2)
      IMPLICIT REAL*8(A-H,O-Z)

      PI=3.141592653589793D0
      ERROR=1.0D-15
      XA=X1-X0
      YA=Y1-Y0
      XB=X2-X0
      YB=Y2-Y0
      PRDIN=XA*XB+YA*YB
      PRDEX=XA*YB-XB*YA
      IF(DABS(PRDIN).LT.ERROR)THEN
        THETA=PI/2.D0
      ELSE
        THETA=DATAN(PRDEX/PRDIN)
*        WRITE(*,*) '[THETA]',THETA
        IF(THETA.LT.0)THEN
          THETA=THETA+PI
        END IF
        IF(THETA.GT.PI)THEN
          THETA=THETA-PI
        END IF
        END IF
      RETURN
      END

