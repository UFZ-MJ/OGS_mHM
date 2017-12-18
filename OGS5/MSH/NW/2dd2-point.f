C     ------   SUBROUTINE PANLYZ  ---------------------------
C     入力点を分析する
C     LAST MODIFIED : 15 JAN 2007
C
      SUBROUTINE PANLYZ(KTJ,NTP,PX,PY,PD,XMAX,XMIN,YMAX,YMIN,DMAX)

      IMPLICIT NONE
      INTEGER KTJ,NTP
      REAL*8 PX(KTJ+3),PY(KTJ+3),PD(KTJ+3),XMAX,XMIN,YMAX,YMIN,DMAX
      INTEGER I
      REAL*8 RAX,RAY

      XMIN=PX(1)
      XMAX=XMIN
      YMIN=PY(1)
      YMAX=YMIN
      DO 50 I=2,NTP
        XMIN=DMIN1(XMIN,PX(I))
        XMAX=DMAX1(XMAX,PX(I))
        YMIN=DMIN1(YMIN,PY(I))
        YMAX=DMAX1(YMAX,PY(I))
   50 CONTINUE
      RAX=XMAX-XMIN
      RAY=YMAX-YMIN
      DMAX=DMAX1(RAX,RAY)

      RETURN
      END
C     ------   SUBROUTINE PNRML  ---------------------------
C     入力点を標準化する
C     LAST MODIFIED : 15 JAN 2007
C
      SUBROUTINE PNRML(KTJ,NTP,PX,PY,XMIN,YMIN,DMAX)

      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION PX(KTJ+3),PY(KTJ+3)

      DO 60 I=1,NTP
        PX(I)=(PX(I)-XMIN)/DMAX
        PY(I)=(PY(I)-YMIN)/DMAX
   60 CONTINUE

      RETURN
      END
C     ------   SUBROUTINE PDNRML  ---------------------------
C     入力点の標準化を解除する
C     LAST MODIFIED : 15 JAN 2007
C
      SUBROUTINE PDNRML(KTJ,NODE,PX,PY,XMIN,YMIN,DMAX)

      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION PX(KTJ+3),PY(KTJ+3)

      DO 80 I=1,NODE
        PX(I)=PX(I)*DMAX+XMIN
        PY(I)=PY(I)*DMAX+YMIN
   80 CONTINUE

      RETURN
      END
C     ------   SUBROUTINE NEAR   ---------------------------
C     近傍点を削除する
C     LAST MODIFIED : 26 JAN 2006
C
      SUBROUTINE NEAR(NTP,NEX,NBK,IBEX,IBNO,IBREAK,NBREAK,PX,PY,
     &                                           DPP,KBD,KTJ)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION IBEX(KBD),IBNO(KBD,KTJ),IBREAK(KBD),NBREAK(KBD,KTJ)
      DIMENSION PX(KTJ+3),PY(KTJ+3)

C     点同士が非常に接近していたら、同一点とみなす
      DO 10 I=1,NTP-1
        DO 20 J=I+1,NTP
          X1=PX(I)
          Y1=PY(I)
          X2=PX(J)
          Y2=PY(J)
          DIF=X1*X1+X2*X2+Y1*Y1+Y2*Y2-2*(X1*X2+Y1*Y2)
          DIF=DSQRT(DABS(DIF))
          IF(DIF.LT.DPP)THEN
*            WRITE(*,*) '[NEAR] REMOVE NEAR POINT : ', I
            CALL REWRITING(NEX,IBEX,IBNO,J,I,KBD,KTJ)
            CALL REWRITING(NBK,IBREAK,NBREAK,J,I,KBD,KTJ)
            DO 30 K=J,NTP-1
              PX(K)=PX(K+1)
              PY(K)=PY(K+1)
              CALL REWRITING(NEX,IBEX,IBNO,K+1,K,KBD,KTJ)
              CALL REWRITING(NBK,IBREAK,NBREAK,K+1,K,KBD,KTJ)
   30       CONTINUE
            NTP=NTP-1
          END IF
   20   CONTINUE
   10 CONTINUE
C     SEARCHES FOR THE DINTANCE OF THE VICINITY AND THE POINT
      DO 100 I=1,NEX
        DO 110 J=1,IBEX(I)
          IP=IBNO(I,J)
          DO 120 K=1,NEX
            IF(I.EQ.K)GO TO 120
            DO 130 L=1,IBEX(K)
              XP=PX(IP)
              YP=PY(IP)
              IA=IBNO(K,L)
              IB=IBNO(K,MOD(L,IBEX(K))+1)
              IF(IP.EQ.IA.OR.IP.EQ.IB)GO TO 130
              CALL MEA(DIS,IA,IB,XP,YP,XX,YY,PX,PY,KTJ)
              DIS=DSQRT(DABS(DIS))
              IF(DIS.LT.DPP)THEN
                PX(IP)=XX
                PY(IP)=YY
                IBEX(K)=IBEX(K)+1
                DO 140 M=IBEX(K),L+2,-1
                  IBNO(K,M)=IBNO(K,M-1)
  140           CONTINUE
                IBNO(K,L+1)=IP
              END IF
  130       CONTINUE
  120     CONTINUE
  110   CONTINUE
  100 CONTINUE
      RETURN
      END
C     ------   SUBROUTINE REWRITING   ----------------------
C     境界情報更新
C     LAST MODIFIED : 26 JAN 2006
C
      SUBROUTINE REWRITING(NUM,IBEX,IBNO,IDIS,IKEY,KBD,KTJ)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION IBEX(KBD),IBNO(KBD,KTJ)

      DO 10 I=1,NUM
      DO 10 J=1,IBEX(I)
        IF(IBNO(I,J).EQ.IDIS)THEN
          IBNO(I,J)=IKEY
        END IF
   10 CONTINUE
      RETURN
      END
