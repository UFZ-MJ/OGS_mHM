C     ------   SUBROUTINE MODEL   --------------------------
C
C     PURPOSE : AUTOMATIC MESH GENERATION FOR
C               2D NON-CONVEX DOMAIN
C               BY THE MODIFIED-DELAUNAY TRIANGULATION
C     LAST MODIFIED : 22 JUN 1990
C
C     'NEX'   -  TOTAL NUMBER OF EXTERIOR BOUNDARIES
C     'NIN'   -  TOTAL NUMBER OF INTERIOR BOUNDARIES
C     'IBEX'  -  NUMBER OF NODES FOR EACH EXTERIOR BOUNDARY
C     'IBIN'  -  NUMBER OF NODES FOR EACH INTERIOR BOUNDARY
C     'IBNO'  -  RELATION BETWEEN BOUNDARY AND NODES
C     'NOB'   -  TOTAL NUMBER OF NODES ON BOUNDARY
C     'NIB'   -  TOTAL NUMBER OF NODES WITHIN BOUNDARY
C     'INDEX' -  FIRST ELEMENT NUMBER FOR EACH DOMAIN
C     'NODE'  -  TOTAL NUMBER OF NODES
C     'NELM'  -  TOTAL NUMBER OF ELEMENTS
C     'PX'    -  X-COORDINATES FOR DATA POINT
C     'PY'    -  Y-COORDINATES FOR DATA POINT
C     'MTJ'   -  RELATION BETWEEN ELEMENT AND NODES
C     'JAC'   -  RELATION BETWEEN ELEMENT AND
C                ADJACENT ELEMENTS
C     'IDM'   -  RELATION BETWEEN ELEMENT AND DOMAIN
C     KBD     -  最大境界数
C     KTJ     -  最大節点数
C     KCM     -  節点に集まる最大要素数
C     DPP     -  最小許容点間隔 etc: 1.0E-7
C     STL     -  標準要素寸法(PD=1の場合、この寸法になる)
C
      SUBROUTINE MODEL(NEX,NBK,IBEX,IBREAK,IBNO,NBREAK,NIDM,NOB,NIB,
     &                NODE,PX,PY,PD,DPP,STL,
     &                NELM,MTJ,JAC,IDM,IERR,KBD,KTJ,KCM)
      IMPLICIT NONE
C     [パラメータ]
C     JNB - 各節点の周辺要素数
C     NEI - 各節点の周辺要素番号
C      PARAMETER(KBD=15,KTJ=10000,KCM=200,KTE=2*KTJ+1)
C     ==arguments==
      INTEGER KTJ,KBD,KCM
      INTEGER NEX,NBK,IBEX(KBD),IBREAK(KBD),IBNO(KBD,KTJ)
      INTEGER NIDM(KBD),NOB,NIB,NODE,NELM
      INTEGER NBREAK(KBD,KTJ)
      INTEGER MTJ(2*KTJ+1,3),JAC(2*KTJ+1,3),IDM(2*KTJ+1),IERR
      REAL*8 PX(KTJ+3),PY(KTJ+3),PD(KTJ+3),DPP,STL
C     ==tmp variables==
      INTEGER I
      INTEGER KTE
      INTEGER NER(KTJ+3)
      INTEGER IFIX(KTJ),ISTACK(KTJ)
      INTEGER IADRES(KTJ+3)
      INTEGER IONL(KTJ),NTP,IA,IB,IC
      REAL*8 XMAX,XMIN,YMAX,YMIN,DMAX

*      WRITE(*,*) '[MODEL] START'
C     INITIALIZATION OUTPUTDATA
      IERR=0
      NODE=0
      NELM=0
      KTE=2*KTJ+1
      CALL ARRAYSET(MTJ, KTE*3, 0)
      CALL ARRAYSET(JAC, KTE*3, 0)
      CALL ARRAYSET(IDM, KTE, 0)
C     INITIALIZATION TEMP VARIABLES
      CALL ARRAYSET(NER, KTJ+3, 0)
      CALL ARRAYSET(IFIX, KTJ, 0)
      CALL ARRAYSET(ISTACK, KTJ, 0)
      CALL ARRAYSET(IADRES, KTJ+3, 0)
      CALL ARRAYSET(IONL, KTJ, 0)
      NTP=0
      IA=0
      IB=0
      IC=0
      XMAX=0.D0
      XMIN=0.D0
      YMAX=0.D0
      YMIN=0.D0
      DMAX=0.D0


C     COMPUTATION OF MAX & MIN COORDS FOR X,Y
C     NORMALIZATION OF X,Y-COORDS OF POINTS
      NTP=NOB+NIB
      CALL PANLYZ(KTJ,NTP,PX,PY,PD,XMAX,XMIN,YMAX,YMIN,DMAX)
      CALL PNRML(KTJ,NTP,PX,PY,XMIN,YMIN,DMAX)
C== modify 2007/3/18 watanabe
      DPP=DSQRT((DPP/DMAX)*(DPP/DMAX))
      DPP=DMAX1(DPP,1.0D-7)
C==
      DO 110 I=1,KTJ+3
        PD(I)=PD(I)*STL/DMAX
  110 CONTINUE

C     DEFINE VERTEX AND ADJACENCY LISTS FOR SUPERTRIANGLE
      NELM=1
      IA=KTJ+1
      IB=KTJ+2
      IC=KTJ+3
      MTJ(1,1)=IA
      MTJ(1,2)=IB
      MTJ(1,3)=IC
      JAC(1,1)=0
      JAC(1,2)=0
      JAC(1,3)=0

C     SET COORDS OF SUPERTRIANGLE
      PX(IA)=-1.23D0
      PY(IA)=-0.50D0
      PX(IB)= 2.23D0
      PY(IB)=-0.50D0
      PX(IC)= 0.50D0
      PY(IC)= 2.50D0

C     CLOSE POINT CHECK
*      WRITE(*,*) '[MODEL] POINT CHECK'
      CALL NEAR(NTP,NEX,NBK,IBEX,IBNO,IBREAK,NBREAK,PX,PY,
     &                                     DPP,KBD,KTJ)

C     ROUGH MESHES
*      WRITE(*,*) '[MODEL] ROUGH MESHES'
      CALL ROUGH(NEX,NBK,IBEX,IBREAK,IBNO,NBREAK,
     &          NTP,NIB,NODE,PX,PY,DPP,NER,NELM,MTJ,JAC,IDM,
     &          IADRES,ISTACK,KBD,KTE,KTJ,KCM,PD,IONL
     &          ,IERR)
      IF(IERR.NE.0) RETURN


*      WRITE(*,*) 'IADRES=',(IADRES(I),I=1,NODE)

*      WRITE(*,*)
*      WRITE(*,*) '== RESULT OF ROUGH ========================'
*      WRITE(*,*) 'NODE = ',NODE
*      WRITE(*,*) 'NELM = ',NELM
*      WRITE(*,*) '==========================================='
*      WRITE(*,*)

*      CALL OUTEX(NODE,NELM,MTJ,JAC,IDM,PX,PY,KTJ)

*      WRITE(*,*) '[MODEL] CHECK BOUNDARY'
      CALL CHKBND(IERR,NELM,NEX,IBEX,IBNO,NIDM,NBK,IBREAK,NBREAK,
     &  MTJ,JAC,IDM,NER,PX,PY,KBD,KTE,KTJ,KCM)
      IF (IERR.NE.0) RETURN
*      WRITE(*,*) '[MODEL] CHECK JAC'
      CALL CHECK(NELM,MTJ,JAC,KTE,IERR)
      IF (IERR.NE.0) RETURN


      CALL ARRAYSET(IFIX, NODE, 1)

C     IDM CHECK
      DO 100 I=1,NELM
        IDM(I)=0
  100 CONTINUE 
*      WRITE(*,*) '[MODEL] IDM CHECK'
      CALL DOMAIN(NELM,NEX,IBEX,IBNO,NIDM,MTJ,JAC,IDM,KBD,KTE,KTJ)

C     FINE MESHES
*      WRITE(*,*) '[MODEL] FINE MESHES'
      CALL NFINE(NODE,PX,PY,PD,NER,IFIX,NELM,MTJ,JAC,IDM,ISTACK,
     &          KTE,KTJ,KCM,KBD,
     &          NEX,IBEX,IBNO,NBK,IBREAK,NBREAK,IONL,DPP,IERR)
      IF (IERR.NE.0) RETURN

*      WRITE(*,*)
*      WRITE(*,*) '== RESULT OF FINE ========================'
*      WRITE(*,*) 'NODE = ',NODE
*      WRITE(*,*) 'NELM = ',NELM
*      WRITE(*,*) '==========================================='
*      WRITE(*,*)
*      CALL OUTEX(NODE,NELM,MTJ,JAC,IDM,PX,PY,KTJ)

C
C     CHECK OF ELEMENTS FORMED
C
      IF(NELM.NE.2*NODE+1)THEN
        IERR = 101
        RETURN
C        WRITE(*,'('' ***ERROR IN SUBROUTINE MODEL***'')')
C        WRITE(*,'('' ***INCORRECT NUMBER OF ELEMENTS***'')')
C        PAUSE
C        STOP
      END IF

C     RELOCATE NODES USING BY THE GRID SMOOTHING METHOD
*      WRITE(*,*) '[MODEL] LAPLAS'
      CALL LAPLAS2(NODE,NELM,IFIX,MTJ,JAC,PX,PY,KTE,KTJ,KCM)
*      CALL LAPLAS(NODE,NELM,IFIX,MTJ,PX,PY,JNB,NEI,KTE,KTJ,
*     &            KCM)

C     REMOVE ALL TRIANGLES OUTSIDE OF BOUNDARY
*      WRITE(*,*) '[MODEL] POST PROCESS'
      CALL REMOVE(NEX,NIDM,NELM,MTJ,JAC,IDM,KBD,KTE,IERR)
      IF(IERR.NE.0)RETURN
*      WRITE(*,*)
*      WRITE(*,*) '== RESULT OF REMOVE ========================'
*      WRITE(*,*) 'NODE = ',NODE
*      WRITE(*,*) 'NELM = ',NELM
*      WRITE(*,*) '==========================================='
*      WRITE(*,*)

C     CHECK OF RESULTS
      CALL CHECK(NELM,MTJ,JAC,KTE,IERR)
      IF (IERR.NE.0) RETURN

C     RESET X,Y-COORDS TO ORIGINAL VALUES
      CALL PDNRML(KTJ,NODE,PX,PY,XMIN,YMIN,DMAX)

      RETURN
      END
C     ------   SUBROUTINE REMOVE   -------------------------
C
C     PURPOSE : REMOVE ALL TRIANGLES OUTSIDE OF BOUNDARY
C
      SUBROUTINE REMOVE(NEX,NIDM,NELM,MTJ,JAC,IDM,KBD,KTE,IERR)
      IMPLICIT NONE
C     ==args==
      INTEGER KBD,KTE
      INTEGER NEX,NIDM(KBD),NELM,MTJ(KTE,3),JAC(KTE,3),IDM(KTE),IERR
C     ==tmp==
      INTEGER I,J,K,L
      INTEGER IZ,INELM,IELM,JELM,KELM,LELM,MKP(3),JKP(3)
      INTEGER INDEX(KBD+1),IKP,KEDG,LEDG

      IZ=0
      INELM=0
      IELM=0
      JELM=0
      KELM=0
      LELM=0
      CALL ARRAYSET(MKP, 3, 0)
      CALL ARRAYSET(JKP, 3, 0)
      CALL ARRAYSET(INDEX, KBD+1, 0)
      IKP=0
      KEDG=0
      LEDG=0

      INDEX(1)=1

      DO 10 I=1,NEX
C== modify 
        IZ=NIDM(I)
C        IZ=I
C==
        INELM=0
        DO 20 J=INDEX(I),NELM
          IF(IDM(J).EQ.IZ)THEN
            IELM=IELM+1
            JELM=J
            INELM=INELM+1
            IF(IELM.NE.JELM)THEN
              DO 30 K=1,3
                MKP(K)=MTJ(IELM,K)
                JKP(K)=JAC(IELM,K)
   30         CONTINUE
              IKP=IDM(IELM)

              DO 40 K=1,3
                KELM=JAC(JELM,K)
                IF(KELM.NE.0)THEN
                  CALL EDGE(KELM,JELM,JAC,KTE,KEDG,IERR)
                  IF (IERR.NE.0) RETURN
                  JAC(KELM,KEDG)=IELM+KTE+1
                END IF
   40         CONTINUE

              DO 50 L=1,3
                LELM=JKP(L)
                IF(LELM.NE.0)THEN
                  CALL EDGE(LELM,IELM,JAC,KTE,LEDG,IERR)
                  IF (IERR.NE.0) RETURN
                  JAC(LELM,LEDG)=JELM+KTE+1
                END IF
   50         CONTINUE

              DO 60 K=1,3
                JAC(IELM,K)=MOD(JAC(IELM,K),KTE+1)
                JAC(JELM,K)=MOD(JAC(JELM,K),KTE+1)
                KELM=JAC(IELM,K)
                LELM=JAC(JELM,K)
                DO 70 L=1,3
                  IF(KELM.NE.0) THEN
                    JAC(KELM,L)=MOD(JAC(KELM,L),KTE+1)
                  END IF
                  IF(LELM.NE.0) THEN
                    JAC(LELM,L)=MOD(JAC(LELM,L),KTE+1)
                  END IF
   70           CONTINUE
   60         CONTINUE

              DO 80 K=1,3
                JKP(K)=JAC(IELM,K)
   80         CONTINUE
              DO 90 K=1,3
                MTJ(IELM,K)=MTJ(JELM,K)
                JAC(IELM,K)=JAC(JELM,K)
                MTJ(JELM,K)=MKP(K)
                JAC(JELM,K)=JKP(K)
   90         CONTINUE
              IDM(IELM)=IDM(JELM)
              IDM(JELM)=IKP

            END IF
          END IF
   20   CONTINUE
        INDEX(I+1)=INDEX(I)+INELM
   10 CONTINUE

      DO 100 I=1,IELM
        DO 100 J=1,3
          IF(JAC(I,J).GT.IELM)THEN
            JAC(I,J)=0
          END IF
  100 CONTINUE

      DO 110 I=IELM+1,NELM
        DO 110 J=1,3
          MTJ(I,J)=0
          JAC(I,J)=0
  110 CONTINUE

      NELM=IELM

      RETURN
      END
C     ------   SUBROUTINE CHECK   --------------------------
C
C     PURPOSE : CHECK FOR ADJACENT ELEMENTS
C     LAST MODIFIED : 22 SEP 1989
C
      SUBROUTINE CHECK(NELM,MTJ,JAC,KTE,IERR)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION MTJ(KTE,3),JAC(KTE,3)
C
      DO 10 I=1,NELM
      DO 10 J=1,3
        IELM=I
        IA=MTJ(I,J)
        IB=MTJ(I,MOD(J,3)+1)
        JELM=JAC(I,J)
        IF(JELM.EQ.0)GO TO 10
        DO 20 K=1,3
          KELM=JAC(JELM,K)
          IF(KELM.EQ.IELM)THEN
            JA=MTJ(JELM,K)
            JB=MTJ(JELM,MOD(K,3)+1)
            IF((IA.EQ.JB).AND.(IB.EQ.JA))GO TO 10
            IERR=102
            RETURN
C          WRITE(*,'('' ***ERROR IN SUBROUTINE CHECK***'')')
C          WRITE(*,'('' ***JAC IS WRONG***'')')
C          PAUSE
C          STOP
          END IF
   20   CONTINUE
        IERR=103
*        WRITE(*,*)'[CHECK] ERROR: IELM=',I
        RETURN
C        WRITE(*,'(''***ERROR IN SUBROUTINE CHCK***'')')
C        WRITE(*,'(''***ADJACENT ELEMENT IS
C     &                               NOT FOUND*** '')')
C        PAUSE
C        STOP
   10 CONTINUE
      RETURN
      END
C     ------   SUBROUTINE DOMAIN   -----------------------
C
C     MODIFIED : 12 JAN 2006
C
      SUBROUTINE DOMAIN(NELM,NEX,IBEX,IBNO,NIDM,MTJ,JAC,IDM,KBD,KTE,KTJ)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION IBEX(KBD),IBNO(KBD,KTJ),NIDM(KBD)
      DIMENSION IADRES(KTJ+3),MTJ(KTE,3),JAC(KTE,3),IDM(KTE)
      DIMENSION ISTA(KTE),ISTACK(KTE),ICSTACK(KTE)

      DO 10 I=1,NEX
        NB=I
        NP=IBEX(NB)
        DO 20 J=1,KTJ+3
          IADRES(J)=0
   20   CONTINUE
        DO 30 J=1,NP
          IADRES(IBNO(NB,J))=J+1
   30   CONTINUE
        DO 40 J=1,NELM
          ISTA(J)=0
          ICSTACK(J)=0
   40   CONTINUE
        ICTOP=0
        DO 50 L=1,NELM
          MA=IADRES(MTJ(L,1))
          MB=IADRES(MTJ(L,2))
          MC=IADRES(MTJ(L,3))
*          IF(L.EQ.46)THEN
*            WRITE(*,*) 'L=',L,':',(MTJ(L,K),K=1,3)
*            WRITE(*,*) 'MABC=',MA,MB,MC
*          ENDIF
          IF((MA-MB).EQ.1 .OR. (MB-MC).EQ.1 .OR. (MC-MA).EQ.1)THEN
*            WRITE(*,*) L
            ISTA(L)=1
            IDM(L)=NB
            ICTOP=ICTOP+1
            ICSTACK(ICTOP)=L
          ELSE IF((MA.EQ.NP+1.AND.MC.EQ.2).OR.
     &      (MB.EQ.NP+1.AND.MA.EQ.2).OR.(MC.EQ.NP+1.AND.MB.EQ.2))THEN
*            WRITE(*,*) L
            ISTA(L)=1
            IDM(L)=NB
            ICTOP=ICTOP+1
            ICSTACK(ICTOP)=L
          END IF
   50   CONTINUE
C
        DO 60 J=1,ICTOP
          ITOP=0
          IS=ICSTACK(J)
*          WRITE(*,*) '[DOMAIN] IS=',IS,':',(MTJ(IS,K),K=1,3)
          MA=IADRES(MTJ(IS,1))
          MB=IADRES(MTJ(IS,2))
          MC=IADRES(MTJ(IS,3))
*          WRITE(*,*) '[DOMAIN] MABC=',MA,MB,MC
          IA=JAC(IS,1)
          MS=MA*MB
*          WRITE(*,*) '[DOMAIN] IA=',IA,':',(MTJ(IA,K),K=1,3)
          IF((MA-MB).NE.1.AND.(MA.NE.2.AND.MB.NE.NP+1).AND.
     &                                         ISTA(IA).NE.1)THEN
            ISTA(IA)=1
            ITOP=ITOP+1
            ISTACK(ITOP)=IA
          END IF
          IB=JAC(IS,2)
          MS=MB*MC
*          WRITE(*,*) '[DOMAIN] IB=',IB,':',(MTJ(IB,K),K=1,3)
          IF((MB-MC).NE.1.AND.(MB.NE.2.AND.MC.NE.NP+1).AND.
     &                                         ISTA(IB).NE.1)THEN
            ISTA(IB)=1
            ITOP=ITOP+1
            ISTACK(ITOP)=IB
          END IF
          IC=JAC(IS,3)
          MS=MC*MA
*          WRITE(*,*) '[DOMAIN] IC=',IC,':',(MTJ(IC,K),K=1,3)
          IF((MC-MA).NE.1.AND.(MC.NE.2.AND.MA.NE.NP+1).AND.
     &                                         ISTA(IC).NE.1)THEN
            ISTA(IC)=1
            ITOP=ITOP+1
            ISTACK(ITOP)=IC
          END IF

   70     IF(ITOP.GT.0)THEN
*            WRITE(*,*) '[DOMAIN] ITOP=',ITOP
            IT=ISTACK(ITOP)
            ITOP=ITOP-1
            IA=JAC(IT,1)
            IB=JAC(IT,2)
            IC=JAC(IT,3)
*            WRITE(*,*) '[DOMAIN] IT=',IT,'-',ISTA(IA),ISTA(IB),ISTA(IC)
*            WRITE(*,*) '[DOMAIN] IT=',IT,':',(MTJ(IT,K),K=1,3)
*            WRITE(*,*) '[DOMAIN] ',IA,IB,IC
            IF(ISTA(IA).NE.1)THEN
              ISTA(IA)=1
              ITOP=ITOP+1
              ISTACK(ITOP)=IA
            END IF
            IF(ISTA(IB).NE.1)THEN
              ISTA(IB)=1
              ITOP=ITOP+1
              ISTACK(ITOP)=IB
            END IF
            IF(ISTA(IC).NE.1)THEN
              ISTA(IC)=1
              ITOP=ITOP+1
              ISTACK(ITOP)=IC
            END IF
            GO TO 70
          END IF
   60   CONTINUE
        DO 80 J=1,NELM
          IF(ISTA(J).EQ.1)THEN
            IDM(J)=NIDM(NB)
          END IF
   80   CONTINUE
   10 CONTINUE
      RETURN
      END
      
