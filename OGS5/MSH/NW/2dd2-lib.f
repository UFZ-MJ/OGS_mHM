C     ------------------------------------------------------
C     PROGRAM : DELAUNAY TRIANGULATION FOR
C               ARBITRARY 2-D DOMAIN
C     LAST MODIFIED: DEC 2006 BY T.KOHARA
C     ------------------------------------------------------
      SUBROUTINE TRIANGULATE(KBD,KTJ,KCM,NEX,IBEX,IBNO,NIDM,
     &  NBK,IBREAK,NBREAK,NOB,NIB,PX,PY,PD,DPP,STL,
     &  NODE,NELM,MTJ,JAC,IDM,IERR)
      IMPLICIT REAL*8(A-H,O-Z)
C      PARAMETER(KBD=15,KTJ=10000)
      DIMENSION IBEX(KBD),IBREAK(KBD),NIDM(KBD)
      DIMENSION IBNO(KBD,KTJ),NBREAK(KBD,KTJ)
      DIMENSION PX(KTJ+3),PY(KTJ+3),PD(KTJ+3)
      DIMENSION MTJ(2*KTJ+1,3),JAC(2*KTJ+1,3),IDM(2*KTJ+1)
C
C     INITIALIZE VARIABLES
C
      IERR = 0

*      WRITE(*,*) '#INPUT DATA'
*      WRITE(*,*) 'KBD=',KBD
*      WRITE(*,*) 'KTJ=',KTJ
*      WRITE(*,*) 'KCM=',KCM
*      WRITE(*,*) 'NEX=',NEX
*      WRITE(*,*) 'IBEX=',(IBEX(I),I=1,NEX)
*      DO I=1,NEX
*        WRITE(*,*) 'IBNO:',(IBNO(I,J),J=1,IBEX(I))
*      ENDDO
*      WRITE(*,*) 'NIDM=',(NIDM(I),I=1,NEX)
*      WRITE(*,*) 'NBK=',NBK
*      WRITE(*,*) 'IBREAK=',(IBREAK(I),I=1,NBK)
*      DO I=1,NBK
*        WRITE(*,*) 'NBREAK:',(NBREAK(I,J),J=1,IBREAK(I))
*      ENDDO
*      WRITE(*,*) 'NOB=',NOB
*      WRITE(*,*) 'NIB=',NIB
*      DO I=1,NOB+NIB
*        WRITE(*,'(3e15.7)') PX(I),PY(I),PD(I)
*      ENDDO
*      WRITE(*,*) 'DPP=',DPP
*      WRITE(*,*) 'STL=',STL

C     SET DEFAULT DENSITY TO UNSET NODES
      DO 30 I=1,NEX
        DO 40 J=1,IBEX(I)
          N1=IBNO(I,J)
          IF(PD(N1).EQ.0.D0) THEN
            PD(N1)=1.D0
          ENDIF
   40   CONTINUE
   30 CONTINUE

C
C     MAIN PART OF THIS PROGRAM
C     MESH GENERATION FOR DOMAIN SURROUNDED BY BOUNDARIES
C
      CALL MODEL(NEX,NBK,IBEX,IBREAK,IBNO,NBREAK,NIDM,NOB,NIB,
     &          NODE,PX,PY,PD,DPP,STL,NELM,MTJ,JAC,IDM,IERR,KBD,KTJ,KCM) 


      RETURN
      END
