C     ------   SUBROUTINE INPUT   --------------------------
C
C     PURPOSE : INPUT OF DATA
C     LAST MODIFIED : 25 JAN 1990
C
      SUBROUTINE INPUT(NEX,NBK,IBEX,IBREAK,IBNO,NBREAK,NIDM,
     &                 NOB,NIB,PX,PY,PD,DPP,STL,KBD,KTJ)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION IBEX(KBD),IBREAK(KBD),NIDM(KBD),PX(KTJ+3),PY(KTJ+3)
      DIMENSION IBNO(KBD,KTJ),PD(KTJ+3),NBREAK(KBD,KTJ)
      CHARACTER FNAME*30
C
      DO 5 I=1,KTJ+3
        PD(I)=0.D0
    5 CONTINUE
      WRITE(*,600)
  600 FORMAT(' INPUT FILE NAME = ? ')
      READ(*,500)FNAME
  500 FORMAT(A30)
      OPEN(8,FILE=FNAME)
      READ(8,*)NEX,NBK
*      READ(8,510)NEX,NBK
      IF(NEX.NE.0)THEN
        READ(8,*)(IBEX(I),I=1,NEX)
        READ(8,*)(NIDM(I),I=1,NEX)
*        READ(8,510)(IBEX(I),I=1,NEX)
*        READ(8,510)(NIDM(I),I=1,NEX)
      END IF
      IF(NBK.NE.0)THEN
        READ(8,510)(IBREAK(I),I=1,NBK)
      ELSE
        READ(8,*)IBREAK(1)
*        READ(8,510)IBREAK(1)
      END IF
      DO 10 I=1,NEX
        READ(8,*)(IBNO(I,J),J=1,IBEX(I))
*        READ(8,510)(IBNO(I,J),J=1,IBEX(I))
   10 CONTINUE
      DO 20 I=1,NBK
        READ(8,*)(NBREAK(I,J),J=1,IBREAK(I))
*        READ(8,510)(NBREAK(I,J),J=1,IBREAK(I))
   20 CONTINUE
      READ(8,*)NOB,NIB
*      READ(8,510)NOB,NIB
      READ(8,*)(PX(I),PY(I),PD(I),I=1,NOB+NIB)
*      READ(8,520)(PX(I),PY(I),PD(I),I=1,NOB+NIB)
C     óvëfê°ñ@ÇÃéwíËÇ≥ÇÍÇƒÇ¢Ç»Ç¢êﬂì_ÇÃóvëfê°ñ@ÇÇPÇ…Ç∑ÇÈ
      DO 30 I=1,NOB+NIB
        IF(PD(I).EQ.0.D0) PD(I)=1.D0
   30 CONTINUE
      
C      DO 30 I=1,NEX
C        DO 40 J=1,IBEX(I)
C          N1=IBNO(I,J)
C          IF(PD(N1).EQ.0.D0) PD(N1)=1.D0
C   40   CONTINUE
C   30 CONTINUE
  510 FORMAT(20I5)
  520 FORMAT(3E15.7)

      CLOSE(8)

C
      WRITE(*,*) 
      WRITE(*,*) '== INPUT DATA ============================='
      WRITE(*,*) 'NUMBER OF CLOSED BOUNDARY = ',NEX
      WRITE(*,*) 'NUMBER LIST OF POINT FOR EACH CLOSED BOUNDARY : '
      WRITE(*,*) (IBEX(I),I=1,NEX)
      WRITE(*,*) 'POINT INDEX LIST FOR EACH CLOSED BOUNDARY : '
      DO I=1,NEX
        WRITE(*,*) (IBNO(I,J),J=1,IBEX(I))
      ENDDO
      WRITE(*,*) 'NUMBER OF POINT FOR EACH OPEN BOUNDARY : '
      WRITE(*,*) (NIDM(I),I=1,NEX)
      WRITE(*,*) 'NUMBER OF OPEN BOUNDARY = ', NBK
      WRITE(*,*) 'NUMBER OF POINT ON BOUNDARY = ', NOB
      WRITE(*,*) 'NUMBER OF POINT INSIDE BOUNDARY = ', NIB
      WRITE(*,*) 'POINT COORDINATE LIST:'
      WRITE(*,520)(PX(I),PY(I),PD(I),I=1,NOB+NIB)
      WRITE(*,*) '==========================================='
      WRITE(*,*) 


*      WRITE(*,610)
*  610 FORMAT(' INTERVAL BETWEEN POINT AND POINT = ? ')
*      READ(*,*)DPP
*  530 FORMAT(E15.7)
*      WRITE(*,620)
*  620 FORMAT(' STANDARD OF ELEMENT LENGTH = ? ')
*      READ(*,*)STL
      DPP=1.0D-7
      STL=1.0D0
  540 FORMAT(E15.7)

      RETURN
      END

C     ------   SUBROUTINE DATA   ---------------------------
C
C     PURPOSE : DATA TO FILE
C     LAST MODIFIED :  5 SEP 1989
C
      SUBROUTINE DATA(NEX,NODE,NELM,MTJ,JAC,IDM,PX,PY,KTJ)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION MTJ(2*KTJ+1,3),JAC(2*KTJ+1,3),IDM(2*KTJ+1)
      DIMENSION PX(KTJ+3),PY(KTJ+3)
      CHARACTER FNAME*30
C
C      WRITE(*,600)
C  600 FORMAT(/,' RESULTS TO DATA FILE ? (0-YES,1-NO) ')
C      READ(*,500)INF
C  500 FORMAT(I5)
      INF=0
      IF(INF.EQ.0)THEN
        WRITE(*,*)
        WRITE(*,*) ' OUTPUT FILE NAME = ? '
        READ(*,510)FNAME
  510   FORMAT(A30)
        OPEN(9,FILE=FNAME)
        WRITE(9,620)NODE,NELM
        WRITE(9,630)((MTJ(I,J),J=1,3),(JAC(I,J),J=1,3),
     &                    IDM(I),I=1,NELM)
        WRITE(9,640)(PX(I),PY(I),I=1,NODE)
  620   FORMAT(2I10)
  630   FORMAT(7I10)
  640   FORMAT(2E20.10)
        CLOSE(9)
      END IF
      RETURN
      END

C     ------   SUBROUTINE OUTEX   ---------------------------
C
C     PURPOSE : DATA TO FILE
C
      SUBROUTINE OUTEX(NODE,NELM,MTJ,JAC,IDM,PX,PY,KTJ)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION MTJ(2*KTJ+1,3),JAC(2*KTJ+1,3),IDM(2*KTJ+1)
      DIMENSION PX(KTJ+3),PY(KTJ+3)
      CHARACTER FNAME*30
C
      INF=0
      WRITE(0,600)
  600 FORMAT(/,' LOG TO DATA FILE ? (0-YES,1-NO) ')
      READ(*,500)INF
  500 FORMAT(I5)
      IF(INF.EQ.0)THEN
        WRITE(*,*)
        WRITE(0,*) ' OUTPUT FILE NAME = ? '
        READ(*,510)FNAME
  510   FORMAT(A30)
        OPEN(9,FILE=FNAME)
        WRITE(9,620) NODE+3,NELM
        WRITE(9,630)(I,(MTJ(I,J),J=1,3),(JAC(I,J),J=1,3),
     &                    IDM(I),I=1,NELM)
        WRITE(9,640)(I,PX(I),PY(I),I=1,NODE)
        WRITE(9,640)(I,PX(I),PY(I),I=KTJ+1,KTJ+3)
  620   FORMAT(2I10)
  630   FORMAT(8I10)
  640   FORMAT(I10,2E20.10)
        CLOSE(9)

        PAUSE
      END IF

      RETURN
      END
