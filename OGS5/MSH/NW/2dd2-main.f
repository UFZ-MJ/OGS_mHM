C
C     ------   MAIN PROGRAM   ------------------------------
C
C     PURPOSE : DELAUNAY TRIANGULATION FOR
C               ARBITRARY 2-D DOMAIN
C     LAST MODIFIED : 30 JUN 1990
C
C     MAX BOUNDARIES : 10
C     MAX POINTS   : 4000
C
C     [�p�����[�^]
C     KBD - �ő勫�E��
C     KTJ - �ő�ߓ_��
C     [�ϐ�]
C     NEX - ���E��
C     NBK - �J���E��
C     IBEX - �e���E�̐ߓ_��
C     IBREAK - �e�J���E�̐ߓ_��
C     IBNO - �e���E���\������ߓ_�ԍ����X�g
C     NBREAK - �e�J���E���\������ߓ_�ԍ����X�g
C     NIDM - �e���E�̗̈�ԍ�
C     NOB - ���E��̑��ߓ_��
C     NIB - �̈�����̑��ߓ_��
C     PX,PY - �ߓ_�̍��W�l
C     PD - �ߓ_�̓_���x
C     MTJ - �v�f�̐ߓ_���
C     JAC - �v�f�̗אڏ��
C     IDM - �v�f�̗̈�ԍ�
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER(KBD=15,KTJ=4500000,KCM=200)
      DIMENSION IBEX(KBD),IBREAK(KBD),NIDM(KBD)
      DIMENSION IBNO(KBD,KTJ),NBREAK(KBD,KTJ)
      DIMENSION PX(KTJ+3),PY(KTJ+3),PD(KTJ+3)
      DIMENSION MTJ(2*KTJ+1,3),JAC(2*KTJ+1,3),IDM(2*KTJ+1)

*      COMMON /CMN1/IBEX,IBREAK,MIDM,IBNO,NBREAK,MTJ,JAC,IDM
*      COMMON /CMN2/PX,PY,PD

C
C     INPUT OF DATA
C
      CALL INPUT(NEX,NBK,IBEX,IBREAK,IBNO,NBREAK,NIDM,
     &           NOB,NIB,PX,PY,PD,DPP,STL,KBD,KTJ)
C
C     MAIN PART OF THIS PROGRAM
C     MESH GENERATION FOR DOMAIN SURROUNDED BY BOUNDARIES
C
C      CALL CPU_TIME(TIME1)

      WRITE(*,*)
      WRITE(*,*) '# MAIN PROGRAM START'
      WRITE(*,*)
      CALL CPU_TIME(TIME1)
      CALL MODEL(NEX,NBK,IBEX,IBREAK,IBNO,NBREAK,NIDM,NOB,NIB,
     &  NODE,PX,PY,PD,DPP,STL,NELM,MTJ,JAC,IDM,IERR,KBD,KTJ,KCM) 
      CALL CPU_TIME(TIME2)

      IF(IERR.NE.0) THEN
        WRITE(*,*) '## ERROR OCCURED!!! -> ERROR_CODE=',IERR
        STOP
      ENDIF
*      WRITE(*,*)
*      WRITE(*,*) '# MAIN PROGRAM FINISH'
*      WRITE(*,*)
      WRITE(*,*) '# TIME = ',TIME2-TIME1

C     CALL CPU_TIME(TIME2)
C     WRITE(*,*)'TIME1=',TIME1
C     WRITE(*,*)'TIME2=',TIME2
C     WRITE(*,*)'TIME=',TIME2-TIME1
      WRITE(*,*)' TOTAL NUMBER OF NODES = ',NODE
      WRITE(*,*)' TOTAL NUMBER OF ELEMENTS = ',NELM
C
C     PRINT RESULTS TO DATA FILE
C
      CALL DATA(NEX,NODE,NELM,MTJ,JAC,IDM,PX,PY,KTJ)
C
      PAUSE
      STOP
      END
