C     ------   SUBROUTINE ROUGH   --------------------------
C
C     PURPOSE : ROUGH TRIANGULATION BY THE MODIFIED-DELAUNAY
C               METHOD
C     [MEMO]
C     IZ=0�̂Ƃ��A�O���̈�BIZ<>0�̂Ƃ��A�����̈�B
C     IADRES�́A�E���BMTJ�͍����B
C
      SUBROUTINE ROUGH(NEX,NBK,IBEX,IBREAK,IBNO,NBREAK,NTP,NIB,
     &                NODE,PX,PY,DPP,NER,NELM,MTJ,JAC,IDM,
     &                IADRES,ISTACK,KBD,KTE,KTJ,KCM,PD,IONL
     &                ,IERR)
      IMPLICIT NONE
C     ==arguments==
      INTEGER KTJ,KTE,KCM,KBD
      INTEGER NEX,NBK,IBEX(KBD),IBREAK(KBD),IBNO(KBD,KTJ)
      INTEGER NBREAK(KBD,KTJ),NTP,NIB,NODE
      INTEGER NELM,MTJ(KTE,3),JAC(KTE,3),IDM(KTE),IADRES(KTJ+3)
      INTEGER ISTACK(KTJ),IONL(KTJ),IERR,NER(KTJ+3)
      REAL*8 PX(KTJ+3),PY(KTJ+3),DPP,PD(KTJ+3)
C     ==temp==
      INTEGER I,J,K,L
      INTEGER LIST(KTJ),ICHECK(KTE),MAP(KTE),IBR(KTE),KV(KTE)
      INTEGER JADRES(KTJ+3,KBD),IAN(KTJ+3),IZ,NB,NP,MA,MB,MC,N1,IS,IELM
      INTEGER MS,MP,IT,IV1
      INTEGER NEI(KCM),NEICNT
      REAL*8 XA,YA,SUP

      CALL ARRAYSET(LIST,KTJ,0)
      CALL ARRAYSET(ICHECK,KTE,0)
      CALL ARRAYSET(MAP,KTE,0)
      CALL ARRAYSET(IBR,KTE,0)
      CALL ARRAYSET(KV,KTE,0)
      CALL ARRAYSET(JADRES,(KTJ+3)*KBD,0)
      CALL ARRAYSET(IAN,KTJ+3,0)
      IZ=0
      NB=0
      NP=0  
      MA=0
      MB=0
      MC=0
      MS=0
      MP=0
      IT=0
      XA=0.D0
      YA=0.D0

C
C     CLOSED BOUNDARIES
C
      DO 10 I=1,NEX
        IZ=0
        NB=I
*        WRITE(*,*)'[ROUGH] CLOSED BOUNDARY No = ',I
        NP=IBEX(NB)
C       �z��LIST�ɂh�Ԗڂ̋��E�_���X�g���Z�b�g����
        DO 20 J=1,NP
          LIST(J)=IBNO(NB,J)
   20   CONTINUE
C       ���E���쐬���Ȃ���A�O�p�`�������s��
*        WRITE(*,*) '[ROUGH] CALL BOUGEN'
        CALL BOUGEN(IZ,NP,LIST,NTP,PX,PY,NER,NELM,ICHECK,MTJ,
     &             JAC,IDM,IADRES,ISTACK,KV,IBR,MAP,NODE,KBD,
     &             KTE,KTJ,KCM,IERR)
        IF(IERR.NE.0)RETURN

C       �쐬�����̈�Ɋ܂܂��v�f�ɗ̈�ԍ���ݒ�
        DO 30 K=1,NELM
          IF(IDM(K).EQ.IZ)THEN
            MA=IADRES(MTJ(K,1))
            MB=IADRES(MTJ(K,2))
            MC=IADRES(MTJ(K,3))
            MS=MA*MB*MC
            MP=(MB-MA)*(MC-MB)*(MA-MC) !�̈�������A�̈�O����
            IF((MS.NE.0).AND.(MP.GT.0))THEN
              IDM(K)=NB
            END IF
          END IF
   30   CONTINUE

*        CALL OUTEX(50,NELM,MTJ,JAC,IDM,PX,PY,KTJ)

   10 CONTINUE

C
C     OPENED BOUDNARIES
C
*      WRITE(*,*)'[ROUGH] GENERATION OF BREAK LINES'
*      CALL OUTEX(50,NELM,MTJ,JAC,IDM,PX,PY,KTJ)

      CALL ARRAYSET(IONL,KTJ,0)

      DO 200 I=1,NBK
        CALL READRES(JADRES,IAN,NBK,IBREAK,NBREAK,KBD,KTJ)
        NB=I
*        WRITE(*,*)'[ROUGH] OPEN BOUNDARY No = ',NB
        NP=IBREAK(NB)
        DO 210 J=1,NP
          LIST(J)=NBREAK(NB,J)
  210   CONTINUE
*        WRITE(*,*) '[ROUGH] LIST=',(LIST(J),J=1,NP)
*        WRITE(*,*) '[ROUGH] CALL BREAKLINE'
        CALL BREAKLINE(IZ,NP,NB,LIST,NTP,PX,PY,NER,NELM,ICHECK,NEX,
     &               IBEX,IBNO,NBK,IBREAK,NBREAK,MTJ,JAC,IDM,JADRES,IAN,
     &               DPP,ISTACK,IONL,KV,IBR,MAP,NODE,KBD,KTE,KTJ,KCM,PD,
     &               IERR)
        IF(IERR.NE.0)RETURN
*          CALL OUTEX(200,NELM,MTJ,JAC,IDM,PX,PY,KTJ)
  200 CONTINUE

C     SET POINT DENSITY ON OPEN BOUNDARY
      DO 250 I=1,NBK
        DO 260 J=1,IBREAK(I)
          N1=NBREAK(I,J)
          IF(PD(N1).NE.0.D0) GO TO 260
          SUP=0.D0
          IS=0
          CALL NEIGH(IERR,KTJ,KTE,KCM,MTJ,JAC,NER,N1,NEI,NEICNT)
          DO 270 K=1,NEICNT
            IELM=NEI(K)
            DO 280 L=1,3
              IV1=MTJ(IELM,L)
              IF(IV1.NE.N1)THEN
                SUP=SUP+PD(IV1)
                IS=IS+1
              END IF
  280       CONTINUE
  270     CONTINUE
          PD(N1)=SUP/IS
  260   CONTINUE
  250 CONTINUE

C
C     REMESH BY POINTS IN BOUNDARY BY THE DELAUNAY
C     TRIANGULATION
C
*      WRITE(*,*)'[ROUGH] REMESH BY POINTS IN BOUNDARY BY THE DELAUNAY'
      CALL ARRAYSET(IADRES,KTJ+3,0)

C     ���̓f�[�^�Ɋ܂܂������_���g�p���āA�O�p�����������Ȃ�
      DO 110 I=1,NIB
        NODE=NODE+1
        XA=PX(NODE)
        YA=PY(NODE)
*        WRITE(*,*)'[ROUGH] INNER POINT NUMBER = ',NODE
        CALL LOCATE(XA,YA,PX,PY,MTJ,JAC,NELM,KTE,KTJ,IT,IERR)
        IF(IERR.NE.0) RETURN
        IZ=IDM(IT)
        CALL DELAUN(IZ,NODE,NODE,NTP,PX,PY,NER,NELM,MTJ,
     &             JAC,IDM,IADRES,ISTACK,KTE,KTJ,KCM,IERR)
*        CALL OUTEX(NEX,NODE,NELM,MTJ,JAC,IDM,PX,PY,KTJ)
        IF(IERR.NE.0)RETURN
  110 CONTINUE

C     CHECK OF NODES GENERATED
      IF(NODE.NE.NTP)THEN
        IERR=201
        RETURN
C        WRITE(*,'('' ***ERROR IN SUBROUTINE ROUGH***'')')
C        WRITE(*,'('' ***INCORRECT NUMBER OF NODES***'')')
C        PAUSE
C        STOP
      END IF
C
      RETURN
      END
C     ------   SUBROUTINE BOUGEN   -------------------------
C
C     PURPOSE : MESH GENERATION FOR DOMAIN SURROUNDED
C               BY BOUNDARIES
C     LAST MODIFIED : 19 MAR 1990
C     [MEMO]
C     ICHECK��1:�����_�A0�F���쐬�_
C     IDT=1�Ȃ炱��BOUGEN�̒��ŐV���ɓ_������Ă�
C
      SUBROUTINE BOUGEN(IZ,NP,LIST,NTP,PX,PY,NER,NELM,ICHECK,
     &                 MTJ,JAC,IDM,IADRES,ISTACK,KV,IBR,MAP,
     &                 NODE,KBD,KTE,KTJ,KCM,IERR)
      IMPLICIT NONE
C     ==args==
      INTEGER KTJ,KTE,KCM,KBD
      INTEGER IZ,NP,LIST(KTJ),NTP,NELM
      INTEGER ICHECK(KTE),MTJ(KTE,3),JAC(KTE,3),IDM(KTE),IADRES(KTJ+3)
      INTEGER ISTACK(KTJ),KV(KTE),IBR(KTE),MAP(KTE),NODE,IERR,NER(KTJ+3)
      REAL*8 PX(KTJ+3),PY(KTJ+3)
C     ==tmp==
      INTEGER I,J,K
      INTEGER IZC(KBD),INP,IDT,IC,IS,KEY,IP,IQ,LOC,IV,NEI(KCM),NEICNT
      REAL*8 XS,YS

C     INITIALIZATION TEMP
      CALL ARRAYSET(IZC,KBD,0)
      INP=0
      IDT=0
      IC=0
      IS=0
      KEY=0
      IP=0
      IQ=0
      LOC=0
      IV=0
      XS=0.D0
      YS=0.D0

      CALL ARRAYSET(IADRES,KTJ+3,0)
C     ���E�_�������Z�b�g
      DO 20 I=1,NP
        IADRES(LIST(I))=I
   20 CONTINUE

C     CREATE THREE TRIANGLES BY INSERTION OF THE FIRST POINT
*      WRITE(*,*) '[BOUGEN] CREATE THREE TRIANGLES'
C     ���E�_���X�g�̍ŏ��ɂ���_��
      IS=LIST(1)
      XS=PX(IS)
      YS=PY(IS)
      IF(ICHECK(IS).NE.1)THEN
C       �V�K�쐬�_�̏ꍇ�A�����O�p�`��T��
        CALL LOCATE(XS,YS,PX,PY,MTJ,JAC,NELM,KTE,KTJ,LOC,IERR)
        IF(IERR.NE.0)RETURN
        IZ=IDM(LOC)
        IDT=1
      ELSE
C       �����_�̏ꍇ�A���ӗv�f���`�F�b�N
        CALL NEIGH(IERR,KTJ,KTE,KCM,MTJ,JAC,NER,IS,NEI,NEICNT)
        DO 70 I=1,NEICNT
          KEY=IDM(NEI(I))
          DO 75 J=1,IC+1
            IF(KEY.EQ.IZC(J))GO TO 70
   75     CONTINUE
          IC=IC+1
          IZC(IC)=KEY
   70   CONTINUE
      END IF

      IF(ICHECK(IS).NE.1)THEN
C       �V�K�쐬�_�̏ꍇ�A�f���[�j�[��������
        ICHECK(IS)=1
        INP=INP+1
*        WRITE(*,*) '[BOUGEN] CALL DELAUN'
        CALL DELAUN(IZ,IS,IS,NTP,PX,PY,NER,NELM,MTJ,JAC,
     &             IDM,IADRES,ISTACK,KTE,KTJ,KCM,IERR)
        IF(IERR.NE.0)RETURN
      END IF

C     LOOP OVER EACH POINT AND CONSTRUCT BOUNDARY EDGES
*      WRITE(*,*) '[BOUGEN] LOOP OVER EACH POINTS'
C     ���E�𐬂���������IP��IQ���ڑ����Ă��邩�`�F�b�N
      DO 30 I=1,NP
*        WRITE(*,*) '[BOUGEN] LINE No =',I
*        CALL OUTEX(50,NELM,MTJ,JAC,IDM,PX,PY,KTJ)
        IP=LIST(MOD(I,NP)+1)
        IQ=LIST(I)
        XS=PX(IP)
        YS=PY(IP)
*        WRITE(*,*) '[BOUGEN] IP=',IP,' IQ=',IQ
        IF(ICHECK(IP).NE.1.AND.IDT.EQ.0)THEN
*          WRITE(*,*) '[BOUGEN] CALL LOCATE'
          CALL LOCATE(XS,YS,PX,PY,MTJ,JAC,NELM,KTE,KTJ,LOC,IERR)
          IF(IERR.NE.0)RETURN
          IZ=IDM(LOC)
C          DO 90 J=1,IC
C            IF(IZ.EQ.IZC(J))THEN
              IDT=1
              GO TO 100
C            ELSE IF(J.EQ.IC)THEN
C              WRITE(*,610)IQ
C  610         FORMAT(' NODE NUMBER = ',I5)
C              WRITE(*,*)   'IP=',IP
C              WRITE(*,'('' ***ERROR IN SUBROUTINE BOUGEN***'')')
C              WRITE(*,'('' ***NEW POINT EXISTS OUTSIDE***'')')
C              STOP
C            END IF
C   90     CONTINUE
        ELSE IF(ICHECK(IP).EQ.1.AND.IDT.EQ.0)THEN
*          WRITE(*,*) '[BOUGEN] JNB(IP)=',JNB(IP)
          CALL NEIGH(IERR,KTJ,KTE,KCM,MTJ,JAC,NER,IP,NEI,NEICNT)
          DO 80 J=1,NEICNT
            KEY=IDM(NEI(J))
            DO 85 K=1,IC+1
              IF(KEY.EQ.IZC(K))GO TO 100
   85       CONTINUE
   80     CONTINUE
          IERR = 202
          RETURN
C          WRITE(*,620)IQ
C  620     FORMAT(' NODE NUMBER = ',I5)
C          WRITE(*,*)   'IP=',IP
C          WRITE(*,'('' ***ERROR IN SUBROUTINE BOUGEN***'')')
C          WRITE(*,'('' ***NEW POINT EXISTS OUTSIDE***'')')
C          PAUSE
C          STOP
        END IF

  100   CONTINUE
*        WRITE(*,*)' NODE NUMBER = ',IQ

        IF((ICHECK(IP).NE.1).AND.(I.NE.NP))THEN
C         IP�𔭐������A�ו������s��
          ICHECK(IP)=1
          INP=INP+1
*          WRITE(*,*) '[BOUGEN] CALL DELAUN'
          CALL DELAUN(IZ,IP,IP,NTP,PX,PY,NER,NELM,MTJ,
     &               JAC,IDM,IADRES,ISTACK,KTE,KTJ,KCM,IERR)
          IF(IERR.NE.0)RETURN
        END IF

C       IQ�̎��ӂ�IP������ꍇ�́A�ڑ�OK�Ƃ���
*        WRITE(*,*)'JNB(IQ) = ',JNB(IQ)
        CALL NEIGH(IERR,KTJ,KTE,KCM,MTJ,JAC,NER,IQ,NEI,NEICNT)
        DO 40 J=1,NEICNT
          DO 40 K=1,3
            IF(MTJ(NEI(J),K).EQ.IP)THEN
*              WRITE(*,*) '[BOUGEN] GOTO 30'
              GO TO 30
            END IF
   40   CONTINUE
C
C       SEARCH FOR TRIANGLES BETWEEN NEW DATA POINT AND
C       OLD DATA POINT
C
*        WRITE(*,*) '[BOUGEN] SEARCH FOR TRIANGLES'
C       IP��IQ�̊Ԃɂ���v�f��T��
        CALL SEARCH(IZ,IP,IQ,NER,NELM,MTJ,JAC,IDM,PX,PY,
     &             ISTACK,IV,KV,IADRES,IBR,KTE,KTJ,KCM,IERR)
        IF (IERR.NE.0) RETURN
C
C       REMOVE TRIANGLES BETWEEN NEW DATA POINT AND 
C       OLD DATA POINTAND SUBDIVIDE
C       THE POLYGON CONSISTED BY THEM
C
C       IP��IQ��ڑ�������悤�ɗv�f����������
*        WRITE(*,*) '[BOUGEN] CALL POLY'
        CALL POLY(IQ,IP,IV,KV,PX,PY,NELM,MTJ,JAC,IDM,NER,
     &            MAP,KTE,KTJ,KCM,IERR)
        IF (IERR.NE.0) RETURN
   30 CONTINUE
      NODE=NODE+INP
C
      RETURN
      END
C     -------   SUBROUTINE DELAUN   ------------------------
C
C     PURPOSE : DELAUNAY TRIANGULATION
C     LAST MODIFIED : 19 MAR 1990
C
C     IS - �ǉ��J�n�_No
C     IG - �ǉ��I���_No
      SUBROUTINE DELAUN(IZ,IS,IG,NTP,PX,PY,NER,NELM,MTJ,
     &                 JAC,IDM,IADRES,ISTACK,KTE,KTJ,KCM,IERR)
      IMPLICIT NONE
C     ==args==
      INTEGER KTJ,KTE,KCM
      INTEGER IZ,IS,IG,NTP,NELM
      INTEGER MTJ(KTE,3),JAC(KTE,3),IDM(KTE),IADRES(KTJ+3),ISTACK(KTJ)
      INTEGER IERR,NER(KTJ+3)
      REAL*8 PX(KTJ+3),PY(KTJ+3)
C     ==functions==
      INTEGER IPUSH
C     ==tmp==
      INTEGER I,ITOP,MAXSTK,IP,IA,IB,IC,IV1,IV2,IV3,MS,IDF,IL,IR,II,IT
      INTEGER IEDGE,IERA,IERB,IERL,ISWAP
      REAL*8 XP,YP

      ITOP=0
      MAXSTK=NTP
      IP=0
      IA=0
      IB=0
      IC=0
      IV1=0
      IV2=0
      IV3=0
      MS=0
      IDF=0
      IL=0
      IR=0
      II=0
      IT=0
      IEDGE=0
      IERA=0
      IERB=0
      IERL=0
      ISWAP=0
      XP=0.D0
      YP=0.D0

      DO 100 I=IS,IG
        IP=I
        XP=PX(IP)
        YP=PY(IP)

C       �_P���܂ގO�p�`��T���A�v�fNo��IT�ɃZ�b�g
        CALL LOCATE(XP,YP,PX,PY,MTJ,JAC,NELM,KTE,KTJ,IT,IERR)
        IF(IERR.NE.0)RETURN

C       �O�p�`IT��_P��3�ɕ����B
        IA=JAC(IT,1)
        IB=JAC(IT,2)
        IC=JAC(IT,3)
        IV1=MTJ(IT,1)
        IV2=MTJ(IT,2)
        IV3=MTJ(IT,3)
        MTJ(IT,1)=IP
        MTJ(IT,2)=IV1
        MTJ(IT,3)=IV2
        JAC(IT,1)=NELM+2
        JAC(IT,2)=IA
        JAC(IT,3)=NELM+1

        NELM=NELM+1
        IDM(NELM)=IZ
        MTJ(NELM,1)=IP
        MTJ(NELM,2)=IV2
        MTJ(NELM,3)=IV3
        JAC(NELM,1)=IT
        JAC(NELM,2)=IB
        JAC(NELM,3)=NELM+1
        NELM=NELM+1
        IDM(NELM)=IZ
        MTJ(NELM,1)=IP
        MTJ(NELM,2)=IV3
        MTJ(NELM,3)=IV1
        JAC(NELM,1)=NELM-1
        JAC(NELM,2)=IC
        JAC(NELM,3)=IT

C       ���ӗv�f�����X�V
*        WRITE(*,*) '[DELAUN] ���ӗv�f�����X�V'
        CALL INCR(IV1,NELM,NER,KTJ,KCM,IERR)
        IF (IERR.NE.0) RETURN
        CALL INCR(IV2,NELM-1,NER,KTJ,KCM,IERR)
        IF (IERR.NE.0) RETURN
        IF(IV3.LE.KTJ)THEN
          IF(NER(IV3).EQ.IT)THEN
            NER(IV3)=NELM-1
          ENDIF
        END IF
        CALL INCR(IV3,NELM,NER,KTJ,KCM,IERR)
        IF (IERR.NE.0) RETURN
        NER(IP)=IT

C       SWAP�����Ώۂ̗v�f�����X�g�A�b�v
*        WRITE(*,*) '[DELAUN] SWAP�����Ώۂ̗v�f�����X�g�A�b�v'
        IF(IA.NE.0)THEN
          MS=IADRES(IV1)*IADRES(IV2)
          IDF=IABS(IADRES(IV1)-IADRES(IV2))
          IF((IDM(IA).EQ.IZ).AND.((MS.EQ.0).OR.
     &                           (IDF.NE.1)))THEN
            ITOP=ITOP+1
            ISTACK(ITOP)=IPUSH(IT,MAXSTK,ITOP,ISTACK,KTJ,IERR)
            IF(IERR.NE.0)RETURN
          END IF
        END IF
        IF(IB.NE.0)THEN
          CALL EDGE(IB,IT,JAC,KTE,IEDGE,IERR)
          IF (IERR.NE.0) RETURN
          JAC(IB,IEDGE)=NELM-1
          MS=IADRES(IV2)*IADRES(IV3)
          IDF=IABS(IADRES(IV2)-IADRES(IV3))
          IF((IDM(IB).EQ.IZ).AND.((MS.EQ.0).OR.
     &                           (IDF.NE.1)))THEN
            ITOP=ITOP+1
            ISTACK(ITOP)=IPUSH(NELM-1,MAXSTK,ITOP,
     &                                     ISTACK,KTJ,IERR)
            IF(IERR.NE.0)RETURN
          END IF
        END IF
        IF(IC.NE.0)THEN
          CALL EDGE(IC,IT,JAC,KTE,IEDGE,IERR)
          IF (IERR.NE.0) RETURN
          JAC(IC,IEDGE)=NELM
          MS=IADRES(IV3)*IADRES(IV1)
          IDF=IABS(IADRES(IV3)-IADRES(IV1))
          IF((IDM(IC).EQ.IZ).AND.((MS.EQ.0).OR.
     &                           (IDF.NE.1)))THEN
            ITOP=ITOP+1
            ISTACK(ITOP)=IPUSH(NELM,MAXSTK,ITOP,ISTACK,KTJ,IERR)
            IF(IERR.NE.0)RETURN
          END IF
        END IF

C       SWAP����
C        WRITE(*,*) '[DELAUN] SWAP'
   50   IF(ITOP.GT.0)THEN
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
C         �v�fIL��IR�̋��L�ӂ�SWAP�Ώۂł��邩���f
          CALL SWAP(PX(IV1),PY(IV1),PX(IV2),PY(IV2),PX(IV3),
     &              PY(IV3),XP,YP,ISWAP)
          IF(ISWAP.EQ.1)THEN
C           SWAP�Ώۂ̏ꍇ�A�Ίp������ւ�
            IA=JAC(IR,IERA)
            IB=JAC(IR,IERB)
            IC=JAC(IL,3)

            MTJ(IL,3)=IV3
            JAC(IL,2)=IA
            JAC(IL,3)=IR

            MTJ(IR,1)=IP
            MTJ(IR,2)=IV3
            MTJ(IR,3)=IV1
            JAC(IR,1)=IL
            JAC(IR,2)=IB
            JAC(IR,3)=IC

            CALL DECR(IV1,IL,MTJ,JAC,NER,KTJ,KTE,KCM,IERR)
            IF (IERR.NE.0) RETURN
            CALL DECR(IV2,IR,MTJ,JAC,NER,KTJ,KTE,KCM,IERR)
            IF (IERR.NE.0) RETURN
            CALL INCR(IP,IR,NER,KTJ,KCM,IERR)
            IF (IERR.NE.0) RETURN
            CALL INCR(IV3,IL,NER,KTJ,KCM,IERR)
            IF (IERR.NE.0) RETURN

C           �V���ɍ��ꂽ�ӂ�SWAP�����Ώۂɓ����
            IF(IA.NE.0)THEN
              CALL EDGE(IA,IR,JAC,KTE,IEDGE,IERR)
              IF (IERR.NE.0) RETURN
              JAC(IA,IEDGE)=IL
              MS=IADRES(IV2)*IADRES(IV3)
              IDF=IABS(IADRES(IV2)-IADRES(IV3))
              IF((IDM(IA).EQ.IZ).AND.((MS.EQ.0).OR.
     &                               (IDF.NE.1)))THEN
                ITOP=ITOP+1
                ISTACK(ITOP)=IPUSH(IL,MAXSTK,ITOP,
     &                                       ISTACK,KTJ,IERR)
                IF(IERR.NE.0)RETURN
              END IF
            END IF
            IF(IB.NE.0)THEN
              MS=IADRES(IV3)*IADRES(IV1)
              IDF=IABS(IADRES(IV3)-IADRES(IV1))
              IF((IDM(IB).EQ.IZ).AND.((MS.EQ.0).OR.
     &                               (IDF.NE.1)))THEN
                ITOP=ITOP+1
                ISTACK(ITOP)=IPUSH(IR,MAXSTK,ITOP,
     &                                       ISTACK,KTJ,IERR)
                IF(IERR.NE.0)RETURN
              END IF
            END IF
            IF(IC.NE.0)THEN
              CALL EDGE(IC,IL,JAC,KTE,IEDGE,IERR)
              IF (IERR.NE.0) RETURN
              JAC(IC,IEDGE)=IR
            END IF
          END IF
          GO TO 50
        END IF
  100 CONTINUE
C
      RETURN
      END
C     -------   SUBROUTINE BDELAUN   ------------------------
C
C     PURPOSE : �O�p�`�̕ӏ�ɓ_��ǉ�����p
C     MEMO:
C     �E���E����������\��������Ƃ��͂�����g��
C     �EIADRES�ł͂Ȃ�JADRES���g�p���Ă���
C
      SUBROUTINE BDELAUN(IZ,IS,IG,NTP,PX,PY,NER,NELM,MTJ,
     &                 JAC,IDM,JADRES,IAN,ISTACK,KBD,KTE,KTJ,KCM,IERR)
      IMPLICIT NONE
C     ==args==
      INTEGER KTJ,KTE,KCM,KBD
      INTEGER IZ,IS,IG,NTP,NELM,MTJ(KTE,3)
      INTEGER JAC(KTE,3),IDM(KTE),JADRES(KTJ+3,KBD),IAN(KTJ+3)
      INTEGER ISTACK(KTJ),IERR,NER(KTJ+3)
      REAL*8 PX(KTJ+3),PY(KTJ+3)
C     ==functions==
      INTEGER IPUSH
C     ==tmp==
      INTEGER I,J,K
      INTEGER ITOP,MAXSTK,IP,IT,IA,IB,IC,IV1,IV2,IV3,II,IDF
      INTEGER KIDF,IL,IR,ISWAP,IERA,IERB,IERL,IEDGE,MS
      REAL*8 XP,YP

      ITOP=0
      MAXSTK=NTP
      IP=0
      IT=0
      IA=0
      IB=0
      IC=0
      IV1=0
      IV2=0
      IV3=0
      II=0
      IDF=0
      KIDF=0
      IL=0
      IR=0
      ISWAP=0
      IERA=0
      IERB=0
      IERL=0
      IEDGE=0
      MS=0
      XP=0.D0
      YP=0.D0

      DO 100 I=IS,IG
        IP=I
        XP=PX(IP)
        YP=PY(IP)

*        WRITE(*,*) '[BDELAUN] IP=',I
*        WRITE(*,*) '[BDELAUN] CALL LOCATE'
        CALL LOCATE(XP,YP,PX,PY,MTJ,JAC,NELM,KTE,KTJ,IT,IERR)
        IF(IERR.NE.0)RETURN

        IA=JAC(IT,1)
        IB=JAC(IT,2)
        IC=JAC(IT,3)
        IV1=MTJ(IT,1)
        IV2=MTJ(IT,2)
        IV3=MTJ(IT,3)
        MTJ(IT,1)=IP
        MTJ(IT,2)=IV1
        MTJ(IT,3)=IV2
        JAC(IT,1)=NELM+2
        JAC(IT,2)=IA
        JAC(IT,3)=NELM+1

        NELM=NELM+1
        IDM(NELM)=IZ
        MTJ(NELM,1)=IP
        MTJ(NELM,2)=IV2
        MTJ(NELM,3)=IV3
        JAC(NELM,1)=IT
        JAC(NELM,2)=IB
        JAC(NELM,3)=NELM+1
        NELM=NELM+1
        IDM(NELM)=IZ
        MTJ(NELM,1)=IP
        MTJ(NELM,2)=IV3
        MTJ(NELM,3)=IV1
        JAC(NELM,1)=NELM-1
        JAC(NELM,2)=IC
        JAC(NELM,3)=IT

        CALL INCR(IV1,NELM,NER,KTJ,KCM,IERR)
        IF (IERR.NE.0) RETURN
        CALL INCR(IV2,NELM-1,NER,KTJ,KCM,IERR)
        IF (IERR.NE.0) RETURN
        IF(IV3.LE.KTJ)THEN
          IF(NER(IV3).EQ.IT)THEN
            NER(IV3)=NELM-1
          ENDIF
        END IF
        CALL INCR(IV3,NELM,NER,KTJ,KCM,IERR)
        IF (IERR.NE.0) RETURN
        NER(IP)=IT

        IF(IA.NE.0)THEN
          MS=JADRES(IV1,1)*JADRES(IV2,1)
          IDF=10
          DO 110 J=1,IAN(IV1)
            DO 110 K=1,IAN(IV2)
              KIDF=IABS(JADRES(IV1,J)-JADRES(IV2,K))
              IF(KIDF.LT.IDF)IDF=KIDF
  110     CONTINUE
          IF((IDM(IB).EQ.IZ).AND.((MS.EQ.0).OR.(IDF.NE.1)))THEN
            ITOP=ITOP+1
            ISTACK(ITOP)=IPUSH(IT,MAXSTK,ITOP,ISTACK,KTJ,IERR)
            IF(IERR.NE.0)RETURN
          END IF
        END IF
        IF(IB.NE.0)THEN
          CALL EDGE(IB,IT,JAC,KTE,IEDGE,IERR)
          IF (IERR.NE.0) RETURN
          JAC(IB,IEDGE)=NELM-1
          MS=JADRES(IV2,1)*JADRES(IV3,1)
          IDF=10
          DO 120 J=1,IAN(IV2)
            DO 120 K=1,IAN(IV3)
              KIDF=IABS(JADRES(IV2,J)-JADRES(IV3,K))
              IF(KIDF.LT.IDF)IDF=KIDF
  120     CONTINUE
          IF((IDM(IB).EQ.IZ).AND.((MS.EQ.0).OR.(IDF.NE.1)))THEN
            ITOP=ITOP+1
            ISTACK(ITOP)=IPUSH(NELM-1,MAXSTK,ITOP,ISTACK,KTJ,IERR)
            IF(IERR.NE.0)RETURN
          END IF
        END IF
        IF(IC.NE.0)THEN
          CALL EDGE(IC,IT,JAC,KTE,IEDGE,IERR)
          IF (IERR.NE.0) RETURN
          JAC(IC,IEDGE)=NELM
          MS=JADRES(IV3,1)*JADRES(IV1,1)
          IDF=10
          DO 130 J=1,IAN(IV3)
            DO 130 K=1,IAN(IV1)
              KIDF=IABS(JADRES(IV3,J)-JADRES(IV1,K))
              IF(KIDF.LT.IDF)IDF=KIDF
  130     CONTINUE
          IF((IDM(IC).EQ.IZ).AND.((MS.EQ.0).OR.(IDF.NE.1)))THEN
            ITOP=ITOP+1
            ISTACK(ITOP)=IPUSH(NELM,MAXSTK,ITOP,ISTACK,KTJ,IERR)
            IF(IERR.NE.0)RETURN
          END IF
        END IF

   50   IF(ITOP.GT.0)THEN
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
*          WRITE(*,*) '[BDELAUN] CALL SWAP'
          CALL SWAP(PX(IV1),PY(IV1),PX(IV2),PY(IV2),PX(IV3),
     &              PY(IV3),XP,YP,ISWAP)
*          WRITE(*,*) '[BDELAUN] ISWAP=',ISWAP
          IF(ISWAP.EQ.1)THEN

            IA=JAC(IR,IERA)
            IB=JAC(IR,IERB)
            IC=JAC(IL,3)

            MTJ(IL,3)=IV3
            JAC(IL,2)=IA
            JAC(IL,3)=IR

            MTJ(IR,1)=IP
            MTJ(IR,2)=IV3
            MTJ(IR,3)=IV1
            JAC(IR,1)=IL
            JAC(IR,2)=IB
            JAC(IR,3)=IC

            CALL DECR(IV1,IL,MTJ,JAC,NER,KTJ,KTE,KCM,IERR)
            IF (IERR.NE.0) RETURN
            CALL DECR(IV2,IR,MTJ,JAC,NER,KTJ,KTE,KCM,IERR)
            IF (IERR.NE.0) RETURN
            CALL INCR(IP,IR,NER,KTJ,KCM,IERR)
            IF (IERR.NE.0) RETURN
            CALL INCR(IV3,IL,NER,KTJ,KCM,IERR)
            IF (IERR.NE.0) RETURN

*            WRITE(*,*) '[BDELAUN] IA=',IA
            IF(IA.NE.0)THEN
*              WRITE(*,*) '[BDELAUN] CALL EDGE'
              CALL EDGE(IA,IR,JAC,KTE,IEDGE,IERR)
              IF (IERR.NE.0) RETURN
              JAC(IA,IEDGE)=IL
              MS=JADRES(IV2,1)*JADRES(IV3,1)
              IDF=10
*              WRITE(*,*) '[BDELAUN] 140 LOOP - ',IAN(IV2)
              DO 140 J=1,IAN(IV2)
                DO 140 K=1,IAN(IV3)
                  KIDF=IABS(JADRES(IV2,J)-JADRES(IV3,K))
                  IF(KIDF.LT.IDF)IDF=KIDF
  140         CONTINUE
              IF((IDM(IA).EQ.IZ).AND.((MS.EQ.0).OR.(IDF.NE.1)))THEN
                ITOP=ITOP+1
                ISTACK(ITOP)=IPUSH(IL,MAXSTK,ITOP,ISTACK,KTJ,IERR)
                IF(IERR.NE.0)RETURN
              END IF
            END IF
            IF(IB.NE.0)THEN
              MS=JADRES(IV3,1)*JADRES(IV1,1)
              IDF=10
*              WRITE(*,*) '[BDELAUN] 150 LOOP - ',IAN(IV3)
              DO 150 J=1,IAN(IV3)
                DO 150 K=1,IAN(IV1)
                  KIDF=IABS(JADRES(IV3,J)-JADRES(IV1,K))
                  IF(KIDF.LT.IDF)IDF=KIDF
  150         CONTINUE
              IF((IDM(IB).EQ.IZ).AND.((MS.EQ.0).OR.(IDF.NE.1)))THEN
                ITOP=ITOP+1
                ISTACK(ITOP)=IPUSH(IR,MAXSTK,ITOP,ISTACK,KTJ,IERR)
                IF(IERR.NE.0)RETURN
              END IF
            END IF
            IF(IC.NE.0)THEN
*              WRITE(*,*) '[BDELAUN] CALL EDGE'
              CALL EDGE(IC,IL,JAC,KTE,IEDGE,IERR)
              IF (IERR.NE.0) RETURN
              JAC(IC,IEDGE)=IR
            END IF
          END IF
          GO TO 50
        END IF
  100 CONTINUE

      RETURN
      END
C     ------   SUBROUTINE SEARCH   -------------------------
C
C     PURPOSE : SEARCH FOR TRIANGLES TO BE MODIFIED
C
C     [DEFINITION]
C     IZ - �̈�ԍ�
C     IP - I+1�_�̐ߓ_�ԍ�
C     IQ - I�_�̐ߓ_�ԍ�
C     JNB - �e�ߓ_�̎��ӗv�f���
C     NEI - �e�ߓ_�̎��ӗv�f�ԍ�
C     NELM - �v�f��
C     IV - i�_��i+1�_�Ԃ̗v�f��
C     KV - i�_��i+1�_�Ԃ̗v�f�ԍ�
C     IBR - ����ԍ�
      SUBROUTINE SEARCH(IZ,IP,IQ,NER,NELM,MTJ,JAC,IDM,PX,PY,
     &                 ISTACK,IV,KV,IADRES,IBR,KTE,KTJ,KCM,IERR)
      IMPLICIT NONE
C     ==args==
      INTEGER KTJ,KTE,KCM
      INTEGER IZ,IP,IQ,NELM,NER(KTJ+3)
      INTEGER MTJ(KTE,3),JAC(KTE,3),IDM(KTE),ISTACK(KTJ),IV,KV(KTE)
      INTEGER IADRES(KTJ+3),IBR(KTE),IERR
      REAL*8 PX(KTJ+3),PY(KTJ+3)
C     ==functions==
      INTEGER IVERT
C     ==tmp==
      INTEGER I,J,II,II1,II2
      INTEGER MSTK,NBR,JA,JB,IA,IB,IDF,MS,IELM,JELM,KELM,IMIN,JR
      INTEGER NEI(KCM),NEICNT,LASTELM
      REAL*8 DMAXX,DMINX,DMAXY,DMINY,DLMAXX,DLMINX,DLMAXY,DLMINY
      REAL*8 TH0,TH1,TH2,THETA
      LOGICAL FLAG

C     INITIALIZATION
      IV=0
      CALL ARRAYSET(KV,NELM,0)
      CALL ARRAYSET(IBR,NELM,0)
      CALL ARRAYSET(ISTACK,KTJ,0)

      MSTK=0
      NBR=0
      JA=0
      JB=0
      IA=0
      IB=0
      IDF=0
      MS=0
      IELM=0
      JELM=0
      KELM=0
      IMIN=0
      JR=0
      NEICNT=0

      DLMAXX=DMAX1(PX(IP),PX(IQ))
      DLMINX=DMIN1(PX(IP),PX(IQ))
      DLMAXY=DMAX1(PY(IP),PY(IQ))
      DLMINY=DMIN1(PY(IP),PY(IQ))

*      WRITE(*,*) '[SEARCH] IZ=',IZ,' IP=',IP,' IQ=',IQ

C     �_IQ�̎��ӗv�f�ɕ���ԍ����Z�b�g
      CALL NEIGH(IERR,KTJ,KTE,KCM,MTJ,JAC,NER,IQ,NEI,NEICNT)
*      WRITE(*,*) '[SEARCH] IQ=',IQ,'-',(NEI(I),I=1,NEICNT)
      DO 30 I=1,NEICNT
        NBR=NBR+1
        IBR(NEI(I))=NBR
   30 CONTINUE

      FLAG=.FALSE.
      DO 40 I=1,NEICNT
        IELM=NEI(I)
*        WRITE(*,*) '[SEARCH] IELM=',IELM,':',(MTJ(IELM,J),J=1,3)
        J=IVERT(IELM,IQ,MTJ,KTE,IERR)
        IF(IERR.NE.0)RETURN
        JA=MOD(J,3)+1
        JB=MOD(JA,3)+1
        IA=MTJ(IELM,JA)
        IB=MTJ(IELM,JB)
        IDF=IABS(IADRES(IA)-IADRES(IB))
        MS=IADRES(IA)*IADRES(IB)
C       �_IA�AIB�����E��̓_�ŘA���ł���ꍇ�͕��򂵂Ȃ�
        IF((IDF.EQ.1).AND.(MS.NE.0))GO TO 40
        JELM=JAC(IELM,JA)
        IF(JELM.EQ.0)GO TO 40
*        WRITE(*,*) '[SEARCH] JELM=',JELM,':',(MTJ(JELM,J),J=1,3)
        NBR=NBR+1
C       ��������쐬
        IBR(JELM)=NBR
C       �_I�̎��ӂɂ���v�fIELM�̗אڗv�fJELM�̐ߓ_���_I+1�ł��邩
        IF(MTJ(JELM,1).EQ.IP 
     &      .OR. MTJ(JELM,2).EQ.IP 
     &      .OR. MTJ(JELM,3).EQ.IP)THEN
          LASTELM=JELM  
          FLAG=.TRUE.
        ENDIF
C        IF(MTJ(JELM,1).EQ.IP)GO TO 80
C        IF(MTJ(JELM,2).EQ.IP)GO TO 80
C        IF(MTJ(JELM,3).EQ.IP)GO TO 80
        MSTK=MSTK+1
        ISTACK(MSTK)=JELM
   40 CONTINUE

      IF(FLAG.EQV..TRUE.)GOTO 80

*      WRITE(*,*) '[SEARCH] IBR='
*      DO I=1,NELM
*        IF(IBR(I).NE.0)THEN
*          WRITE(*,*) I,'-',(MTJ(I,J),J=1,3)
*        ENDIF
*      ENDDO

C     �_I�̎��ӗv�f�̗אڗv�f�œ_I+1��������Ȃ������ꍇ�A����ɂ���
C     �אڗv�f��T��
*      WRITE(*,*) '[SEARCH] SEARCH NEIGHBOR'
   50 CONTINUE
      IF(MSTK.EQ.0)THEN
        IERR = 203
        RETURN
C        WRITE(*,'('' ***ERROR IN SUBROUTINE SEARCH***'')')
C        WRITE(*,'('' ***ISTACK IS EMPTY***'')')
C        PAUSE
C        STOP
      END IF
      IELM=ISTACK(1)
      MSTK=MSTK-1
      DO 60 I=1,MSTK
        ISTACK(I)=ISTACK(I+1)
   60 CONTINUE
      ISTACK(MSTK+1)=0
*      WRITE(*,*) '[SEARCH] IELM=',IELM,':',(MTJ(IELM,J),J=1,3)
      DO 70 J=1,3
        JA=MOD(J,3)+1
        IA=MTJ(IELM,J)
        IB=MTJ(IELM,JA)
*        WRITE(*,*) '[SEARCH] IA=',IA,', IB=',IB
        IDF=IABS(IADRES(IA)-IADRES(IB))
        MS=IADRES(IA)*IADRES(IB)
        IF((IDF.EQ.1).AND.(MS.NE.0))GO TO 70
        JELM=JAC(IELM,J)
*        WRITE(*,*) '[SEARCH] JELM=',JELM,':',(MTJ(JELM,K),K=1,3)
        IF(JELM.EQ.0)GO TO 70
        IF(IBR(JELM).NE.0)GO TO 70
        NBR=NBR+1
C       ��������쐬
        IBR(JELM)=NBR
        IF(MTJ(JELM,1).EQ.IP 
     &      .OR. MTJ(JELM,2).EQ.IP 
     &      .OR. MTJ(JELM,3).EQ.IP)THEN
          IF(MTJ(JELM,1).EQ.IP) II=1
          IF(MTJ(JELM,2).EQ.IP) II=2
          IF(MTJ(JELM,3).EQ.IP) II=3
          II1=MTJ(JELM,MOD(II,3)+1)
          II2=MTJ(JELM,MOD(MOD(II,3)+1,3)+1)

C== add by watanabe 07/08/29
C         �Y���O�p�`�Ɛ���PQ�̌����`�F�b�N 
          TH0 = THETA(PX(IP),PY(IP)
     &        ,PX(II1),PY(II1),PX(II2),PY(II2))
          TH1 = THETA(PX(IP),PY(IP)
     &        ,PX(IQ),PY(IQ),PX(II1),PY(II1))
          TH2 = THETA(PX(IP),PY(IP)
     &        ,PX(IQ),PY(IQ),PX(II2),PY(II2))

*          WRITE(*,*) '[SEARCH] ',IP,II1,II2,IQ
*          WRITE(*,*) '[SEARCH] TH0=',TH0,' TH1=',TH1,' TH2=',TH2
          IF(TH1.LT.1E-10 .OR. TH2.LT.1E-10)THEN
C           �ӂƏd�Ȃ�ꍇ
            GO TO 70
          ELSEIF(TH1.GT.TH0 .OR. TH2.GT.TH0)THEN
C           �������Ă��Ȃ��ꍇ
            GO TO 70
          ELSE
C           �������Ă���ꍇ
            LASTELM=JELM  
            GO TO 80
          ENDIF
C==
        ENDIF
C        IF(MTJ(JELM,1).EQ.IP)GO TO 80
C        IF(MTJ(JELM,2).EQ.IP)GO TO 80
C        IF(MTJ(JELM,3).EQ.IP)GO TO 80
        MSTK=MSTK+1
        ISTACK(MSTK)=JELM
   70 CONTINUE
      GO TO 50


C     POINT HAS BEEN FOUND
C     �ߓ_I+1�𔭌������ꍇ�A�_I+1����_I�܂ł̍ŒZ�o�H�𕪊����
C     ���p���Čv�Z����
   80 CONTINUE
*      WRITE(*,*) '[SEARCH] IBR='
*      DO I=1,NELM
*        IF(IBR(I).NE.0)THEN
*          WRITE(*,*) I,'-',IBR(I),'-',(MTJ(I,J),J=1,3)
*        ENDIF
*      ENDDO
*      WRITE(*,*) '[SEARCH] POINT HAS BEEN FOUND'

      JELM=LASTELM

   85 CONTINUE
      IV=IV+1
      KV(IV)=JELM
*      WRITE(*,*)'[SEARCH] JELM=',JELM
      IF(IBR(JELM).LE.NEICNT)GO TO 100
      IMIN=KTE+1
*      WRITE(*,*)'[SEARCH] IMIN=',IMIN
      DO 90 J=1,3
        JR=JAC(JELM,J)
        IF(JR.EQ.0)GO TO 90
        IF(IBR(JR).EQ.0)GO TO 90
        JA=MOD(J,3)+1
        IA=MTJ(JELM,J)
        IB=MTJ(JELM,JA)
        IDF=IABS(IADRES(IA)-IADRES(IB))
        MS=IADRES(IA)*IADRES(IB)
*        WRITE(*,*)'[SEARCH] IA=',IA,' IB=',IB,' IDF=',IDF,',MS=',MS
        IF((IDF.EQ.1).AND.(MS.NE.0))GO TO 90
*        WRITE(*,*)'[SEARCH] IBR(JR)=',IBR(JR),' IMIN=',IMIN
        IF(IBR(JR).LT.IMIN)THEN

C==       �O�p�`������PQ�ƌ������Ă��邩�`�F�b�N���鏈����ǉ��B07/04/12
C==       ����MAX-MIN�ł������f���Ă��Ȃ��B�ӂ̌����v�Z��ǉ����ׂ�
          DMAXX=DMAX1(PX(MTJ(JR,1)),PX(MTJ(JR,2)),PX(MTJ(JR,3)))
          DMINX=DMIN1(PX(MTJ(JR,1)),PX(MTJ(JR,2)),PX(MTJ(JR,3)))
          DMAXY=DMAX1(PY(MTJ(JR,1)),PY(MTJ(JR,2)),PY(MTJ(JR,3)))
          DMINY=DMIN1(PY(MTJ(JR,1)),PY(MTJ(JR,2)),PY(MTJ(JR,3)))
          
          IF((DMAXX.LT.DLMINX .OR. DLMAXX.LT.DMINX).OR.
     &        (DMAXY.LT.DLMINY .OR. DLMAXY.LT.DMINY)) THEN
            GOTO 90
          ENDIF
          KELM=JR
          IMIN=IBR(JR)
        END IF
   90 CONTINUE
      IF(IMIN.EQ.KTE+1)THEN
        IERR=251
        RETURN
C        WRITE(*,'('' ***ERROR IN SUBROUTINE SEARCH***'')')
C        WRITE(*,'('' ***BRANCH IS CLOSED***'')')
C        PAUSE
C        STOP
      END IF
      IF(JELM.EQ.KELM)THEN
        IERR=252
        RETURN
      ENDIF
      JELM=KELM
      GO TO 85

  100 CONTINUE

      RETURN
      END
C     ------   SUBROUTINE POLY   ---------------------------
C
C     PURPOSE : SUBDIVIDE GIVEN POLYGON INTO FINE ELEMENT
C     ARGUMENTS:
C       IQ - �_i
C       IP - �_i+1
C       IV - �_i�Ɠ_i+1�Ƃ̊ԂɈʒu����v�f��
C       KV - �_i�Ɠ_i+1�Ƃ̊ԂɈʒu����v�f�ԍ�
      SUBROUTINE POLY(IQ,IP,IV,KV,PX,PY,NELM,MTJ,JAC,IDM,
     &               NER,MAP,KTE,KTJ,KCM,IERR)
      IMPLICIT NONE
C     ==args==
      INTEGER KTJ,KTE,KCM
      INTEGER IQ,IP,IV,KV(KTE),NELM,MTJ(KTE,3),JAC(KTE,3),NER(KTJ+3)
      INTEGER IDM(KTE),MAP(KTE),IERR
      REAL*8 PX(KTJ+3),PY(KTJ+3)
C     ==tmp==
      INTEGER I,J,K,L
      INTEGER LTE
      PARAMETER(LTE=100)
      INTEGER NSRA(LTE+2),NSRB(LTE+2),IENA(LTE,3),IENB(LTE,3)
      INTEGER IHEN(2*LTE+1,2),JHEN(2*LTE+1),IAD(2*LTE+1),JEEA(LTE,3)
      INTEGER JEEB(LTE,3),JSTACK(LTE+2),NEI(KCM),NEICNT
      INTEGER IA,JA,IR,JR,IEDGE,IV1,IV2,LTJ,LHN,NPA,NPB,IVX,IELM,NTA
      INTEGER NTB,JELM,IPS,IPG,KELM,IVA,IVB,NLIST(LTE),NLCNT

      CALL ARRAYSET(NSRA,LTE+2,0)
      CALL ARRAYSET(NSRB,LTE+2,0)
      CALL ARRAYSET(IENA,LTE*3,0)
      CALL ARRAYSET(IENB,LTE*3,0)
      CALL ARRAYSET(IHEN,(2*LTE+1)*2,0)
      CALL ARRAYSET(JHEN,2*LTE+1,0)
      CALL ARRAYSET(IAD,2*LTE+1,0)
      CALL ARRAYSET(JEEA,LTE*3,0)
      CALL ARRAYSET(JEEB,LTE*3,0)
      CALL ARRAYSET(JSTACK,LTE+2,0)
      CALL ARRAYSET(NLIST,LTE,0)
      IA=0
      JA=0
      IR=0
      JR=0
      IEDGE=0
      IV1=0
      IV2=0
      LTJ=0
      LHN=0
      NPA=0
      NPB=0
      IVX=0
      IELM=0
      JELM=0
      KELM=0
      NTA=0
      NTB=0
      IPS=0
      IPG=0
      IVA=0
      IVB=0
      NLCNT=0

*      WRITE(*,*) '[POLY] IV=',IV
      IF(IV.EQ.2)THEN
C       �_i�Ɠ_i+1�Ɉʒu����v�f����2�̏ꍇ�A�Ίp�������ւ�
        CALL EDGE(KV(1),KV(2),JAC,KTE,IA,IERR)
        IF (IERR.NE.0) RETURN
        CALL EDGE(KV(2),KV(1),JAC,KTE,JA,IERR)
        IF (IERR.NE.0) RETURN

        IR=JAC(KV(1),MOD(IA,3)+1)
        JR=JAC(KV(2),MOD(JA,3)+1)

        MTJ(KV(1),MOD(IA,3)+1)=IQ
        JAC(KV(1),IA)=JR
        JAC(KV(1),MOD(IA,3)+1)=KV(2)

        MTJ(KV(2),MOD(JA,3)+1)=IP
        JAC(KV(2),JA)=IR
        JAC(KV(2),MOD(JA,3)+1)=KV(1)

        IF(IR.NE.0)THEN
          CALL EDGE(IR,KV(1),JAC,KTE,IEDGE,IERR)
          IF (IERR.NE.0) RETURN
          JAC(IR,IEDGE)=KV(2)
        END IF
        IF(JR.NE.0)THEN
          CALL EDGE(JR,KV(2),JAC,KTE,IEDGE,IERR)
          IF (IERR.NE.0) RETURN
          JAC(JR,IEDGE)=KV(1)
        END IF

        IV1=MTJ(KV(1),IA)
        IV2=MTJ(KV(2),JA)
        CALL DECR(IV1,KV(2),MTJ,JAC,NER,KTJ,KTE,KCM,IERR)
        IF (IERR.NE.0) RETURN
        CALL DECR(IV2,KV(1),MTJ,JAC,NER,KTJ,KTE,KCM,IERR)
        IF (IERR.NE.0) RETURN
        CALL INCR(IQ,KV(1),NER,KTJ,KCM,IERR)
        IF (IERR.NE.0) RETURN
        CALL INCR(IP,KV(2),NER,KTJ,KCM,IERR)
        IF (IERR.NE.0) RETURN
      ELSE
C       �_i�Ɠ_i+1�̊Ԃɗv�f��3�ȏ゠��ꍇ�A�֌W����O�p�`�Q��
C       �_i�Ai+1������
C       ���ŗ̈�A,B��2�ɕ����A���ꂼ��̗̈�ŎO�p�`�������������B
C       
*        CALL OUTEX(200,NELM,MTJ,JAC,IDM,PX,PY,KTJ)

        LTJ=LTE+2
        LHN=2*LTE+1
        NPA=0
        NPB=0
        CALL ARRAYSET(NSRA, LTJ, 0)
        CALL ARRAYSET(NSRB, LTJ, 0)
        CALL ARRAYSET(MAP, NELM, 0)
        DO 30 I=1,IV
          MAP(KV(I))=1
   30   CONTINUE

        NLCNT=0
        DO I=1,IV
          DO 31 J=1,3
            IELM=JAC(KV(I),J)
            IF(IELM.EQ.0) GOTO 31
            DO 32 K=1,NLCNT
              IF(NLIST(K).EQ.IELM)GOTO 31
   32       CONTINUE
            NLCNT=NLCNT+1
            NLIST(NLCNT)=IELM
   31     CONTINUE
        ENDDO

        DO 40 I=1,IV
          IELM=KV(I)
          DO 40 J=1,3
            IVX=MTJ(IELM,J)
            CALL DECR(IVX,IELM,MTJ,JAC,NER,KTJ,KTE,KCM,IERR)
            IF (IERR.NE.0) RETURN
   40   CONTINUE

C       PICK UP NODES WHICH SURROUND GIVEN POLYGON
*        WRITE(*,*) '[POLY] CALL PICK'
        CALL PICK(IQ,IP,IV,KV,MTJ,JAC,MAP,NPA,NPB,NSRA,NSRB,
     &           KTE,KTJ,LTJ,IERR)
        IF (IERR.NE.0) RETURN
        

C       DIVIDE GIVEN POLYGON INTO TWO DOMAINS
C       AND APPLY THE MODIFIED-DELAUNAY TO EACH DOMAIN
*        WRITE(*,*) '[POLY] CALL SUBDIV'
        CALL SUBDIV(NPA,NSRA,PX,PY,NTA,IENA,JEEA,IHEN,JHEN,
     &             IAD,JSTACK,KTJ,LTE,LTJ,LHN,IERR)
        IF (IERR.NE.0) RETURN
*        CALL CHECK(NELM,MTJ,JAC,KTE,IERR)
*        WRITE(*,*) '[POLY] CALL SUBDIV'
        CALL SUBDIV(NPB,NSRB,PX,PY,NTB,IENB,JEEB,IHEN,JHEN,
     &             IAD,JSTACK,KTJ,LTE,LTJ,LHN,IERR)
        IF (IERR.NE.0) RETURN
*        CALL CHECK(NELM,MTJ,JAC,KTE,IERR)

C       CHECK OF ELEMENTS SUBDIVIDED
        IF(IV.NE.NTA+NTB)THEN
          IERR = 205
          RETURN
C          WRITE(*,'('' ***ERROR IN SUBROUTINE POLY***'')')
C          WRITE(*,'('' ***INCORRECT NUMBER OF
C     &                             ELEMENTS FORMED***'')')
C          PAUSE
C          STOP
        END IF

C       MODIFICATION OF MTJ,JAC,JNB,NEI
*        WRITE(*,*) '[POLY] MODIFICATION OF MTJ,JAC'
*        WRITE(*,*) '[POLY] KV-',(KV(I),I=1,IV)
        DO 50 I=1,IV
          DO 50 J=1,3
            JAC(KV(I),J)=0
   50   CONTINUE

        DO 60 I=1,NTA
          IELM=KV(I)
          DO 70 J=1,3
            MTJ(IELM,J)=IENA(I,J)
            IF(JEEA(I,J).NE.0)THEN
              JAC(IELM,J)=KV(JEEA(I,J))
            END IF
   70     CONTINUE
   60   CONTINUE

        DO 80 I=1,NTB
          IELM=KV(NTA+I)
          DO 90 J=1,3
            MTJ(IELM,J)=IENB(I,J)
            IF(JEEB(I,J).NE.0)THEN
              JAC(IELM,J)=KV(NTA+JEEB(I,J))
            END IF
   90     CONTINUE
   80   CONTINUE

        DO 100 I=1,IV
          IELM=KV(I)
          DO 110 J=1,3
            IVX=MTJ(IELM,J)
            CALL INCR(IVX,IELM,NER,KTJ,KCM,IERR)
            IF (IERR.NE.0) RETURN
  110     CONTINUE
  100   CONTINUE

*        WRITE(*,*) '[POLY] KV:',(KV(I),I=1,IV)

        DO 120 I=1,IV
          IELM=KV(I)
*          WRITE(*,*) '[POLY] I/IV=',I,'/',IV,' KV(I)=',KV(I)
          DO 130 J=1,3
            JELM=JAC(IELM,J)
*            WRITE(*,*) '[POLY] J=',J,' JELM=',JELM
            IF(JELM.NE.0)GO TO 130
            IPS=MTJ(IELM,J)
            IPG=MTJ(IELM,MOD(J,3)+1)
            NEICNT=0
            DO K=1,NLCNT
              DO L=1,3
                IF(MTJ(NLIST(K),L).EQ.IPG)THEN
                  NEICNT=NEICNT+1
                  NEI(NEICNT)=NLIST(K)
                ENDIF
              ENDDO
            ENDDO
*            WRITE(*,*) '[POLY] IPG=',IPG,' NER(1)=',NER(1)
*            WRITE(*,*) '[POLY] NEI=',(NEI(K),K=1,NEICNT)
            DO 140 K=1,NEICNT
              KELM=NEI(K)
              DO 150 L=1,3
                IVA=MTJ(KELM,L)
                IVB=MTJ(KELM,MOD(L,3)+1)
                IF((IVA.EQ.IPG).AND.(IVB.EQ.IPS))THEN
                  JAC(IELM,J)=KELM
                  JAC(KELM,L)=IELM
                  GO TO 130
                END IF
  150         CONTINUE
  140       CONTINUE
  130     CONTINUE
  120   CONTINUE
*        CALL OUTEX(100,NELM,MTJ,JAC,IDM,PX,PY,KTJ)
*        CALL CHECK(NELM,MTJ,JAC,KTE,IERR)
      END IF

      RETURN
      END
C     ------   SUBROUTINE BREAKLINE   ----------------------
C
C     PURPOSE : 
C
      SUBROUTINE BREAKLINE(IZ,NP,NNB,LIST,NTP,PX,PY,NER,NELM,ICHECK,
     &               NEX,IBEX,IBNO,NBK,IBREAK,NBREAK,MTJ,JAC,IDM,JADRES,
     &           IAN,DPP,ISTACK,IONL,KV,IBR,MAP,NODE,KBD,KTE,KTJ,KCM,PD,
     &           IERR)
      IMPLICIT NONE
C     ==args==
      INTEGER KTJ,KTE,KBD,KCM
      INTEGER IZ,NP,NNB,LIST(KTJ),NTP
      INTEGER NELM,ICHECK(KTE),NEX,IBEX(KBD),IBNO(KBD,KTJ),NBK
      INTEGER IBREAK(KBD),NBREAK(KBD,KTJ),MTJ(KTE,3),JAC(KTE,3),IDM(KTE)
      INTEGER JADRES(KTJ+3,KBD),IAN(KTJ+3),ISTACK(KTJ),IONL(KTJ),KV(KTE)
      INTEGER IBR(KTE),MAP(KTE),NODE,IERR,NER(KTJ+3)
      REAL*8 PX(KTJ+3),PY(KTJ+3),PD(KTJ+3),DPP
C     ==temp==
      INTEGER I,J,K
      INTEGER IFIX(KTJ),IADRES(KTJ+3),INP,IS,LOC,IP,IQ,IN,IFI,IW,IELM
      INTEGER JELM,IED,JED,IV,NEI(KCM),NEICNT
      REAL*8 XS,YS

C     INITIALIZATION
      CALL ARRAYSET(IFIX, KTJ, 0)
      CALL ARRAYSET(IADRES, KTJ+3, 0)
      INP=0
      IS=0
      LOC=0
      IP=0
      IQ=0
      IN=0
      IFI=0
      IW=0
      IELM=0
      JELM=0
      IED=0
      JED=0
      IV=0
      XS=0.D0
      YS=0.D0

      DO 20 I=1,NP
        IADRES(LIST(I))=I
   20 CONTINUE
      CALL ARRAYSET(IFIX, NODE, 1)

C     CREATE THREE TRIANGLES BY INSERTION OF THE FIRST POINT
*      WRITE(*,*) '[BREAKLINE] CREATE THREE TRIANGLES'
      IS=LIST(1)
      XS=PX(IS)
      YS=PY(IS)
*      WRITE(*,*)'[BREAKLINE] IS=',IS
      CALL LOCATE(XS,YS,PX,PY,MTJ,JAC,NELM,KTE,KTJ,LOC,IERR)
      IF(IERR.NE.0)RETURN
      IZ=IDM(LOC)
C     ���̂P�s��IS�����ɍ���Ă���_���ǂ����̃`�F�b�N
      IF(ICHECK(IS).NE.1)THEN
        CALL DISTANCE(NEX,LOC,IS,XS,YS,PX,PY,PD,IBEX,IBNO,
     &              IONL,DPP,MTJ,KBD,KTE,KTJ)
        ICHECK(IS)=1
        INP=INP+1
*        WRITE(*,*)'[BREAKLINE] IONL(IS)=',IONL(IS)
        IF(IONL(IS).EQ.0)THEN
C         �ӏ�łȂ��ꍇ�A
          CALL BDELAUN(IZ,IS,IS,NTP,PX,PY,NER,NELM,MTJ,JAC,
     &                 IDM,JADRES,IAN,ISTACK,KBD,KTE,KTJ,KCM,IERR)
          IF(IERR.NE.0)RETURN
        ELSE
C         �ӏ�̏ꍇ
          CALL WLOCATE(XS,YS,PX,PY,MTJ,JAC,NELM,KTE,KTJ,IELM,JELM,IERR)
          IF(IERR.NE.0)RETURN
          CALL EDGE(IELM,JELM,JAC,KTE,IED,IERR)
          IF (IERR.NE.0) RETURN
          CALL EDGE(JELM,IELM,JAC,KTE,JED,IERR)
          IF (IERR.NE.0) RETURN
          CALL REMESH(IELM,IED,JELM,JED,IS,PX,PY,NER,
     &                 JADRES,IAN,IFIX,NELM,MTJ,JAC,IDM,ISTACK,KTE,KTJ,
     &                 KCM,KBD,NBK,IBREAK,NBREAK,1,IERR)
          IF(IERR.NE.0)RETURN
        END IF
      END IF

C     LOOP OVER EACH POINT AND CONSTRUCT BOUNDARY EDGES
*      WRITE(*,*) '[BREAKLINE]  LOOP OVER EACH POINT'
      I=0
   30 CONTINUE
        I=I+1
        IP=LIST(MOD(I,NP)+1)
        IQ=LIST(I)
        XS=PX(IP)
        YS=PY(IP)
*        WRITE(*,*)'[BREAKLINE] IP=',IP,' IQ=',IQ,' I=',I

*        IF(IP.EQ.19.AND.IQ.EQ.18)THEN
*          CALL OUTEX(100,NELM,MTJ,JAC,IDM,PX,PY,KTJ)
*        ENDIF
*        WRITE(*,*)'[BREAKLINE] CALL CHECK'
*        CALL CHECK(NELM,MTJ,JAC,KTE,IERR)
        IF(IERR.NE.0) RETURN
*        IF(I.EQ.12) THEN
*          WRITE(*,*)'[BREAKLINE] LIST=',(LIST(J),J=1,NP)
*          WRITE(*,*)'[BREAKLINE] XS=',XS,' YS=',YS
*        ENDIF
C       �_P���܂ގO�p�`��T��
*        WRITE(*,*) '[BREAKLINE] CALL LOCATE'
        CALL LOCATE(XS,YS,PX,PY,MTJ,JAC,NELM,KTE,KTJ,LOC,IERR)
        IF(IERR.NE.0)RETURN
*        IF(I.EQ.8) THEN
*          WRITE(*,*) '[BREAKLINE] LOC=',LOC
*        ENDIF
        IF(ICHECK(IP).NE.1)THEN
C       �ǉ��_�������ӂɋ߂��ꍇ�́A�����ӏ�Ɉړ�
*        WRITE(*,*) '[BREAKLINE] CALL DISTANCE'
        CALL DISTANCE(NEX,LOC,IP,XS,YS,PX,PY,PD,IBEX,IBNO,
     &              IONL,DPP,MTJ,KBD,KTE,KTJ)
        ENDIF
        IN=I
        IFI=NEX
        IW=0
*        IF(I.EQ.12) THEN
*          WRITE(*,*) '[BREAKLINE] IN=',I,' IFI=',IFI
*        ENDIF

*        WRITE(*,*) '[BREAKLINE] CALL INTERSECTION1'
        CALL INTERSECTION(NP,IN,IP,IQ,PX,PY,NTP,IFI,NEX,NNB,IBEX,IBNO,
     &                  IBREAK,NBREAK,LIST,IONL,DPP,JADRES,IADRES,
     &                  IAN,KTJ,KBD,IW,PD,KTE,NELM,MTJ,JAC,IDM)
        IFI=NNB-1
        IW=1
*        IF(I.EQ.12) THEN
*          WRITE(*,*) '[BREAKLINE] IFI=',IFI
*        ENDIF
*        WRITE(*,*) '[BREAKLINE] CALL INTERSECTION2'
        CALL INTERSECTION(NP,IN,IP,IQ,PX,PY,NTP,IFI,NBK,NNB,IBREAK,
     &                NBREAK,IBREAK,NBREAK,LIST,IONL,DPP,JADRES,IADRES,
     &                IAN,KTJ,KBD,IW,PD,KTE,NELM,MTJ,JAC,IDM)
        XS=PX(IP)
        YS=PY(IP)
*        WRITE(*,*) '[BREAKLINE] CALL LOCATE'
        CALL LOCATE(XS,YS,PX,PY,MTJ,JAC,NELM,KTE,KTJ,LOC,IERR)
        IF(IERR.NE.0)RETURN
        IF(IONL(IP).NE.1)THEN
          IZ=IDM(LOC)
        END IF

        IF((ICHECK(IP).NE.1).AND.(I.NE.NP))THEN
          ICHECK(IP)=1
          INP=INP+1
          IF(IONL(IP).EQ.0)THEN
*            WRITE(*,*) '[BREAKLINE] CALL BDELAUN'
            CALL BDELAUN(IZ,IP,IP,NTP,PX,PY,NER,NELM,MTJ,
     &               JAC,IDM,JADRES,IAN,ISTACK,KBD,KTE,KTJ,KCM,IERR)
            IF(IERR.NE.0)RETURN
*            WRITE(*,*)'[BREAKLINE] CALL CHECK'
*            CALL CHECK(NELM,MTJ,JAC,KTE,IERR)
          ELSE
*            WRITE(*,*) '[BREAKLINE] CALL WLOCATE'
            CALL WLOCATE(XS,YS,PX,PY,MTJ,JAC,NELM,KTE,KTJ,IELM,JELM
     &                  ,IERR)
*            WRITE(*,*) '[BREAKLINE] IELM=',IELM,', JELM=',JELM
            IF(IERR.NE.0)RETURN
*            WRITE(*,*) '[BREAKLINE] CALL EDGE'
            CALL EDGE(IELM,JELM,JAC,KTE,IED,IERR)
            IF (IERR.NE.0) RETURN
*            WRITE(*,*) '[BREAKLINE] CALL EDGE'
            CALL EDGE(JELM,IELM,JAC,KTE,JED,IERR)
            IF (IERR.NE.0) RETURN
*            WRITE(*,*) '[BREAKLINE] CALL REMESH'
            CALL REMESH(IELM,IED,JELM,JED,IP,PX,PY,NER,
     &                 JADRES,IAN,IFIX,NELM,MTJ,JAC,IDM,ISTACK,KTE,KTJ,
     &                 KCM,KBD,NBK,IBREAK,NBREAK,1,IERR)
            IF(IERR.NE.0)RETURN
*            WRITE(*,*)'[BREAKLINE] CALL CHECK'
*            CALL CHECK(NELM,MTJ,JAC,KTE,IERR)
*            WRITE(*,*) '[BREAKLINE] CALL READRES'
            CALL READRES(JADRES,IAN,NBK,IBREAK,NBREAK,KBD,KTJ)
          END IF
        END IF

*        IF(IP.EQ.19.AND.IQ.EQ.18)THEN
*          CALL OUTEX(100,NELM,MTJ,JAC,IDM,PX,PY,KTJ)
*        ENDIF

*        WRITE(*,*) '[BREAKLINE] 40 LOOP - ',JNB(IQ)

        CALL NEIGH(IERR,KTJ,KTE,KCM,MTJ,JAC,NER,IQ,NEI,NEICNT)
        DO 40 J=1,NEICNT
          DO 40 K=1,3
            IF(MTJ(NEI(J),K).EQ.IP) GO TO 50
   40   CONTINUE


C       SEARCH FOR TRIANGLES BETWEEN NEW DATA POINT AND
C       OLD DATA POINT
*        WRITE(*,*) '[BREAKLINE] CALL SEARCH'
        CALL SEARCH(IZ,IP,IQ,NER,NELM,MTJ,JAC,IDM,PX,PY,
     &             ISTACK,IV,KV,IADRES,IBR,KTE,KTJ,KCM,IERR)
        IF (IERR.NE.0) RETURN
*        WRITE(*,*) '[BREAKLINE] IV=',IV
*        DO J=1,IV
*          WRITE(*,*) '[BREAKLINE] ',KV(J),':',(MTJ(KV(J),K),K=1,3)
*        ENDDO

C       REMOVE TRIANGLES BETWEEN NEW DATA POINT AND 
C       OLD DATA POINTAND SUBDIVIDE
C       THE POLYGON CONSISTED BY THEM
*        WRITE(*,*) '[BREAKLINE] CALL POLY'
        CALL POLY(IQ,IP,IV,KV,PX,PY,NELM,MTJ,JAC,IDM,NER,
     &            MAP,KTE,KTJ,KCM,IERR)
        IF (IERR.NE.0) RETURN

   50 CONTINUE
      IF(I.LT.IBREAK(NNB)-1) GO TO 30

*      WRITE(*,610)' NODE NUMBER = ',LIST(NP)
      NODE=NODE+INP

      RETURN
      END
C     ------   SUBROUTINE READRES   ------------------------
C
C     PURPOSE : 
C     ARGUMENTS:
C       JADRES - �ߓ_-���E�_�����}�b�s���O�z��
C       IAN    - �ߓ_-���E�_�Ƃ��ė��p����Ă����
C       NBK    - �J���E��
C       IBREAK - �e�J���E���\������_��
C       NBREAK - �e�J���E���\������_�ԍ����X�g
C
      SUBROUTINE READRES(JADRES,IAN,NBK,IBREAK,NBREAK,KBD,KTJ)
      IMPLICIT NONE
C     ==args==
      INTEGER KTJ,KBD
      INTEGER JADRES(KTJ+3,KBD),IAN(KTJ+3),NBK,IBREAK(KBD)
      INTEGER NBREAK(KBD,KTJ)
C     ==temp==
      INTEGER I,J
      INTEGER NP,JJ,IBK,NB,JB,LIST(KTJ)

      CALL ARRAYSET(IAN,KTJ+3,0)
      CALL ARRAYSET(JADRES,(KTJ+3)*KBD,0)

      NP=0
      JJ=0
      IBK=0
      NB=0
      JB=0
      CALL ARRAYSET(LIST,KTJ,0)

      DO 20 I=1,NBK
        IBK=NP+1
        NB=I
        NP=IBREAK(NB)+IBK
        DO 30 J=IBK,NP-1
          JB=J+1-IBK
          JJ=JJ+1
          LIST(JJ)=NBREAK(NB,JB)
   30   CONTINUE
        JJ=JJ+1
        LIST(JJ)=1
   20 CONTINUE

      DO 40 I=1,NP
        IAN(LIST(I))=IAN(LIST(I))+1
        JADRES(LIST(I),IAN(LIST(I)))=I+1
   40 CONTINUE

      DO 50 I=1,KBD
        JADRES(1,I)=0
   50 CONTINUE

      RETURN
      END
C     ------   SUBROUTINE DISTANCE   -----------------
C
C     PURPOSE : 
C       �����ӂɋ߂��_�����̕ӏ�Ɉړ�������
C     OUT:
C       XP,YP,IONL,PX,PY,PD,IBEX,IBNO
      SUBROUTINE DISTANCE(NEX,LOC,IS,XP,YP,PX,PY,PD,
     &                   IBEX,IBNO,IONL,DPP,MTJ,KBD,KTE,KTJ)
C      SUBROUTINE DISTANCE(NEX,NODE,INP,LOC,IS,NNB,XP,YP,PX,PY,PD,
C     &                   IBEX,IBNO,IBREAK,NBREAK,LIST,IONL,DPP,NTP,
C     &                   MTJ,KBD,KTE,KTJ)
      IMPLICIT NONE
C     ==args==
      INTEGER KTE,KTJ,KBD
      INTEGER NEX,LOC,IS,IBEX(KBD),IBNO(KBD,KTJ),IONL(KTJ),MTJ(KTE,3)
      REAL*8 XP,YP,PX(KTJ+3),PY(KTJ+3),PD(KTJ+3),DPP
C     ==temp==
      INTEGER I,J,K
      INTEGER IA,IB,IC,ID,IE,IBC,NB,NPN
      REAL*8 DNEA,PXP,PYP,DIS,XX,YY

      IA=MTJ(LOC,1)
      IB=MTJ(LOC,2)
      IC=MTJ(LOC,3)
      ID=0
      IE=0
      IBC=0
      NB=0
      NPN=0
      DNEA=1.0D0
      PXP=0.0D0
      PYP=0.0D0
      DIS=0.0D0
      XX=0.0D0
      YY=0.0D0

C     �v�fLOC��3�ӂ̂����A�_P�Ɉ�ԋ߂��ӂ�T���A�������v�Z����
C     ��AB�Ɠ_P�̋������v�Z
      CALL MEA(DIS,IA,IB,XP,YP,XX,YY,PX,PY,KTJ)
      IF(DIS.LT.DNEA)THEN
        DNEA=DIS
        ID=IA
        IE=IB
        PXP=XX
        PYP=YY
      END IF
C     ��BC�Ɠ_P�̋������v�Z
      CALL MEA(DIS,IB,IC,XP,YP,XX,YY,PX,PY,KTJ)
      IF(DIS.LT.DNEA)THEN
        DNEA=DIS
        ID=IB
        IE=IC
        PXP=XX
        PYP=YY
      END IF
C     ��CA�Ɠ_P�̋������v�Z
      CALL MEA(DIS,IC,IA,XP,YP,XX,YY,PX,PY,KTJ)
      IF(DIS.LT.DNEA)THEN
        DNEA=DIS
        ID=IC
        IE=IA
        PXP=XX
        PYP=YY
      END IF

      IBC = 0
      IF(DSQRT(DNEA).LT.DPP)THEN
C       �w�肳�ꂽ�ߓ_�Ԋu��菬�����ꍇ�A
        XP=PXP
        YP=PYP
        PX(IS)=PXP
        PY(IS)=PYP
        IONL(IS)=1  !�ӏ�ɓ_������
        IBC=1
C       �ӏ�ɓ_��ǉ��B�_�̈ʒu�͓_���x����v�Z
        IF(ID.EQ.IA)THEN
          PD(IS)=(PD(IA)+PD(IB))*0.5D0
        ELSE IF(ID.EQ.IB)THEN
          PD(IS)=(PD(IB)+PD(IC))*0.5D0
        ELSE IF(ID.EQ.IC)THEN
          PD(IS)=(PD(IC)+PD(IA))*0.5D0
        ELSE
C          WRITE(*,'('' ***ERROR IN SUBROUTINE NEAR***'')')
C         WRITE(*,'('' ***INCORRECT NUMBER OF
C     &                             ELEMENTS FORMED***'')')
        END IF
      ELSE
        IONL(IS)=0
      END IF

C     �ӏ�̏ꍇ�A���E�����X�V
      IF(IBC.EQ.1)THEN
        DO 30 I=1,NEX
          NB=IBEX(I)
          DO 40 J=1,NB
            IF((IBNO(I,J).EQ.ID.AND.IBNO(I,MOD(J,NB)+1).EQ.IE).OR.
     &        (IBNO(I,J).EQ.IE.AND.IBNO(I,MOD(J,NB)+1).EQ.ID))THEN
              IBEX(I)=IBEX(I)+1
              NPN=IBEX(I)
              DO 50 K=NPN,J+2,-1
               IBNO(I,K)=IBNO(I,K-1)
   50         CONTINUE
              IBNO(I,J+1)=IS
            END IF
   40     CONTINUE
   30   CONTINUE
      END IF

      RETURN
      END
C     ------   SUBROUTINE INTERSECTION   -----------------
C 
C     PURPOSE : 
C       �J���E�Ƒ��̋��E���������Ă��邩�`�F�b�N���A�������Ă���ꍇ�A
C       ���E�����X�V����B
C
C     ARGUMENTS:
C       IN  - 
C       NTP - 
C       IFI - ���E��?
C       IW  - 
      SUBROUTINE INTERSECTION(NP,IN,IP,IQ,PX,PY,NTP,IFI,NEX,NNB,IBEX,
     &                        IBNO,IBREAK,NBREAK,LIST,IONL,DPP,JADRES,
     &                        IADRES,IAN,KTJ,KBD,IW,PD
     &                        ,KTE,NELM,MTJ,JAC,IDM)
      IMPLICIT NONE
C     == arguments ==
      INTEGER KTJ,KBD
      INTEGER NP,IN,IP,IQ,NTP,IFI,NEX,NNB,IBEX(KBD),IBNO(KBD,KTJ)
      INTEGER IBREAK(KBD),NBREAK(KBD,KTJ),LIST(KTJ),IONL(KTJ)
      INTEGER JADRES(KTJ+3,KBD),IADRES(KTJ+3),IAN(KTJ+3),IW
      REAL*8 PX(KTJ+3),PY(KTJ+3),DPP,PD(KTJ+3)
      INTEGER KTE,NELM,MTJ(KTE,3),JAC(KTE,3),IDM(KTE)
C     == temp variables ==
      INTEGER I,J
      INTEGER IWC,IPED,IDEN,IHORA,IHORC,NB,NNP,IR,IS
      REAL*8 XMAX,XMIN,YMAX,YMIN,A,B,C,D,XX,YY

      INTEGER IDEBUG
      IF(IP.EQ.42.AND.IQ.EQ.50) THEN
        IDEBUG=1
      ELSE
        IDEBUG=0
      ENDIF

      IWC=0
      IPED=0
      IDEN=0
      IHORA=0
      IHORC=0
      NB=0
      NNP=0
      IR=0
      IS=0
      XMAX=0.D0
      XMIN=0.D0
      YMAX=0.D0
      YMIN=0.D0
      A=0.D0
      B=0.D0
      C=0.D0
      D=0.D0
      XX=0.D0
      YY=0.D0

  300 CONTINUE

*      WRITE(*,*) '[INTERSECTION] 300 CONTINUE'

      XMAX=DMAX1(PX(IP),PX(IQ))
      XMIN=DMIN1(PX(IP),PX(IQ))
      YMAX=DMAX1(PY(IP),PY(IQ))
      YMIN=DMIN1(PY(IP),PY(IQ))

C     ��IP-IQ�̊֐����擾
      CALL EQUATION(IP,IQ,PX,PY,A,B,IDEN,IHORA,KTJ)
      IF(IDEN.EQ.1)THEN
C       Y���ɕ��s�ȏꍇ
        GO TO 100
      END IF

C     ��IP-IQ��Y���ɕ��s�łȂ��ꍇ�̌����v�Z

C     IWC BreakLine�͕��Ȃ��̂ŋ��E����J�̃��[�v���P���炳�Ȃ�
C     �Ƃ����Ȃ�
C     ���ł�NNP���P���炳�Ȃ��Ƃ����Ȃ����A���̈ʒu���ƕ����񌸂炵��
C     ���܂��̂łP�񌸂炵�����ǂ����̃`�F�b�N�悤�ɗ^����B
      IWC=0

*      WRITE(*,*) '[INTERSECTION] LOOP1'
      DO 10 I=1,IFI
        NB=I
        NNP=IBEX(NB)
        IF(IW.EQ.1.AND.IWC.EQ.0)THEN
          NNP=NNP-1
          IWC=1
        END IF
        DO 20 J=1,NNP
          IR=IBNO(NB,MOD(J,IBEX(NB))+1)
          IS=IBNO(NB,J)
          CALL EQUATION(IR,IS,PX,PY,C,D,IDEN,IHORC,KTJ)
          IF(IDEN.EQ.1) GO TO 200

          XX=(D-B)/(A-C)
          IF(((PX(IR).LE.XX.AND.PX(IS).GE.XX).OR.
     &        (PX(IS).LE.XX.AND.PX(IR).GE.XX)).AND.
     &              (XX.LE.XMAX.AND.XX.GE.XMIN))THEN
            YY=A*XX+B
            CALL INCREASES(NP,IN,IP,I,J,XX,YY,PX,PY,NTP,NNB,IPED,IONL,
     &                       DPP,IBEX,IBNO,NEX,IBREAK,NBREAK,
     &                       LIST,JADRES,IADRES,IAN,KTJ,KBD,IW,PD)
            IF(IPED.EQ.1)THEN
              IPED=0
              GO TO 20
            END IF
            GO TO 300
          END IF
          GO TO 20

  200     CONTINUE
          IF(C.LE.XMAX.AND.C.GE.XMIN)THEN
            YY=A*C+B
            IF(IHORA.NE.1.AND.YY.LT.YMAX.AND.YY.GT.YMIN.AND.
     &                     (PY(IR).LE.YY.AND.PY(IS).GE.YY).OR.
     &                     (PY(IS).LE.YY.AND.PY(IR).GE.YY))THEN
              CALL INCREASES(NP,IN,IP,I,J,C,YY,PX,PY,NTP,NNB,IPED,IONL,
     &                       DPP,IBEX,IBNO,NEX,IBREAK,NBREAK,
     &                       LIST,JADRES,IADRES,IAN,KTJ,KBD,IW,PD)
              IF(IPED.EQ.1)THEN
                IPED=0
                GO TO 20
              END IF
              GO TO 300
            ELSE IF(IHORA.EQ.1)THEN
              IF((PY(IR).LE.YY.AND.PY(IS).GE.YY).OR.
     &               (PY(IS).LE.YY.AND.PY(IR).GE.YY))THEN
                CALL INCREASES(NP,IN,IP,I,J,C,YY,PX,PY,NTP,NNB,IPED,
     &                       IONL,DPP,IBEX,IBNO,NEX,IBREAK,NBREAK,
     &                       LIST,JADRES,IADRES,IAN,KTJ,KBD,IW,PD)
                IF(IPED.EQ.1)THEN
                  IPED=0
                  GO TO 20
                END IF
                GO TO 300
              END IF
            END IF
          END IF
   20   CONTINUE
   10 CONTINUE
      RETURN


  100 CONTINUE
C     ��IP-IQ��Y���ɕ��s�ȏꍇ�̌����v�Z
*      IF(IDEBUG.EQ.1) THEN
*        WRITE(*,*)'[INTERSECTION] LIST=',(LIST(J),J=1,NP)
*        CALL OUTEX(NEX,100,NELM,MTJ,JAC,IDM,PX,PY,KTJ)
*      ENDIF
*      WRITE(*,*) '[INTERSECTION] LOOP2: IFI=',IFI
      DO 110 I=1,IFI
        NB=I  ! Closed Boundary No.
        NNP=IBEX(NB)  ! Numer of vertex on this boundary
        IF(IW.EQ.1.AND.IWC.EQ.0)THEN
          NNP=NNP-1
          IWC=1
        END IF
*        WRITE(*,*) '[INTERSECTION] NB=',NB,' NNP=',NNP,' IWC=',IWC
        DO 120 J=1,NNP
          IR=IBNO(NB,MOD(J,IBEX(NB))+1)
          IS=IBNO(NB,J)
*          WRITE(*,*) '[INTERSECTION] CALL EQUATION: IR=',IR,' IS=',IS
          CALL EQUATION(IR,IS,PX,PY,C,D,IDEN,IHORC,KTJ)
C         ��IR-IS��Y���ɕ��s�łȂ��ꍇ�̂݁A�����v�Z
          IF(IDEN.NE.1)THEN 
            YY=C*A+D
            IF(IHORC.NE.1) THEN
C               ��IR-IS��X���ɕ��s�łȂ��ꍇ
              IF (YY.LT.YMAX.AND.YY.GT.YMIN.AND.
     &                   ((PX(IR).LE.A.AND.PX(IS).GE.A).OR.
     &                    (PX(IS).LE.A.AND.PX(IR).GE.A)))THEN
*                WRITE(*,*) '[INTERSECTION] CALL INCREASES1'
                CALL INCREASES(NP,IN,IP,I,J,A,YY,PX,PY,NTP,NNB,IPED,
     &                       IONL,DPP,IBEX,IBNO,NEX,IBREAK,NBREAK,
     &                       LIST,JADRES,IADRES,IAN,KTJ,KBD,IW,PD)
                IF(IPED.EQ.1)THEN
                  IPED=0
                  GO TO 120
                END IF
                GO TO 300
              ENDIF
            ELSE IF(IHORC.EQ.1)THEN
C              ��IR-IS��X���ɕ��s�̏ꍇ
              IF((YY.LE.YMAX.AND.YY.GE.YMIN).AND.
     &         ((PX(IR).LE.A.AND.PX(IS).GE.A).OR.
     &          (PX(IS).LE.A.AND.PX(IR).GE.A)))THEN
                IF(IDEBUG.EQ.1) THEN
*                  WRITE(*,*) 'YY=',YY,' YMAX=',YMAX,' YMIN=',YMIN
*                  WRITE(*,*) 'PX(IR)=',PX(IR),' PX(IS)=',PX(IS),' A=',A
                ENDIF

*                WRITE(*,*) '[INTERSECTION] CALL INCREASES2'
                CALL INCREASES(NP,IN,IP,I,J,A,YY,PX,PY,NTP,NNB,IPED,
     &                       IONL,DPP,IBEX,IBNO,NEX,IBREAK,NBREAK,
     &                       LIST,JADRES,IADRES,IAN,KTJ,KBD,IW,PD)
                IF(IPED.EQ.1)THEN
                  IPED=0
                  GO TO 120
                END IF
                GO TO 300
              END IF
            END IF
          END IF
  120   CONTINUE
  110 CONTINUE

      RETURN

      END
C     ------   SUBROUTINE INCREASES   --------------------
C
C     PURPOSE : 
C       �w�肳�ꂽ�J���E�����̋��E�ƌ�������ꍇ�ɁA�J���E�𕪂���B
C     ARGUMENTS:
C       NP
C       IN
C       IP  - Vertex P
C       IEX - Closed Boundary No.
C       INO - Vertex No. of Closed Boundary
C       XX  -
C       YY
C       NNB - Numer of vertex on this boundary
C       IPED - 
C
      SUBROUTINE INCREASES(NP,IN,IP,IEX,INO,XX,YY,PX,PY,NTP,NNB,IPED,
     &                     IONL,DPP,IBEX,IBNO,NEX,IBREAK,NBREAK,LIST,
     &                     JADRES,IADRES,IAN,KTJ,KBD,IW,PD)
*      IMPLICIT REAL*8(A-H,O-Z)
      IMPLICIT NONE
C     == args ==
      INTEGER KTJ,KBD
      INTEGER NP,IN,IP,IEX,INO,NTP,NNB,IPED,IONL(KTJ),IBEX(KBD)
      INTEGER IBNO(KBD,KTJ),NEX,IBREAK(KBD),NBREAK(KBD,KTJ),LIST(KTJ)
      REAL*8 XX,YY,PX(KTJ+3),PY(KTJ+3),DPP,PD(KTJ+3)
      INTEGER JADRES(KTJ+3,KBD),IADRES(KTJ+3),IAN(KTJ+3),IW
C     ==========
      INTEGER I,J,K,NPN,IK,IS,IR1,IB,IJ,IR
      REAL*8 CX,CY,DIF

      NPN=0
      IK=0
      IS=0
      IR1=0
      IB=0
      IJ=0
      IR=0
      CX=0.0D0
      CY=0.0D0
      DIF=0.0D0

*      WRITE(*,*) '[INCREASE] IP=',IP,' IW=',IW,' XX=',XX,' YY=',YY
      IF(IW.EQ.1)GO TO 200

C     ���e�ŏ��_�Ԋu��菬�����_���`�F�b�N
      DO 100 I=1,NTP
        CX=PX(I)
        CY=PY(I)
        DIF=CX*CX+XX*XX+CY*CY+YY*YY-2.0D0*(CX*XX+CY*YY)
        DIF=DSQRT(DABS(DIF))
*        WRITE(*,*) '[INCREASE] CX=',CX,' CY=',CY,' DIF=',DIF,'DPP=',DPP
        
        IF(DIF.LT.DPP)THEN
*          WRITE(*,*) '[INCREASE] I=',I,' DIF=',DIF
          IF(NBREAK(NNB,IN).NE.I.AND.NBREAK(NNB,IN+1).NE.I)THEN
C           �J���E�����č\�z
            IBREAK(NNB)=IBREAK(NNB)+1
            NPN=IBREAK(NNB)
            DO 110 K=NPN,IN+2,-1
              NBREAK(NNB,K)=NBREAK(NNB,K-1)
              LIST(K)=LIST(K-1)
  110       CONTINUE
            NBREAK(NNB,IN+1)=I
            LIST(IN+1)=I
            DO 120 J=1,KTJ+3
              IADRES(J)=0
  120       CONTINUE
            DO 130 J=1,NP+1
              IADRES(LIST(J))=J
  130       CONTINUE
            NP=NP+1
            IP=LIST(MOD(IN,NP)+1)
          ELSE
            IPED=1
          END IF
*          WRITE(*,*) '[INCREASES] SAME POINT'
          RETURN
        END IF
  100 CONTINUE

C     #�J���E
      NTP=NTP+1
      IBREAK(NNB)=IBREAK(NNB)+1
      NPN=IBREAK(NNB)
      DO 10 K=NPN,IN+2,-1
        NBREAK(NNB,K)=NBREAK(NNB,K-1)
        LIST(K)=LIST(K-1)
   10 CONTINUE
      NBREAK(NNB,IN+1)=NTP
      LIST(IN+1)=NTP
      PX(NTP)=XX
      PY(NTP)=YY
      IONL(NTP)=1
      DO 20 I=1,KTJ+3
        IADRES(I)=0
   20 CONTINUE
      DO 30 I=1,NP+1
        IADRES(LIST(I))=I
   30 CONTINUE
      NP=NP+1
      IP=LIST(MOD(IN,NP)+1)
C     #���E
      IBEX(IEX)=IBEX(IEX)+1
      IK=IBEX(IEX)
      DO 50 I=IK,INO+2,-1
        IBNO(IEX,I)=IBNO(IEX,I-1)
   50 CONTINUE
      IBNO(IEX,INO+1)=NTP
C     SEARCHES FOR THE SHARED LINE
      IS=IBNO(IEX,INO)
      IR1=MOD(INO+2,IBEX(IEX))
      IF(IR1.EQ.0) IR1=IBEX(IEX)
      IR=IBNO(IEX,IR1)
      PD(NTP)=(PD(IS)+PD(IR))*0.5D0
      DO 60 I=1,NEX
        IF(I.EQ.IEX) GO TO 60
        IB=IBEX(I)
        DO 70 J=1,IB
          IJ=MOD(J,IB)+1
          IF(IBNO(I,IJ).EQ.IS.AND.IBNO(I,J).EQ.IR)THEN
            IBEX(I)=IBEX(I)+1
            IK=IBEX(I)
            DO 80 K=IK,IJ+1,-1
              IBNO(I,K)=IBNO(I,K-1)
   80       CONTINUE
            IBNO(I,IJ)=NTP
          END IF
   70   CONTINUE
   60 CONTINUE
*      WRITE(*,*) '[INCREASES] RETURN CLOSED BOUNDARY'
      RETURN


  200 CONTINUE
C     #�J���E
      DO 210 I=1,NTP
        CX=PX(I)
        CY=PY(I)
        DIF=CX*CX+XX*XX+CY*CY+YY*YY-2*(CX*XX+CY*YY)
        DIF=DSQRT(DABS(DIF))
        IF(DIF.LT.DPP)THEN
          IF(NBREAK(NNB,IN).NE.I.AND.NBREAK(NNB,IN+1).NE.I)THEN
            IBREAK(NNB)=IBREAK(NNB)+1
            NPN=IBREAK(NNB)
            DO 220 K=NPN,IN+2,-1
              NBREAK(NNB,K)=NBREAK(NNB,K-1)
              LIST(K)=LIST(K-1)
  220       CONTINUE
            NBREAK(NNB,IN+1)=I
            LIST(IN+1)=I
            NP=NP+1
            IP=LIST(MOD(IN,NP)+1)
          ELSE
            IPED=1
          END IF
*          WRITE(*,*) '[INCREASES] RETURN OPEN BOUNDARY1'
          RETURN
        END IF
  210 CONTINUE
      NTP=NTP+1
      IBREAK(NNB)=IBREAK(NNB)+1
      NPN=IBREAK(NNB)
      DO 230 K=NPN,IN+2,-1
        NBREAK(NNB,K)=NBREAK(NNB,K-1)
        LIST(K)=LIST(K-1)
  230 CONTINUE
      NBREAK(NNB,IN+1)=NTP
      LIST(IN+1)=NTP
      PX(NTP)=XX
      PY(NTP)=YY
      IONL(NTP)=2
      NP=NP+1
      IP=LIST(MOD(IN,NP)+1)
C     IBEX IS CHANGED
      IBEX(IEX)=IBEX(IEX)+1
      IK=IBEX(IEX)
      DO 250 I=IK,INO+2,-1
      IBNO(IEX,I)=IBNO(IEX,I-1)
  250 CONTINUE
      IBNO(IEX,INO+1)=NTP
C     SEARCHES FOR THE SHARED LINE
      IS=IBNO(IEX,INO)
      IR=IBNO(IEX,INO+2)
c      PD(NTP)=(PD(IS)+PD(IR))/2
      DO 260 I=1,NEX
        IF(I.EQ.IEX) GO TO 260
        IB=IBEX(I)
        DO 270 J=1,IB
          IJ=MOD(J,IB)+1
          IF(IBNO(I,IJ).EQ.IS.AND.IBNO(I,J).EQ.IR)THEN
            IBEX(I)=IBEX(I)+1
            IK=IBEX(I)
            DO 280 K=IK,IJ+1,-1
              IBNO(I,K)=IBNO(I,K-1)
  280       CONTINUE
            IBNO(I,IJ)=NTP
          END IF
  270   CONTINUE
  260 CONTINUE

      RETURN
      END
C     ------   SUBROUTINE PICK   ---------------------------
C
C     PURPOSE : PICK UP POINTS ON BOUNDARY OF GIVEN POLYGON
C
      SUBROUTINE PICK(IQ,IP,IV,KV,MTJ,JAC,MAP,NPA,NPB,NSRA,
     &               NSRB,KTE,KTJ,LTJ,IERR)
      IMPLICIT NONE
      INTEGER KTE,KTJ,LTJ
      INTEGER IQ,IP,IV,KV(KTE),MTJ(KTE,3),JAC(KTE,3),MAP(KTE)
      INTEGER NPA,NPB,NSRA(LTJ),NSRB(LTJ),IERR
C     == functions
      INTEGER IVERT
C     == temp
      INTEGER I
      INTEGER IVX,JVX,JELM

*      WRITE(*,*) '[PICK] P=',IP,' Q=',IQ
C
C     DOMAIN A
C
*      WRITE(*,*) '[PICK] DOMAIN A'
      NPA=1
      NSRA(NPA)=IP
*      WRITE(*,*) '[PICK] CALL IVERT 1'
      IVX=IVERT(KV(1),IP,MTJ,KTE,IERR)
      IF(IERR.NE.0)RETURN
*      WRITE(*,*) '[PICK] IVX=',IVX,MTJ(KV(1),MOD(IVX,3)+1)
      NPA=NPA+1
      NSRA(NPA)=MTJ(KV(1),MOD(IVX,3)+1)
      DO 10 I=2,IV-1
*        WRITE(*,*) '[PICK] CALL IVERT 2'
        JVX=IVERT(KV(I),NSRA(NPA),MTJ,KTE,IERR)
        IF(IERR.NE.0)RETURN
        JELM=JAC(KV(I),JVX)
*        WRITE(*,*) '[PICK] JELM=',JELM,' MAP(JELM)=',MAP(JELM)
*        WRITE(*,*) '[PICK] MTJ(JELM)=',(MTJ(JELM,J),J=1,3)
        IF(JELM.EQ.0)THEN
C         �אڗv�f���Ȃ��ꍇ
          NPA=NPA+1
          NSRA(NPA)=MTJ(KV(I),MOD(JVX,3)+1)
        ELSE IF(MAP(JELM).EQ.0)THEN
C         �אڗv�f���̈�A�Ɋ܂܂�Ȃ��ꍇ
          NPA=NPA+1
          NSRA(NPA)=MTJ(KV(I),MOD(JVX,3)+1)
        END IF
   10 CONTINUE
      NPA=NPA+1
      NSRA(NPA)=IQ
C
C     DOMAIN B
C
*      WRITE(*,*) '[PICK] DOMAIN B'
      NPB=1
      NSRB(NPB)=IQ
*      WRITE(*,*) '[PICK] CALL IVERT 3'
      IVX=IVERT(KV(IV),IQ,MTJ,KTE,IERR)
      IF(IERR.NE.0)RETURN
      NPB=NPB+1
      NSRB(NPB)=MTJ(KV(IV),MOD(IVX,3)+1)
      DO 20 I=IV-1,2,-1
*        WRITE(*,*) '[PICK] CALL IVERT 4'
        JVX=IVERT(KV(I),NSRB(NPB),MTJ,KTE,IERR)
        IF(IERR.NE.0)RETURN
        JELM=JAC(KV(I),JVX)
        IF(JELM.EQ.0)THEN
          NPB=NPB+1
          NSRB(NPB)=MTJ(KV(I),MOD(JVX,3)+1)
        ELSE IF(MAP(JELM).EQ.0)THEN
          NPB=NPB+1
          NSRB(NPB)=MTJ(KV(I),MOD(JVX,3)+1)
        END IF
   20 CONTINUE
      NPB=NPB+1
      NSRB(NPB)=IP

C     CHECK OF NODES WHICH IS PICKED UP
      IF(IV.NE.NPA+NPB-4)THEN
        IERR = 206
        RETURN
C        WRITE(*,'('' ***ERROR IN SUBROUTINE PICK***'')')
C        WRITE(*,'('' ***INCORRECT NUMBER OF NODES***'')')
C        PAUSE
C        STOP
      END IF

*      WRITE(*,*) '[PICK] NPA=',NPA,' NPB=',NPB
*      WRITE(*,*) '[PICK] NSRA=',(NSRA(I),I=1,NPA)
*      WRITE(*,*) '[PICK] NSRB=',(NSRB(I),I=1,NPB)

      RETURN
      END
C     ------   SUBROUTINE SUBDIV   -------------------------
C
C     PURPOSE : SUBDIVIDE GIVEN POLYGON USING
C               BY THE MODIFIED-DELAUNAY
C               ���p�`�̈���O�p�`�ɕ�������B�������@�𗘗p�B
C
C     'NTE' -  NUMBER OF ELEMENTS IN DOMAIN
C     'IEN' -  RELATION BETWEEN ELEMENT AND NODES IN DOMAIN
C     'JEE' -  RELATION BETWEEN ELEMENT AND
C              ADJACENT ELEMENTS IN DOMAIN
C
C     NPL   -   �^����ꂽ�̈�̋��E��ߓ_��
C     NSR   -   �^����ꂽ�̈�̋��E��ߓ_�ԍ�
C     NTE   -   �ĕ�����̗̈���v�f��
C     IEN   -   �ĕ�����̗v�f-�ߓ_�֌W
C     JEE   -   �ĕ�����̗v�f-�אڗv�f�֌W
C
C     IX    -   �ĕ�����̕Ӑ�
C     IHEN  -   �ĕ�����̕�-�ߓ_�֌W
C     JHEN  -   �ĕ�����̕�-�v�f�֌W
C
      SUBROUTINE SUBDIV(NPL,NSR,PX,PY,NTE,IEN,JEE,IHEN,JHEN,
     &                 IAD,JSTACK,KTJ,LTE,LTJ,LHN,IERR)
      IMPLICIT NONE
C     ==args==
      INTEGER KTJ,LTE,LTJ,LHN
      INTEGER NPL,NSR(LTJ),NTE,IEN(LTE,3),JEE(LTE,3),IHEN(LHN,2)
      INTEGER JHEN(LHN),IAD(LHN),JSTACK(LTJ),IERR
      REAL*8 PX(KTJ+3),PY(KTJ+3)
C     ==temp==
      INTEGER I,J,K
      INTEGER NNPL,NBS,IA,IB,IC,IX,NNPLPRE
      REAL*8 XA,YA,XB,YB,SEE

C     INITIALIZATION
      NTE=0
      CALL ARRAYSET(IEN,LTE*3,0)
      CALL ARRAYSET(JEE,LTE*3,0)

      NNPL=NPL
      NBS=0
      IA=0
      IB=0
      IC=0
      IX=0
      XA=0.0D0
      YA=0.0D0
      XB=0.0D0
      YB=0.0D0
      SEE=0.0D0

C     COMPUTATION OF IEN
*      WRITE(*,*) '[SUBDIV] COMPUTATION OF IEN'
   20 CONTINUE
*      WRITE(*,*) '[SUBDIV] NNPL=',NNPL
      IF(NNPL.GE.3)THEN
        NBS=1
        NNPLPRE=NNPL
   30   CONTINUE
*        WRITE(*,*) '[SUBDIV] NBS=',NBS
        IF(NBS.LE.NNPL-1)THEN
          IA=NSR(NBS)
          IB=NSR(NBS+1)
          IC=NSR(MOD(NBS+1,NNPL)+1)
*          WRITE(*,*) '[SUBDIV] IA,IB,IC=',IA,IB,IC
          XA=PX(IB)-PX(IA)
          YA=PY(IB)-PY(IA)
          XB=PX(IC)-PX(IA)
          YB=PY(IC)-PY(IA)
          SEE=XA*YB-XB*YA
*          WRITE(*,*) '[SUBDIV] PA=  ',PX(IA),'  ',PY(IA)
*          WRITE(*,*) '[SUBDIV] PB=  ',PX(IB),'  ',PY(IB)
*          WRITE(*,*) '[SUBDIV] PC=  ',PX(IC),'  ',PY(IC)
*          WRITE(*,*) '[SUBDIV] SEE=',SEE
          IF(SEE.GT.1.0D-15)THEN
            NTE=NTE+1
            IEN(NTE,1)=IA
            IEN(NTE,2)=IB
            IEN(NTE,3)=IC
            NNPL=NNPL-1
            DO 40 I=NBS+1,NNPL
              NSR(I)=NSR(I+1)
   40       CONTINUE
          END IF
          NBS=NBS+1
          GO TO 30
        ELSEIF(NBS.EQ.NNPL)THEN
          IF(NNPL.EQ.NNPLPRE)THEN
            IERR=2131
            RETURN
          ENDIF
        END IF
        GO TO 20
      END IF
C
C     COMPUTATION OF JEE
C
*      WRITE(*,*) '[SUBDIV] COMPUTATION OF JEE'
      IX=0
      DO 50 I=1,LHN
        IHEN(I,1)=0
        IHEN(I,2)=0
        JHEN(I)=0
        IAD(I)=0
   50 CONTINUE

      DO 60 I=1,NTE
        DO 70 J=1,3
        IA=IEN(I,J)
        IB=IEN(I,MOD(J,3)+1)
        DO 80 K=1,IX
          IF((IHEN(K,1).EQ.IB).AND.(IHEN(K,2).EQ.IA))THEN
            JEE(I,J)=JHEN(K)
            JEE(JHEN(K),IAD(K))=I
            GO TO 70
          END IF
   80     CONTINUE
        IX=IX+1
        JHEN(IX)=I
        IAD(IX)=J
        IHEN(IX,1)=IA
        IHEN(IX,2)=IB
   70   CONTINUE
   60 CONTINUE
C
C     APPLY THE MODIFIED-DELAUNAY TO NEW ELEMENTS
C
*      WRITE(*,*) '[SUBDIV] CALL LAWSON'
      CALL LAWSON(NTE,IEN,JEE,NPL,PX,PY,JSTACK,KTJ,LTE,LTJ,IERR)
      IF (IERR.NE.0) RETURN

      RETURN
      END
C     ------   SUBROUTINE LAWSON   -------------------------
C
C     PURPOSE : APPLY LAWSON'S SWAPPING ALGORITHM
C               TO NEW ELEMENTS
C     LAST MODIFIED :  4 JUN 1990
C
      SUBROUTINE LAWSON(NTE,IEN,JEE,NPL,PX,PY,JSTACK,KTJ,
     &                 LTE,LTJ,IERR)
      IMPLICIT NONE
C     ==args==
      INTEGER KTJ,LTE,LTJ
      INTEGER NTE,IEN(LTE,3),JEE(LTE,3),NPL,JSTACK(LTJ),IERR
      REAL*8 PX(KTJ+3),PY(KTJ+3)
C     ==functions==
      INTEGER IPUSH
C     ==temp==
      INTEGER I,J
      INTEGER ITOP,MAXSTK,NCOUNT,IELM,IL,JL1,JL2,JL3,IR
      INTEGER IV1,IV2,IV3,IV4,JR1,JR2,JR3,ISWAP,IA,IB,IEDGE
      REAL*8 XX,YY

      ITOP=0
      MAXSTK=NPL
      NCOUNT=0
      IELM=0
      IL=0
      JL1=0
      JL2=0
      JL3=0
      IR=0
      IV1=0
      IV2=0
      IV3=0
      IV4=0
      JR1=0
      JR2=0
      JR3=0
      ISWAP=0
      IA=0
      IB=0
      IEDGE=0
      XX=0.0D0
      YY=0.0D0

      DO 10 I=1,NTE
        IELM=I
        ITOP=ITOP+1
        JSTACK(ITOP)=IPUSH(IELM,MAXSTK,ITOP,JSTACK,LTJ,IERR)
        IF(IERR.NE.0)RETURN
   10 CONTINUE

   20 IF(ITOP.GT.0)THEN
        NCOUNT=NCOUNT+1
        IF(NCOUNT.GT.LTE)THEN
          IERR = 2141
          RETURN
C          WRITE(*,'('' ***ERROR IN SUBROUTINE LAWSON***'')')
C          WRITE(*,'('' ***NON-CONVERGENCE***'')')
C          PAUSE
C          STOP
        END IF
        IL=JSTACK(ITOP)
        ITOP=ITOP-1
        DO 30 J=1,3
          JL1=J
          JL2=MOD(JL1,3)+1
          JL3=MOD(JL2,3)+1
          IR=JEE(IL,JL1)
          IF(IR.EQ.0)GO TO 30
          IV1=IEN(IL,JL1)
          IV2=IEN(IL,JL2)
          IV3=IEN(IL,JL3)
          XX=PX(IEN(IL,JL3))
          YY=PY(IEN(IL,JL3))
          CALL EDGE(IR,IL,JEE,LTE,JR1,IERR)
          IF (IERR.NE.0) RETURN
          JR2=MOD(JR1,3)+1
          JR3=MOD(JR2,3)+1
          IV4=IEN(IR,JR3)
          CALL SWAP(PX(IV2),PY(IV2),PX(IV1),PY(IV1),PX(IV4),
     &              PY(IV4),XX,YY,ISWAP)
          IF(ISWAP.EQ.1)THEN
            IA=JEE(IL,JL2)
            IB=JEE(IR,JR2)

            IEN(IL,JL2)=IV4
            JEE(IL,JL1)=IB
            JEE(IL,JL2)=IR

            IEN(IR,JR2)=IV3
            JEE(IR,JR1)=IA
            JEE(IR,JR2)=IL

            IF(IA.NE.0)THEN
              CALL EDGE(IA,IL,JEE,LTE,IEDGE,IERR)
              IF (IERR.NE.0) RETURN
              JEE(IA,IEDGE)=IR
            END IF
            IF(IB.NE.0)THEN
              CALL EDGE(IB,IR,JEE,LTE,IEDGE,IERR)
              IF (IERR.NE.0) RETURN
              JEE(IB,IEDGE)=IL
            END IF
            ITOP=ITOP+1
            JSTACK(ITOP)=IPUSH(IL,MAXSTK,ITOP,JSTACK,LTJ,IERR)
            IF(IERR.NE.0)RETURN
            GO TO 20
          END IF
   30   CONTINUE
        GO TO 20
      END IF
      RETURN
      END

      
