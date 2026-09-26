C     STANDALONE TEST FOR KTRBSC - Basic Bending Triangle
C     Uses IOPT=0 so KTRBSC does full computation including coord setup
C     ====================================================================
      PROGRAM TESTKTRBSC
      IMPLICIT NONE
C
C     COMMON blocks - EXACTLY matching test_ktrplt.f (which works)
C
      REAL             ECPT(100)
      INTEGER          NECPT(100)
      EQUIVALENCE     (NECPT(1), ECPT(1))
      COMMON /SMA1ET/ ECPT
C
      REAL             CONSTS(5)
      COMMON /CONDAS/ CONSTS
C
      INTEGER          MATID, INFLAG
      REAL             ELTEMP, STRESS, SINTH, COSTH
      COMMON /MATIN / MATID, INFLAG, ELTEMP, STRESS, SINTH, COSTH
C
      REAL             G11, G12, G13, G22, G23, G33, RHO
      REAL             ALPHA1, ALPHA2, ALP12, TSUB0, GSUBE
      REAL             SIGTEN, SIGCOM, SIGSHE
      REAL             G2X211, G2X212, G2X222
      COMMON /MATOUT/ G11, G12, G13, G22, G23, G33, RHO,
     1                ALPHA1, ALPHA2, ALP12, TSUB0, GSUBE,
     2                SIGTEN, SIGCOM, SIGSHE,
     3                G2X211, G2X212, G2X222
C
      REAL             DUM1(10)
      INTEGER          IFKGG
      REAL             DUM2(1)
      INTEGER          IF4GG
      REAL             DUM3(23)
      COMMON /SMA1IO/ DUM1, IFKGG, DUM2, IF4GG, DUM3
C
      INTEGER          IOPT4, K4GGSW, NPVT
      REAL             DUMCL(7)
      INTEGER          LINK(10), IDETCK
      REAL             DODET
      INTEGER          NOGO
      COMMON /SMA1CL/ IOPT4, K4GGSW, NPVT, DUMCL, LINK, IDETCK,
     1                DODET, NOGO
C
      DOUBLE PRECISION DPWORK(400)
      COMMON /SMA1DP/ DPWORK
C
      INTEGER          MACH
      COMMON /MACHIN/ MACH
C
C     Capture buffer for SMA1B (same as test_ktrplt.f)
      DOUBLE PRECISION KSTORE(6,6,4,4)
      INTEGER          CURPVT, NPVTS
      COMMON /CAPTUREP/ KSTORE, CURPVT, NPVTS
C
C     Local
      INTEGER          I, J, BI, BJ, BLK
      DOUBLE PRECISION KSU(9,9)
C
C     Initialize (same as test_ktrplt.f)
      MACH = 1
      CONSTS(1) = 3.14159265358979
      CONSTS(2) = 0.0
      CONSTS(3) = 0.0
      CONSTS(4) = 3.14159265358979 / 180.0
      CONSTS(5) = 0.0
      IOPT4 = 0
      K4GGSW = 0
      NOGO = 0
      IDETCK = 0
      DODET = 0.0
      IFKGG = 0
      IF4GG = 0
      NPVTS = 3
      GSUBE = 0.0
      CURPVT = 1
C
      DO 5 I = 1, 100
        ECPT(I) = 0.0
    5 CONTINUE
      DO 6 I = 1, 400
        DPWORK(I) = 0.0D0
    6 CONTINUE
      DO 7 I = 1, 4
        DO 7 J = 1, 4
          DO 7 BI = 1, 6
            DO 7 BJ = 1, 6
              KSTORE(BI,BJ,I,J) = 0.0D0
    7 CONTINUE
C
C     ECPT for KTRBSC:
C     Sub-triangle with A=(1/3,1/3,0) B=(0,0,0) C=(1,0,0)
C     This is sub-triangle 1 of the right-angle plate (0,0)-(1,0)-(0,1)
      NECPT(1) = 1
      NECPT(2) = 101
      NECPT(3) = 102
      NECPT(4) = 103
      ECPT(5) = 0.0
      NECPT(6) = 1
      ECPT(7) = 8.3333333E-8
      NECPT(8) = 0
      ECPT(9) = 0.0
      ECPT(10) = 0.0
      ECPT(11) = -0.005
      ECPT(12) = 0.005
C     Node A: centroid (1/3, 1/3, 0)
      NECPT(13) = 0
      ECPT(14) = 0.333333
      ECPT(15) = 0.333333
      ECPT(16) = 0.0
C     Node B: (0, 0, 0)
      NECPT(17) = 0
      ECPT(18) = 0.0
      ECPT(19) = 0.0
      ECPT(20) = 0.0
C     Node C: (1, 0, 0)
      NECPT(21) = 0
      ECPT(22) = 1.0
      ECPT(23) = 0.0
      ECPT(24) = 0.0
      ECPT(25) = 20.0
C
C     Set NPVT to node A so KTRBSC knows pivot
      NPVT = 101
C
C     Pre-compute XSUBB, XSUBC, YSUBC (local coords) and place in SMA1DP
C     For IOPT=1, KTRBSC expects these pre-set by the caller.
C     SMA1DP layout: A(225)+PROD9(9)+TEMP9(9) = 243, then:
C       XSUBB=DPWORK(244), XSUBC=DPWORK(245), YSUBC=DPWORK(246)
C     A=(1/3,1/3,0) B=(0,0,0) C=(1,0,0):
C       XSUBB=|B-A|=sqrt(2)/3=0.471405
C       i=(B-A)/XSUBB=(-0.707107,-0.707107,0)
C       XSUBC=i.(C-A)=-0.235702
C       YSUBC=|i x (C-A)|=0.707107
      DPWORK(244) = 0.4714045D0
      DPWORK(245) = -0.2357023D0
      DPWORK(246) = 0.7071068D0
C
      WRITE(6,*) '=== KTRBSC Standalone Test (IOPT=1) ==='
      WRITE(6,*) 'A=(0.333,0.333,0) B=(0,0,0) C=(1,0,0)'
      WRITE(6,*) 'XSUBB=0.4714 XSUBC=-0.2357 YSUBC=0.7071'
      WRITE(6,*) 'E=200e9, NU=0.3, I=8.333e-8'
C
C     Call KTRBSC with IOPT=1 (9x9 K-super-U only)
      CALL KTRBSC(1)
C
      IF (NOGO .EQ. 1) THEN
        WRITE(6,*) 'ERROR: KTRBSC failed (NOGO=1)'
        STOP
      ENDIF
C
C     Extract 9x9 from DPWORK(1:81)
C     9 blocks of 3x3, each stored by rows
      DO 40 BI = 1, 3
        DO 40 BJ = 1, 3
          BLK = 9*(3*(BI-1) + (BJ-1))
          DO 35 I = 1, 3
            DO 35 J = 1, 3
              KSU(3*(BI-1)+I, 3*(BJ-1)+J) = DPWORK(BLK+(I-1)*3+J)
   35     CONTINUE
   40 CONTINUE
C
      WRITE(6,*) ''
      WRITE(6,*) '9x9 K-super-U:'
      DO 50 I = 1, 9
        WRITE(6,200) (KSU(I,J), J=1,9)
   50 CONTINUE
  200 FORMAT(9(1X,E13.6))
C
      WRITE(6,*) ''
      WRITE(6,*) '=== DONE ==='
      END
