C=======================================================================
C     Standalone driver for KTRMEM (B3, CTRIA3 triangular membrane).
C
C     Links the UNMODIFIED nastran/NASTRAN-95/mis/ktrmem.f together with
C     the real gmmatd.f, and prints the assembled stiffness matrix. This
C     output becomes benchmark B3's reference, replacing the analytical
C     CST matrix that was previously computed in Python.
C
C     Test case (identical to the one stated in the prompt):
C       nodes (0,0), (2,0), (1,1.5)   E = 2.1e11   nu = 0.3   t = 0.01
C
C     KTRMEM is called once per pivot grid point. Each call fills the row
C     of 3x3 blocks belonging to that pivot, via SMA1B. The stubs
C     accumulate them into a 9x9 (3 nodes x 3 translational DOF).
C
C     Note on dimensions: NASTRAN works in 3 translational DOF per node,
C     so the natural output is 9x9 with the z rows/columns identically
C     zero for a flat membrane. The 6x6 in-plane sub-matrix (ux,uy per
C     node) is printed separately, since that is the form the benchmark
C     compares.
C=======================================================================

      PROGRAM TSTKTRMEM
      IMPLICIT NONE
C
C     ---- COMMON blocks KTRMEM expects -------------------------------
C
      INTEGER          MECPT(1), NGRID(3), MATID1
      REAL             ANGLE, T, FMU, DUMMY1
      REAL             X1, Y1, Z1, DUMMY2, X2, Y2, Z2, DUMMY3
      REAL             X3, Y3, Z3, DUMB(80)
      COMMON /SMA1ET/  MECPT, NGRID, ANGLE, MATID1, T, FMU, DUMMY1,
     1                 X1, Y1, Z1, DUMMY2, X2, Y2, Z2, DUMMY3,
     2                 X3, Y3, Z3, DUMB
C
      DOUBLE PRECISION KIJ(36), C(18), E(9), TEMPAR(27), TI(9), TEMP
      DOUBLE PRECISION XSUBB, XSUBC, YSUBC, VOL, REELMU, DELTA, FLAMDA
      DOUBLE PRECISION THETA
      INTEGER          KA, NPOINT, NSAVE, DUMMY(382)
      COMMON /SMA1DP/  KIJ, C, E, TEMPAR, TI, TEMP,
     1                 XSUBB, XSUBC, YSUBC, VOL, REELMU, DELTA, FLAMDA,
     2                 THETA, KA, NPOINT, NSAVE, DUMMY
C
      LOGICAL          HEAT
      COMMON /SMA1HT/  HEAT
C
      INTEGER          MATID, INFLAG
      REAL             ELTEMP, STRESS, SINTH, COSTH
      COMMON /MATIN /  MATID, INFLAG, ELTEMP, STRESS, SINTH, COSTH
C
      REAL             G11, G12, G13, G22, G23, G33, RHO
      REAL             ALPHA1, ALPHA2, ALP12, TSUB0, GSUBE
      REAL             SIGTEN, SIGCOM, SIGSHE
      REAL             G2X211, G2X212, G2X222
      COMMON /MATOUT/  G11, G12, G13, G22, G23, G33, RHO,
     1                 ALPHA1, ALPHA2, ALP12, TSUB0, GSUBE,
     2                 SIGTEN, SIGCOM, SIGSHE,
     3                 G2X211, G2X212, G2X222
C
      INTEGER          IOPT4, K4GGSW, NPVT, ISKP, NOGOO
      COMMON /SMA1CL/  IOPT4, K4GGSW, NPVT, ISKP, NOGOO
C
      INTEGER          DUM1, IFKGG, DUM2, IF4GG, DUM3
      COMMON /SMA1IO/  DUM1, IFKGG, DUM2, IF4GG, DUM3
C
      INTEGER          SYSBUF, IOUT, NOGO
      COMMON /SYSTEM/  SYSBUF, IOUT, NOGO
C
      REAL             CONSTS(5)
      COMMON /CONDAS/  CONSTS
C
      REAL             MATBUF(4)
      COMMON /HMTOUT/  MATBUF
C
C     ---- capture area shared with the SMA1B stub --------------------
C
      DOUBLE PRECISION KSTORE(9,9)
      INTEGER          NPVTG
      COMMON /KCAP  /  KSTORE, NPVTG
C
C     ---- ECPT overlay ----------------------------------------------
C
      REAL             ECPT(21)
      INTEGER          NECPT(21)
      EQUIVALENCE      (ECPT(1), NECPT(1))
C
      INTEGER          I, J, IPVT, NTYPE
      DOUBLE PRECISION KSMALL(6,6)
      DOUBLE PRECISION ROWSUM, ASYM
C
C     ================================================================
C     Initialise everything KTRMEM reads
C     ================================================================
C
      DO 5 I = 1, 9
        DO 4 J = 1, 9
          KSTORE(I,J) = 0.0D0
    4   CONTINUE
    5 CONTINUE
C
      HEAT   = .FALSE.
      NOGO   = 0
      NOGOO  = 0
      IOPT4  = 0
      K4GGSW = 0
      ISKP   = 0
      IFKGG  = 1
      IF4GG  = 2
      SYSBUF = 1024
      IOUT   = 6
      INFLAG = 2
      STRESS = 0.0
      SINTH  = 0.0
      COSTH  = 1.0
      ELTEMP = 0.0
C
C     CONSTS(4) is DEGRA (degrees to radians) via EQUIVALENCE in KTRMEM
      CONSTS(1) = 0.0
      CONSTS(2) = 0.0
      CONSTS(3) = 0.0
      CONSTS(4) = 3.1415926535897932 / 180.0
      CONSTS(5) = 0.0
C
      DO 6 I = 1, 4
        MATBUF(I) = 0.0
    6 CONTINUE
C
C     ---- element connection and property table ----------------------
C     Layout is documented in ktrmem.f lines 40-62.
C
      NECPT(1)  = 1001
      NECPT(2)  = 101
      NECPT(3)  = 102
      NECPT(4)  = 103
      ECPT(5)   = 0.0
      NECPT(6)  = 1
      ECPT(7)   = 0.01
      ECPT(8)   = 0.0
      NECPT(9)  = 0
      ECPT(10)  = 0.0
      ECPT(11)  = 0.0
      ECPT(12)  = 0.0
      NECPT(13) = 0
      ECPT(14)  = 2.0
      ECPT(15)  = 0.0
      ECPT(16)  = 0.0
      NECPT(17) = 0
      ECPT(18)  = 1.0
      ECPT(19)  = 1.5
      ECPT(20)  = 0.0
      ECPT(21)  = 0.0
C
C     Copy into /SMA1ET/, which is what KTRMEM actually reads
      DO 8 I = 1, 21
        CALL PUTECP (I, ECPT(I), NECPT(I))
    8 CONTINUE
C
      WRITE (6,*) 'KTRMEM Test Driver - CTRIA3 membrane (B3)'
      WRITE (6,*) '=========================================='
      WRITE (6,*) 'Nodes: (0,0) (2,0) (1,1.5)'
      WRITE (6,*) 'E = 2.1e11, nu = 0.3, t = 0.01'
      WRITE (6,*) ''
C
C     ================================================================
C     One call per pivot grid point
C     ================================================================
C
      NTYPE = 0
      DO 100 IPVT = 1, 3
        NPVT  = 100 + IPVT
        NPVTG = NPVT
        WRITE (6,910) NPVT
  910   FORMAT ('  Calling KTRMEM with NPVT =', I6)
        CALL KTRMEM (NTYPE)
        IF (NOGO .NE. 0) THEN
          WRITE (6,*) '  ERROR: NOGO set after pivot', IPVT
        END IF
  100 CONTINUE
C
C     ================================================================
C     Report
C     ================================================================
C
      WRITE (6,*) ''
      WRITE (6,*) 'BEGIN_FULL_9X9'
      DO 200 I = 1, 9
        WRITE (6,920) (KSTORE(I,J), J = 1, 9)
  200 CONTINUE
  920 FORMAT (9(1X,E16.9))
      WRITE (6,*) 'END_FULL_9X9'
C
C     In-plane 6x6: drop the z DOF of each node (rows/cols 3,6,9)
C
      DO 320 I = 1, 3
        DO 310 J = 1, 3
          KSMALL(2*I-1, 2*J-1) = KSTORE(3*I-2, 3*J-2)
          KSMALL(2*I-1, 2*J  ) = KSTORE(3*I-2, 3*J-1)
          KSMALL(2*I  , 2*J-1) = KSTORE(3*I-1, 3*J-2)
          KSMALL(2*I  , 2*J  ) = KSTORE(3*I-1, 3*J-1)
  310   CONTINUE
  320 CONTINUE
C
      WRITE (6,*) ''
      WRITE (6,*) 'BEGIN_INPLANE_6X6'
      DO 400 I = 1, 6
        WRITE (6,930) (KSMALL(I,J), J = 1, 6)
  400 CONTINUE
  930 FORMAT (6(1X,E16.9))
      WRITE (6,*) 'END_INPLANE_6X6'
C
C     Sanity checks
C
      WRITE (6,*) ''
      WRITE (6,*) 'Checks:'
      ASYM = 0.0D0
      DO 500 I = 1, 6
        DO 490 J = I+1, 6
          IF (DABS(KSMALL(I,J)-KSMALL(J,I)) .GT. ASYM)
     1        ASYM = DABS(KSMALL(I,J)-KSMALL(J,I))
  490   CONTINUE
  500 CONTINUE
      WRITE (6,940) ASYM
  940 FORMAT ('  max asymmetry (6x6) = ', E16.9)
C
      DO 600 I = 1, 6
        ROWSUM = 0.0D0
        DO 590 J = 1, 6
          ROWSUM = ROWSUM + KSMALL(I,J)
  590   CONTINUE
        WRITE (6,950) I, ROWSUM
  600 CONTINUE
  950 FORMAT ('  row ', I2, ' sum (rigid body) = ', E16.9)
C
      WRITE (6,*) ''
      WRITE (6,*) 'Computation completed.'
      STOP
      END


C=======================================================================
C     Copy one ECPT word into /SMA1ET/. Done through a subroutine so the
C     REAL/INTEGER duality of the ECPT is preserved without relying on
C     an EQUIVALENCE across COMMON boundaries.
C=======================================================================
      SUBROUTINE PUTECP (IDX, RVAL, IVAL)
      IMPLICIT NONE
      INTEGER IDX, IVAL
      REAL    RVAL
C
      REAL             SBUF(101)
      INTEGER          IBUF(101)
      COMMON /SMA1ET/  SBUF
      EQUIVALENCE      (SBUF(1), IBUF(1))
C
C     Words 1,2,3,4,6,9,13,17 are integers in the ECPT; the rest real.
      IF (IDX .EQ. 1 .OR. IDX .EQ. 2 .OR. IDX .EQ. 3 .OR.
     1    IDX .EQ. 4 .OR. IDX .EQ. 6 .OR. IDX .EQ. 9 .OR.
     2    IDX .EQ. 13 .OR. IDX .EQ. 17) THEN
        IBUF(IDX) = IVAL
      ELSE
        SBUF(IDX) = RVAL
      END IF
      RETURN
      END
