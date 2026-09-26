C=======================================================================
C     Standalone driver for KQDMEM (B4, CQUAD4 quadrilateral membrane).
C
C     Links the UNMODIFIED kqdmem.f and ktrmem.f with the real gmmatd.f
C     and prints the assembled stiffness matrix, which becomes benchmark
C     B4's reference in place of the analytical two-CST assembly
C     previously computed in Python.
C
C     Test case (identical to the one stated in the prompt):
C       nodes (0,0) (2,0) (2,1.5) (0,1.5)
C       E = 200e9   nu = 0.3   t = 0.01
C
C     Why this matters for B4 specifically. KQDMEM does NOT split the
C     quad into two triangles the way the analytical reference does. It
C     uses the map
C         DATA M / 1,2,4,  2,3,1,  3,4,2,  4,1,3 /
C     i.e. FOUR overlapping corner triangles, halves the thickness
C         ECPT(7) = ECPT(8)/2.0
C     and for each pivot skips the one triangle that does not contain it
C     (JNOT), so THREE sub-triangles at t/2 contribute per pivot.
C     The analytical reference sums TWO triangles at full t. These are
C     different elements, which is why B4's residuals clustered on
C     discrete rational ratios rather than looking like noise.
C
C     KQDMEM is called once per pivot. Each call emits four 3x3 blocks
C     through SMA1B, so four pivots fill the 12x12.
C=======================================================================

      PROGRAM TSTKQDMEM
      IMPLICIT NONE
C
C     ---- COMMON blocks KQDMEM and KTRMEM expect ---------------------
C
      REAL             CONSTS(5)
      COMMON /CONDAS/  CONSTS
C
      LOGICAL          HEAT
      COMMON /SMA1HT/  HEAT
C
      REAL             SECPT(100)
      COMMON /SMA1ET/  SECPT
C
      INTEGER          DUM1(10), IFKGG, DUM2(1), IF4GG, DUM3(23)
      COMMON /SMA1IO/  DUM1, IFKGG, DUM2, IF4GG, DUM3
C
      INTEGER          IOPT4, K4GGSW, NPVT, DUMCL(7), LINK(10)
      INTEGER          IDETCK, DODET, NOGO
      COMMON /SMA1CL/  IOPT4, K4GGSW, NPVT, DUMCL, LINK, IDETCK,
     1                 DODET, NOGO
C
C     Sized to cover KQDMEM's own view of this block, which is the
C     largest of the participating routines (~4016 bytes).
      DOUBLE PRECISION DPBLK(600)
      COMMON /SMA1DP/  DPBLK
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
      REAL             MATBUF(7)
      COMMON /HMTOUT/  MATBUF
C
      INTEGER          SYSBUF, IOUT, NOGOS
      COMMON /SYSTEM/  SYSBUF, IOUT, NOGOS
C
      REAL             BSKIP, VOLUME, SURFAC
      COMMON /BLANK /  BSKIP, VOLUME, SURFAC
C
      LOGICAL          HYDRO
      COMMON /HYDROE/  HYDRO
C
C     ---- capture area shared with the SMA1B stub --------------------
C
      DOUBLE PRECISION KSTORE(12,12)
      INTEGER          NPVTG
      COMMON /KCAP  /  KSTORE, NPVTG
C
      INTEGER          I, J, IPVT
      DOUBLE PRECISION KSMALL(8,8)
      DOUBLE PRECISION ASYM, ROWSUM
C
C     ================================================================
C     Initialise
C     ================================================================
C
      DO 5 I = 1, 12
        DO 4 J = 1, 12
          KSTORE(I,J) = 0.0D0
    4   CONTINUE
    5 CONTINUE
C
      DO 6 I = 1, 600
        DPBLK(I) = 0.0D0
    6 CONTINUE
C
      HEAT   = .FALSE.
      HYDRO  = .FALSE.
      NOGO   = 0
      NOGOS  = 0
      IOPT4  = 0
      K4GGSW = 0
      IFKGG  = 1
      IF4GG  = 2
      SYSBUF = 1024
      IOUT   = 6
      INFLAG = 2
      STRESS = 0.0
      SINTH  = 0.0
      COSTH  = 1.0
      ELTEMP = 0.0
      IDETCK = 0
      DODET  = 0
      BSKIP  = 0.0
      VOLUME = 0.0
      SURFAC = 0.0
C
C     CONSTS(4) is DEGRA, reached through EQUIVALENCE in both routines
      CONSTS(1) = 0.0
      CONSTS(2) = 0.0
      CONSTS(3) = 0.0
      CONSTS(4) = 3.1415926535897932 / 180.0
      CONSTS(5) = 0.0
C
      DO 7 I = 1, 7
        MATBUF(I) = 0.0
    7 CONTINUE
C
      WRITE (6,*) 'KQDMEM Test Driver - CQUAD4 membrane (B4)'
      WRITE (6,*) '========================================='
      WRITE (6,*) 'Nodes: (0,0) (2,0) (2,1.5) (0,1.5)'
      WRITE (6,*) 'E = 200e9, nu = 0.3, t = 0.01'
      WRITE (6,*) 'NASTRAN: 4 overlapping corner triangles at t/2,'
      WRITE (6,*) '         3 of them per pivot (JNOT skips one).'
      WRITE (6,*) ''
C
C     ================================================================
C     One call per pivot. The ECPT MUST be rebuilt every time: KQDMEM
C     overwrites words 5..8 and 21 for KTRMEM's benefit and rewrites
C     NECPT(2..4) plus the coordinate slots once per sub-triangle.
C     ================================================================
C
      DO 100 IPVT = 1, 4
        CALL SETECP
        NPVT  = 100 + IPVT
        NPVTG = NPVT
        WRITE (6,910) NPVT
  910   FORMAT ('  Calling KQDMEM with NPVT =', I6)
        CALL KQDMEM
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
      WRITE (6,*) 'BEGIN_FULL_12X12'
      DO 200 I = 1, 12
        WRITE (6,920) (KSTORE(I,J), J = 1, 12)
  200 CONTINUE
  920 FORMAT (12(1X,E16.9))
      WRITE (6,*) 'END_FULL_12X12'
C
C     In-plane 8x8: drop the z DOF of each node (rows/cols 3,6,9,12)
C
      DO 320 I = 1, 4
        DO 310 J = 1, 4
          KSMALL(2*I-1, 2*J-1) = KSTORE(3*I-2, 3*J-2)
          KSMALL(2*I-1, 2*J  ) = KSTORE(3*I-2, 3*J-1)
          KSMALL(2*I  , 2*J-1) = KSTORE(3*I-1, 3*J-2)
          KSMALL(2*I  , 2*J  ) = KSTORE(3*I-1, 3*J-1)
  310   CONTINUE
  320 CONTINUE
C
      WRITE (6,*) ''
      WRITE (6,*) 'BEGIN_FULL_8X8'
      DO 400 I = 1, 8
        WRITE (6,930) (KSMALL(I,J), J = 1, 8)
  400 CONTINUE
  930 FORMAT (8(1X,E16.9))
      WRITE (6,*) 'END_FULL_8X8'
C
      WRITE (6,*) ''
      WRITE (6,*) 'Checks:'
      ASYM = 0.0D0
      DO 500 I = 1, 8
        DO 490 J = I+1, 8
          IF (DABS(KSMALL(I,J)-KSMALL(J,I)) .GT. ASYM)
     1        ASYM = DABS(KSMALL(I,J)-KSMALL(J,I))
  490   CONTINUE
  500 CONTINUE
      WRITE (6,940) ASYM
  940 FORMAT ('  max asymmetry (8x8) = ', E16.9)
C
      DO 600 I = 1, 8
        ROWSUM = 0.0D0
        DO 590 J = 1, 8
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
C     Build the ECPT in /SMA1ET/. Layout documented in kqdmem.f lines
C     35-62 (the left-hand column, "ECPT RECEIVED BY KQDMEM"):
C       1 EL ID | 2-5 grid A,B,C,D | 6 THETA | 7 MATID | 8 T | 9 NSM
C       10 CSID1 | 11-13 X1,Y1,Z1 | 14 CSID2 | 15-17 X2,Y2,Z2
C       18 CSID3 | 19-21 X3,Y3,Z3 | 22 CSID4 | 23-25 X4,Y4,Z4
C       26 ELTEMP
C     INTEGER words: 1, 2, 3, 4, 5, 7, 10, 14, 18, 22.
C
C     Note KQDMEM tests IF (NECPT(8) .EQ. 1) to detect a TRAPRG heat
C     ring. Word 8 holds T = 0.01 as a REAL, whose bit pattern is not
C     the integer 1, so HRING stays .FALSE. as intended.
C=======================================================================
      SUBROUTINE SETECP
      IMPLICIT NONE
C
      REAL             SECPT(100)
      INTEGER          IECPT(100)
      COMMON /SMA1ET/  SECPT
      EQUIVALENCE      (SECPT(1), IECPT(1))
C
      INTEGER I
C
      DO 10 I = 1, 100
        SECPT(I) = 0.0
   10 CONTINUE
C
      IECPT(1)  = 1001
      IECPT(2)  = 101
      IECPT(3)  = 102
      IECPT(4)  = 103
      IECPT(5)  = 104
      SECPT(6)  = 0.0
      IECPT(7)  = 1
      SECPT(8)  = 0.01
      SECPT(9)  = 0.0
C
      IECPT(10) = 0
      SECPT(11) = 0.0
      SECPT(12) = 0.0
      SECPT(13) = 0.0
C
      IECPT(14) = 0
      SECPT(15) = 2.0
      SECPT(16) = 0.0
      SECPT(17) = 0.0
C
      IECPT(18) = 0
      SECPT(19) = 2.0
      SECPT(20) = 1.5
      SECPT(21) = 0.0
C
      IECPT(22) = 0
      SECPT(23) = 0.0
      SECPT(24) = 1.5
      SECPT(25) = 0.0
C
      SECPT(26) = 0.0
      RETURN
      END
