C=======================================================================
C     Standalone driver for KROD (B1, CROD axial + torsional rod).
C
C     Links the UNMODIFIED nastran/NASTRAN-95/mis/krod.f with the real
C     gmmatd.f and prints the assembled 12x12 stiffness matrix, which
C     becomes benchmark B1's reference in place of the analytical
C     EA/L + GJ/L matrix previously computed in Python.
C
C     Test case (identical to the one stated in the prompt):
C       node A (0,0,0)  node B (2,0,0)
C       E = 200e9  A = 0.01  G = 76.923e9  J = 5e-6
C
C     KROD is called once per pivot grid point. Each call emits two 6x6
C     blocks through SMA1B: (pivot,pivot) then (pivot,far). Two pivots
C     therefore fill the whole 12x12.
C
C     IMPORTANT: KROD swaps IECPT(2) and IECPT(3) in place when the pivot
C     is grid B. Because the ECPT lives in COMMON, it must be rebuilt
C     before every call or the second call sees corrupted connectivity.
C=======================================================================

      PROGRAM TSTKROD
      IMPLICIT NONE
C
C     ---- COMMON blocks KROD expects ---------------------------------
C
      INTEGER          ICOM
      COMMON /BLANK /  ICOM
C
      INTEGER          ISYS
      COMMON /SYSTEM/  ISYS
C
      INTEGER          IFCSTM, IFMPT, IFDIT, IDUM1, IFECPT, IGECPT
      INTEGER          IFGPCT, IGGPCT, IFGEI, IGGEI, IFKGG, IGKGG
      INTEGER          IF4GG, IG4GG, IFGPST, IGGPST, INRW, OUTRW
      INTEGER          CLSNRW, CLSRW, NEOR, EOR, MCBKGG(7), MCB4GG(7)
      COMMON /SMA1IO/  IFCSTM, IFMPT, IFDIT, IDUM1, IFECPT, IGECPT,
     1                 IFGPCT, IGGPCT, IFGEI, IGGEI, IFKGG, IGKGG,
     2                 IF4GG, IG4GG, IFGPST, IGGPST, INRW, OUTRW,
     3                 CLSNRW, CLSRW, NEOR, EOR, MCBKGG, MCB4GG
C
      INTEGER          ICSTM, NCSTM, IGPCT, NGPCT, IPOINT, NPOINT
      INTEGER          I6X6K, N6X6K, I6X64, N6X64
      COMMON /SMA1BK/  ICSTM, NCSTM, IGPCT, NGPCT, IPOINT, NPOINT,
     1                 I6X6K, N6X6K, I6X64, N6X64
C
      INTEGER          IOPT4, K4GGSW, NPVT, LEFT, FROWIC, LROWIC
      INTEGER          NROWSC, TNROWS, JMAX, NLINKS, LINK(10), IDETCK
      INTEGER          DODET, NOGO
      COMMON /SMA1CL/  IOPT4, K4GGSW, NPVT, LEFT, FROWIC, LROWIC,
     1                 NROWSC, TNROWS, JMAX, NLINKS, LINK, IDETCK,
     2                 DODET, NOGO
C
      LOGICAL          HEAT
      COMMON /SMA1HT/  HEAT
C
      REAL             SECPT(100)
      COMMON /SMA1ET/  SECPT
C
      INTEGER          MATIDC, MATFLG
      REAL             ELTEMP, STRESS, SINTH, COSTH
      COMMON /MATIN /  MATIDC, MATFLG, ELTEMP, STRESS, SINTH, COSTH
C
      REAL             E, G, NU, RHO, ALPHA, TSUBO, GSUBE
      REAL             SIGT, SIGC, SIGS
      COMMON /MATOUT/  E, G, NU, RHO, ALPHA, TSUBO, GSUBE,
     1                 SIGT, SIGC, SIGS
C
      REAL             FK
      COMMON /HMTOUT/  FK
C
      DOUBLE PRECISION XX, YY, ZZ, XL, XN(3), DSCL, DSCR, DAMPC
      DOUBLE PRECISION D(18), KE(36), TI(9), DUMDP(227)
      COMMON /SMA1DP/  XX, YY, ZZ, XL, XN, DSCL, DSCR, DAMPC,
     1                 D, KE, TI, DUMDP
C
C     ---- capture area shared with the SMA1B stub --------------------
C
      DOUBLE PRECISION KSTORE(12,12)
      INTEGER          NPVTG
      COMMON /KCAP  /  KSTORE, NPVTG
C
      INTEGER          I, J, IPVT
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
      HEAT   = .FALSE.
      NOGO   = 0
      IOPT4  = 0
      K4GGSW = 0
      IFKGG  = 11
      IF4GG  = 13
      ISYS   = 1024
      ICOM   = 0
      MATFLG = 1
      STRESS = 0.0
      SINTH  = 0.0
      COSTH  = 1.0
      ELTEMP = 0.0
      FK     = 0.0
      DODET  = 0
      IDETCK = 0
C
      WRITE (6,*) 'KROD Test Driver - CROD rod element (B1)'
      WRITE (6,*) '========================================'
      WRITE (6,*) 'Node A (0,0,0)  Node B (2,0,0)'
      WRITE (6,*) 'E = 200e9, A = 0.01, G = 76.923e9, J = 5e-6'
      WRITE (6,*) ''
C
C     ================================================================
C     One call per pivot. ECPT rebuilt each time (KROD swaps words 2/3).
C     ================================================================
C
      DO 100 IPVT = 1, 2
        CALL SETECP
        NPVT  = 100 + IPVT
        NPVTG = NPVT
        WRITE (6,910) NPVT
  910   FORMAT ('  Calling KROD with NPVT =', I6)
        CALL KROD
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
      WRITE (6,*) ''
      WRITE (6,*) 'Checks:'
      ASYM = 0.0D0
      DO 500 I = 1, 12
        DO 490 J = I+1, 12
          IF (DABS(KSTORE(I,J)-KSTORE(J,I)) .GT. ASYM)
     1        ASYM = DABS(KSTORE(I,J)-KSTORE(J,I))
  490   CONTINUE
  500 CONTINUE
      WRITE (6,940) ASYM
  940 FORMAT ('  max asymmetry (12x12) = ', E16.9)
C
      DO 600 I = 1, 12
        ROWSUM = 0.0D0
        DO 590 J = 1, 12
          ROWSUM = ROWSUM + KSTORE(I,J)
  590   CONTINUE
        WRITE (6,950) I, ROWSUM
  600 CONTINUE
  950 FORMAT ('  row ', I2, ' sum (rigid body) = ', E16.9)
C
      WRITE (6,*) ''
      WRITE (6,*) 'Theory: A*E/L = 1.0e9, J*G/L = 1.923075e5'
      WRITE (6,*) ''
      WRITE (6,*) 'Computation completed.'
      STOP
      END


C=======================================================================
C     Build the ECPT in /SMA1ET/. Layout documented in krod.f lines 6-24.
C     Words 1,2,3,4,9,13 are INTEGER; the rest REAL.
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
      IECPT(4)  = 1
      SECPT(5)  = 0.01
      SECPT(6)  = 5.0E-6
      SECPT(7)  = 0.0
      SECPT(8)  = 0.0
      IECPT(9)  = 0
      SECPT(10) = 0.0
      SECPT(11) = 0.0
      SECPT(12) = 0.0
      IECPT(13) = 0
      SECPT(14) = 2.0
      SECPT(15) = 0.0
      SECPT(16) = 0.0
      SECPT(17) = 0.0
      RETURN
      END
