C=======================================================================
C     Standalone driver for KBAR (B2, CBAR beam).
C
C     Links the UNMODIFIED nastran/NASTRAN-95/mis/kbar.f with the real
C     gmmatd.f and prints the assembled 12x12 stiffness matrix, which
C     becomes benchmark B2's reference in place of the analytical
C     Euler-Bernoulli matrix previously computed in Python.
C
C     Test case (identical to the one stated in the prompt):
C       node A (0,0,0)  node B (2,0,0)   L = 2
C       E = 200e9  G = 76.923e9  A = 0.01
C       I1 = I2 = 8.333e-6   J = 1.667e-5
C
C     Deliberately neutralised so the comparison isolates the beam
C     formulation itself:
C       pin flags = 0        no DOF released at either end
C       offsets   = 0        elastic axis through the grid points
C       I12       = 0        no product of inertia coupling
C       K1 = K2   = 0        no transverse shear correction, which sends
C                            KBAR down the statement-210/230 branch
C                            R = 12*E*I/L**3, i.e. Euler-Bernoulli
C     These are exactly the features the analytical reference lacks, so
C     leaving them non-zero would compare two different elements.
C
C     KBAR is called once per pivot; each call runs two internal passes
C     emitting (pivot,pivot) then (pivot,far). Two pivots fill the 12x12.
C=======================================================================

      PROGRAM TSTKBAR
      IMPLICIT NONE
C
C     ---- COMMON blocks KBAR expects ---------------------------------
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
C     100 words: KBAR's own EQUIVALENCE (ECPT(71),DP(1)) uses the tail of
C     this block as double-precision scratch, so 42 words is not enough.
      REAL             SECPT(100)
      COMMON /SMA1ET/  SECPT
C
      DOUBLE PRECISION KE(144), KEP(144), DELA(6), DELB(6)
      COMMON /SMA1DP/  KE, KEP, DELA, DELB
C
      INTEGER          MATIDC, MATFLG
      REAL             ELTEMP, STRESS
      COMMON /MATIN /  MATIDC, MATFLG, ELTEMP, STRESS
C
      REAL             E, G, NU, RHO, ALPHA, TSUBO, GSUBE
      REAL             SIGT, SIGC, SIGS
      COMMON /MATOUT/  E, G, NU, RHO, ALPHA, TSUBO, GSUBE,
     1                 SIGT, SIGC, SIGS
C
      REAL             FK
      COMMON /HMTOUT/  FK
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
      MATFLG = 1
      STRESS = 0.0
      ELTEMP = 0.0
      FK     = 0.0
      DODET  = 0
      IDETCK = 0
C
      WRITE (6,*) 'KBAR Test Driver - CBAR beam element (B2)'
      WRITE (6,*) '========================================='
      WRITE (6,*) 'Node A (0,0,0)  Node B (2,0,0)  L = 2'
      WRITE (6,*) 'E = 200e9, G = 76.923e9, A = 0.01'
      WRITE (6,*) 'I1 = I2 = 8.333e-6, J = 1.667e-5'
      WRITE (6,*) 'pin flags = 0, offsets = 0, I12 = 0, K1 = K2 = 0'
      WRITE (6,*) ''
C
C     ================================================================
C     One call per pivot. ECPT rebuilt each time as a precaution.
C     ================================================================
C
      DO 100 IPVT = 1, 2
        CALL SETECP
        NPVT  = 100 + IPVT
        NPVTG = NPVT
        WRITE (6,910) NPVT
  910   FORMAT ('  Calling KBAR with NPVT =', I6)
        CALL KBAR
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
  950 FORMAT ('  row ', I2, ' sum = ', E16.9)
C
      WRITE (6,*) ''
      WRITE (6,*) 'Theory: AE/L = 1.0e9   GJ/L = 6.411532e5'
      WRITE (6,*) '        12EI/L**3 = 2.4999e6'
      WRITE (6,*) '        6EI/L**2  = 2.4999e6  (equal at L=2)'
      WRITE (6,*) '        4EI/L     = 3.3332e6'
      WRITE (6,*) '        2EI/L     = 1.6666e6'
      WRITE (6,*) ''
      WRITE (6,*) 'Computation completed.'
      STOP
      END


C=======================================================================
C     Build the ECPT in /SMA1ET/. Layout documented in kbar.f lines 7-48:
C       1 IELID | 2,3 ISILNO | 4-6 SMALLV | 7 ICSSV | 8,9 IPINFL
C       10-12 ZA | 13-15 ZB | 16 IMATID | 17 A | 18 I1 | 19 I2 | 20 FJ
C       21 NSM | 22 FE | 23-30 C1..G2 | 31 K1 | 32 K2 | 33 I12
C       34 MCSIDA | 35-37 GPA | 38 MCSIDB | 39-41 GPB | 42 TEMPEL
C     INTEGER words: 1, 2, 3, 7, 8, 9, 16, 34, 38.
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
C
C     Reference vector: must not be parallel to the element axis. It
C     defines the plane of principal axis 1.
      SECPT(4)  = 0.0
      SECPT(5)  = 1.0
      SECPT(6)  = 0.0
      IECPT(7)  = 0
C
C     Pin flags off
      IECPT(8)  = 0
      IECPT(9)  = 0
C
C     Offsets zero (words 10-15 already zeroed above)
C
      IECPT(16) = 1
      SECPT(17) = 0.01
      SECPT(18) = 8.333E-6
      SECPT(19) = 8.333E-6
      SECPT(20) = 1.667E-5
      SECPT(21) = 0.0
      SECPT(22) = 0.0
C
C     Stress recovery coefficients C1..G2 (words 23-30) stay zero
C
C     Shear area factors zero -> no shear correction
      SECPT(31) = 0.0
      SECPT(32) = 0.0
C
C     Product of inertia zero
      SECPT(33) = 0.0
C
      IECPT(34) = 0
      SECPT(35) = 0.0
      SECPT(36) = 0.0
      SECPT(37) = 0.0
      IECPT(38) = 0
      SECPT(39) = 2.0
      SECPT(40) = 0.0
      SECPT(41) = 0.0
      SECPT(42) = 0.0
      RETURN
      END
