C     TEST DRIVER FOR KSOLID - WEDGE ELEMENT (6 nodes, ITYPE=1)
C     ===========================================================
C     This program calls KSOLID(1) which decomposes a wedge into
C     3 tetrahedra and calls KTETRA for each. The full 18x18 stiffness
C     matrix is assembled from the sub-block outputs.
C
C     Test geometry (triangular prism / wedge):
C       Node 1: (0, 0, 0)    (base triangle)
C       Node 2: (1, 0, 0)
C       Node 3: (0, 1, 0)
C       Node 4: (0, 0, 1)    (top triangle)
C       Node 5: (1, 0, 1)
C       Node 6: (0, 1, 1)
C
C     Material: E=200e9, NU=0.3
C
      PROGRAM TESTWEDGE
      IMPLICIT NONE
C
C     COMMON blocks matching KSOLID/KTETRA expectations
C
      REAL             ECPT(100)
      INTEGER          NECPT(100)
      EQUIVALENCE     (NECPT(1), ECPT(1))
      COMMON /SMA1ET/ ECPT
C
      DOUBLE PRECISION C(72), G(36), H(16), TEMP(12), T(9)
      DOUBLE PRECISION CT(18), GCT(18), KIJ(36), HDETER, TEMP1
      INTEGER          NGPT, DIREC, KOUNT
      DOUBLE PRECISION TVOL
      COMMON /SMA1DP/ C, G, H, TEMP, T, CT, GCT, KIJ, HDETER, TEMP1,
     1                NGPT, DIREC, KOUNT, TVOL
C
      LOGICAL          HEAT
      COMMON /SMA1HT/ HEAT
C
      LOGICAL          HYDRO
      COMMON /HYDROE/ HYDRO
C
      INTEGER          MATID, INFLAG
      REAL             ELTEMP
      COMMON /MATIN / MATID, INFLAG, ELTEMP
C
      REAL             E, GG, NU, RHO, ALPHA, TSUB0, GSUBE
      REAL             SIGT, SIGC, SIGS
      COMMON /MATOUT/ E, GG, NU, RHO, ALPHA, TSUB0, GSUBE,
     1                SIGT, SIGC, SIGS
C
      REAL             MATBUF(7)
      COMMON /HMTOUT/ MATBUF
C
      INTEGER          IOPT4, K4GGSW, NPVT, ISKP(19)
      LOGICAL          NOGOO
      COMMON /SMA1CL/ IOPT4, K4GGSW, NPVT, ISKP, NOGOO
C
      REAL             DUM1(10), IFKGG, DUM2(1), IF4GG, DUM3(23)
      COMMON /SMA1IO/ DUM1, IFKGG, DUM2, IF4GG, DUM3
C
      INTEGER          SYSBUF, OUT
      LOGICAL          NOGO
      COMMON /SYSTEM/ SYSBUF, OUT, NOGO
C
      INTEGER          MACH
      COMMON /MACHIN/ MACH
C
      REAL             BSKIP(16), VOLUME, SURFAC
      COMMON /BLANK / BSKIP, VOLUME, SURFAC
C
      CHARACTER        UFM*23
      COMMON /XMSSG / UFM
C
C     Local variables for assembling the 18x18 matrix
C
      DOUBLE PRECISION KFULL(18,18)
      INTEGER          I, J, IPVT
      DOUBLE PRECISION ROWSUM, TRACE
      REAL             ECPT_SAVE(100)
C
C     Initialize system commons
C
      UFM = ' *** USER FATAL MESSAGE'
      MACH = 1
      OUT = 6
      SYSBUF = 1
      NOGO = .FALSE.
      HEAT = .FALSE.
      HYDRO = .FALSE.
      NOGOO = .FALSE.
      IOPT4 = 0
      K4GGSW = 0
      VOLUME = 0.0
      SURFAC = 0.0
C
C     Zero out SMA1IO
C
      DO 5 I = 1, 10
        DUM1(I) = 0.0
    5 CONTINUE
      IFKGG = 0.0
      DUM2(1) = 0.0
      IF4GG = 0.0
      DO 6 I = 1, 23
        DUM3(I) = 0.0
    6 CONTINUE
C
C     Set material properties
C
      E = 200.0E9
      NU = 0.3
      GG = E / (2.0*(1.0+NU))
      RHO = 7800.0
      ALPHA = 0.0
      TSUB0 = 0.0
      GSUBE = 0.0
      SIGT = 0.0
      SIGC = 0.0
      SIGS = 0.0
C
C     Fill ECPT for wedge (ITYPE=1, 6 nodes)
C     Layout from ksolid.f comments:
C       ECPT(1) = EL ID
C       ECPT(2) = MAT-ID
C       ECPT(3-8) = GRID-1 through GRID-6 (SIL IDs)
C       ECPT(9)  = CSID-1
C       ECPT(10-12) = X1, Y1, Z1
C       ECPT(13) = CSID-2
C       ECPT(14-16) = X2, Y2, Z2
C       ECPT(17) = CSID-3
C       ECPT(18-20) = X3, Y3, Z3
C       ECPT(21) = CSID-4
C       ECPT(22-24) = X4, Y4, Z4
C       ECPT(25) = CSID-5
C       ECPT(26-28) = X5, Y5, Z5
C       ECPT(29) = CSID-6
C       ECPT(30-32) = X6, Y6, Z6
C       ECPT(33) = ELTEMP
C
C     Zero ECPT first
      DO 10 I = 1, 100
        ECPT(I) = 0.0
   10 CONTINUE
C
C     Element ID and Material ID
      NECPT(1) = 1
      NECPT(2) = 1
C
C     Grid SIL IDs (nodes 1-6)
      NECPT(3) = 101
      NECPT(4) = 102
      NECPT(5) = 103
      NECPT(6) = 104
      NECPT(7) = 105
      NECPT(8) = 106
C
C     Node 1: (0,0,0)
      NECPT(9) = 0
      ECPT(10) = 0.0
      ECPT(11) = 0.0
      ECPT(12) = 0.0
C
C     Node 2: (1,0,0)
      NECPT(13) = 0
      ECPT(14) = 1.0
      ECPT(15) = 0.0
      ECPT(16) = 0.0
C
C     Node 3: (0,1,0)
      NECPT(17) = 0
      ECPT(18) = 0.0
      ECPT(19) = 1.0
      ECPT(20) = 0.0
C
C     Node 4: (0,0,1)
      NECPT(21) = 0
      ECPT(22) = 0.0
      ECPT(23) = 0.0
      ECPT(24) = 1.0
C
C     Node 5: (1,0,1)
      NECPT(25) = 0
      ECPT(26) = 1.0
      ECPT(27) = 0.0
      ECPT(28) = 1.0
C
C     Node 6: (0,1,1)
      NECPT(29) = 0
      ECPT(30) = 0.0
      ECPT(31) = 1.0
      ECPT(32) = 1.0
C
C     Element temperature
      ECPT(33) = 20.0
C
C     Zero the full stiffness matrix
C
      DO 20 I = 1, 18
        DO 15 J = 1, 18
          KFULL(I,J) = 0.0D0
   15   CONTINUE
   20 CONTINUE
C
C     Save the ECPT (KSOLID modifies it internally)
C
      DO 22 I = 1, 100
        ECPT_SAVE(I) = ECPT(I)
   22 CONTINUE
C
C     Call KSOLID 6 times, once for each grid point as NPVT
C
      WRITE(6,*) 'KSOLID Test Driver - Wedge Element (ITYPE=1)'
      WRITE(6,*) '=============================================='
      WRITE(6,*) 'E = 200e9, NU = 0.3, G = ', GG
      WRITE(6,*) 'ITYPE = 1 (3 tetrahedra decomposition)'
      WRITE(6,*) ''
      WRITE(6,*) 'Nodes:'
      WRITE(6,*) '  1: (0,0,0)  2: (1,0,0)  3: (0,1,0)'
      WRITE(6,*) '  4: (0,0,1)  5: (1,0,1)  6: (0,1,1)'
      WRITE(6,*) ''
C
      DO 100 IPVT = 1, 6
C
C       Restore ECPT (KSOLID modifies it)
        DO 24 I = 1, 100
          ECPT(I) = ECPT_SAVE(I)
   24   CONTINUE
C
C       Set the pivot grid point
        NPVT = 100 + IPVT
C
        WRITE(6,*) '  Calling KSOLID with NPVT =', NPVT,
     1             ' (node', IPVT, ')'
C
        CALL KSOLID(1)
C
        IF (NOGOO) THEN
          WRITE(6,*) '  ERROR: NOGOO set after pivot', IPVT
          STOP
        END IF
  100 CONTINUE
C
C     Retrieve the assembled matrix from our SMA1B capture
C
      CALL GET_KMATRIX_WEDGE(KFULL)
C
C     Print summary information
C
      WRITE(6,*) ''
      WRITE(6,*) 'Full 18x18 Stiffness Matrix (diagonal entries):'
      WRITE(6,*) '================================================'
      DO 200 I = 1, 18
        WRITE(6,300) I, I, KFULL(I,I)
  200 CONTINUE
  300 FORMAT('  K(',I2,',',I2,') = ',E16.9)
C
C     Print row sums (rigid body mode check - should be ~0)
C
      WRITE(6,*) ''
      WRITE(6,*) 'Row sums (rigid body mode check):'
      DO 400 I = 1, 18
        ROWSUM = 0.0D0
        DO 350 J = 1, 18
          ROWSUM = ROWSUM + KFULL(I,J)
  350   CONTINUE
        WRITE(6,450) I, ROWSUM
  400 CONTINUE
  450 FORMAT('  Row ',I2,' sum = ',E16.9)
C
C     Print trace
C
      TRACE = 0.0D0
      DO 500 I = 1, 18
        TRACE = TRACE + KFULL(I,I)
  500 CONTINUE
      WRITE(6,*) ''
      WRITE(6,550) TRACE
  550 FORMAT('  Matrix trace = ',E16.9)
C
C     Print full matrix (first 9x9 block for readability)
C
      WRITE(6,*) ''
      WRITE(6,*) 'Top-left 9x9 block of stiffness matrix:'
      WRITE(6,*) '========================================='
      DO 600 I = 1, 9
        WRITE(6,650) (KFULL(I,J), J=1,9)
  600 CONTINUE
  650 FORMAT(9(1X,E11.4))
C
C     Symmetry check
C
      WRITE(6,*) ''
      WRITE(6,*) 'Symmetry check (max |K(i,j)-K(j,i)|):'
      ROWSUM = 0.0D0
      DO 700 I = 1, 18
        DO 690 J = I+1, 18
          IF (ABS(KFULL(I,J) - KFULL(J,I)) .GT. ROWSUM) THEN
            ROWSUM = ABS(KFULL(I,J) - KFULL(J,I))
          END IF
  690   CONTINUE
  700 CONTINUE
      WRITE(6,750) ROWSUM
  750 FORMAT('  Max asymmetry = ',E16.9)
C
C     Print the FULL 18x18 matrix, one row per line, between explicit
C     delimiters so the harness can extract it without ambiguity.
C     Added 2026-08-19: the stored reference previously covered only the
C     top-left 9x9 block, so generated code that correctly emitted the
C     full 18x18 was recorded as a dimension mismatch, i.e. as a numerical
C     failure. See 论文一_实验清单.md, "验证契约缺陷".
C
      WRITE(6,*) ''
      WRITE(6,*) 'BEGIN_FULL_18X18'
      DO 800 I = 1, 18
        WRITE(6,850) (KFULL(I,J), J=1,18)
  800 CONTINUE
  850 FORMAT(18(1X,E16.9))
      WRITE(6,*) 'END_FULL_18X18'
C
      IF (NOGOO) THEN
        WRITE(6,*) 'WARNING: NOGOO flag was set - errors occurred!'
      ELSE
        WRITE(6,*) ''
        WRITE(6,*) 'Computation completed successfully.'
      END IF
C
      STOP
      END
