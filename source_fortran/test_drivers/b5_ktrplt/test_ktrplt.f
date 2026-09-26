C     TEST DRIVER FOR KTRPLT - TRIANGULAR PLATE BENDING ELEMENT (B5)
C     ================================================================
C     This program calls KTRPLT 3 times (once per grid point as NPVT)
C     to assemble the full 18x18 stiffness matrix (6 DOF per node in
C     SMA1B format), then extracts the 9x9 bending portion.
C
C     Test geometry (right-angled triangle in XY plane):
C       Node A: (0, 0, 0)
C       Node B: (1, 0, 0)
C       Node C: (0, 1, 0)
C
C     Material: E=200e9, NU=0.3
C     Plate properties: thickness t=0.01
C       I = t^3/12 = 8.333e-8
C       T2 = 0 (no transverse shear)
C       Z1 = -0.005, Z2 = 0.005
C
C     DOF layout in SMA1B 6x6 blocks:
C       Position 3 = w (transverse displacement)
C       Position 4 = theta_x (rotation about x)
C       Position 5 = theta_y (rotation about y)
C       Positions 1,2,6 = zero for pure bending
C
      PROGRAM TESTKTRPLT
      IMPLICIT NONE
C
C     COMMON blocks matching KTRPLT/KTRBSC expectations
C
      REAL             ECPT(100)
      INTEGER          NECPT(100)
      EQUIVALENCE     (NECPT(1), ECPT(1))
      COMMON /SMA1ET/ ECPT
C
C     CONDAS - constants (DEGRA = degrees to radians at position 4)
      REAL             CONSTS(5)
      COMMON /CONDAS/ CONSTS
C
C     MATIN/MATOUT - material interface
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
C     SMA1IO - I/O control
      REAL             DUM1(10)
      INTEGER          IFKGG
      REAL             DUM2(1)
      INTEGER          IF4GG
      REAL             DUM3(23)
      COMMON /SMA1IO/ DUM1, IFKGG, DUM2, IF4GG, DUM3
C
C     SMA1CL - control parameters
      INTEGER          IOPT4, K4GGSW, NPVT
      REAL             DUMCL(7)
      INTEGER          LINK(10), IDETCK
      REAL             DODET
      INTEGER          NOGO
      COMMON /SMA1CL/ IOPT4, K4GGSW, NPVT, DUMCL, LINK, IDETCK,
     1                DODET, NOGO
C
C     SMA1DP - double precision work area (must be large enough)
      DOUBLE PRECISION DPWORK(400)
      COMMON /SMA1DP/ DPWORK
C
C     MACHIN - machine type
      INTEGER          MACH
      COMMON /MACHIN/ MACH
C
C     Capture buffer for SMA1B output
      DOUBLE PRECISION KSTORE(6,6,4,4)
      INTEGER          CURPVT, NPVTS
      COMMON /CAPTUREP/ KSTORE, CURPVT, NPVTS
C
C     Local variables
      DOUBLE PRECISION KFULL(18,18)
      INTEGER          IPVT, I, J, II, JJ
      DOUBLE PRECISION K9(9,9)
      DOUBLE PRECISION ROWSUM, TRACE, SYMMERR, MAXSYM
C
C     Initialize machine and constants
      MACH = 1
      CONSTS(1) = 3.14159265358979
      CONSTS(2) = 0.0
      CONSTS(3) = 0.0
      CONSTS(4) = 3.14159265358979 / 180.0
      CONSTS(5) = 0.0
C
C     Initialize control variables
      IOPT4 = 0
      K4GGSW = 0
      NOGO = 0
      IDETCK = 0
      DODET = 0.0
      IFKGG = 0
      IF4GG = 0
      NPVTS = 3
C
      DO 5 I = 1, 10
        DUM1(I) = 0.0
    5 CONTINUE
      DUM2(1) = 0.0
      DO 6 I = 1, 23
        DUM3(I) = 0.0
    6 CONTINUE
      DO 7 I = 1, 7
        DUMCL(I) = 0.0
    7 CONTINUE
      DO 8 I = 1, 10
        LINK(I) = 0
    8 CONTINUE
C
C     Set material properties
      GSUBE = 0.0
C
C     Fill ECPT for TRPLT (triangular plate bending)
      DO 10 I = 1, 100
        ECPT(I) = 0.0
   10 CONTINUE
C
C     ECPT layout for TRPLT (right-angled triangle):
      NECPT(1) = 1          ! Element ID
      NECPT(2) = 101        ! Grid Point A (SIL)
      NECPT(3) = 102        ! Grid Point B (SIL)
      NECPT(4) = 103        ! Grid Point C (SIL)
      ECPT(5) = 0.0         ! THETA (material angle, degrees)
      NECPT(6) = 1          ! Material ID 1
      ECPT(7) = 8.3333333E-8  ! I = t^3/12, t=0.01
      NECPT(8) = 0          ! Material ID 2 (0=no shear)
      ECPT(9) = 0.0         ! T2 (transverse shear, 0=none)
      ECPT(10) = 0.0        ! Non-structural mass
      ECPT(11) = -0.005     ! Z1 = -t/2
      ECPT(12) = 0.005      ! Z2 = +t/2
C     Node A: (0, 0, 0), CSID=0
      NECPT(13) = 0
      ECPT(14) = 0.0
      ECPT(15) = 0.0
      ECPT(16) = 0.0
C     Node B: (1, 0, 0), CSID=0
      NECPT(17) = 0
      ECPT(18) = 1.0
      ECPT(19) = 0.0
      ECPT(20) = 0.0
C     Node C: (0, 1, 0), CSID=0
      NECPT(21) = 0
      ECPT(22) = 0.0
      ECPT(23) = 1.0
      ECPT(24) = 0.0
C     Element temperature
      ECPT(25) = 20.0
C
C     Zero capture buffer
      DO 15 I = 1, 4
        DO 15 J = 1, 4
          DO 15 II = 1, 6
            DO 15 JJ = 1, 6
              KSTORE(II, JJ, I, J) = 0.0D0
   15 CONTINUE
C
C     Zero full matrix
      DO 20 I = 1, 18
        DO 20 J = 1, 18
          KFULL(I, J) = 0.0D0
   20 CONTINUE
C
      WRITE(6,*) 'KTRPLT Test Driver - Right-Angled Triangle Plate'
      WRITE(6,*) '================================================='
      WRITE(6,*) 'E = 200e9, NU = 0.3, t = 0.01'
      WRITE(6,*) 'I = t^3/12 = ', ECPT(7)
      WRITE(6,*) 'Node A: (0, 0, 0)'
      WRITE(6,*) 'Node B: (1, 0, 0)'
      WRITE(6,*) 'Node C: (0, 1, 0)'
      WRITE(6,*) 'THETA = 0 degrees, T2 = 0 (no shear)'
      WRITE(6,*) ''
C
C     Call KTRPLT 3 times, once for each grid point as pivot
      DO 100 IPVT = 1, 3
        NPVT = NECPT(IPVT + 1)
        CURPVT = IPVT
C
C       Zero the workspace before each call
        DO 30 I = 1, 400
          DPWORK(I) = 0.0D0
   30   CONTINUE
C
        WRITE(6,*) '  Calling KTRPLT with pivot =', IPVT,
     1             ' (SIL=', NPVT, ')'
C
        CALL KTRPLT
C
        IF (NOGO .NE. 0) THEN
          WRITE(6,*) '  ERROR: NOGO flag set after pivot', IPVT
          STOP
        END IF
  100 CONTINUE
C
C     Retrieve the assembled 18x18 matrix from SMA1B capture
      CALL GET_KMATRIX_PLATE(KFULL, 18)
C
C     Extract the 9x9 plate bending matrix
C     Bending DOFs are at positions 3,4,5 of each 6-DOF block
C     (w, theta_x, theta_y)
      DO 200 I = 1, 3
        DO 200 J = 1, 3
          DO 200 II = 1, 3
            DO 200 JJ = 1, 3
              K9(3*(I-1)+II, 3*(J-1)+JJ) =
     1            KFULL(6*(I-1)+II+2, 6*(J-1)+JJ+2)
  200 CONTINUE
C
      WRITE(6,*) ''
      WRITE(6,*) '9x9 Plate Bending Stiffness Matrix:'
      WRITE(6,*) 'DOF: (wA,thxA,thyA, wB,thxB,thyB, wC,thxC,thyC)'
      WRITE(6,*) '================================================='
      DO 300 I = 1, 9
        WRITE(6,310) (K9(I,J), J=1,9)
  300 CONTINUE
  310 FORMAT(9(1X,E14.7))
C
C     Print by 3x3 blocks
      WRITE(6,*) ''
      WRITE(6,*) '3x3 Sub-blocks (Kij for node i, node j):'
      DO 400 I = 1, 3
        DO 400 J = 1, 3
          WRITE(6,410) I, J
          DO 350 II = 1, 3
            WRITE(6,420) (K9(3*(I-1)+II, 3*(J-1)+JJ), JJ=1,3)
  350     CONTINUE
  400 CONTINUE
  410 FORMAT('  K(',I1,',',I1,'):')
  420 FORMAT(3(2X,E16.9))
C
C     Print full 6x6 blocks as output by SMA1B
      WRITE(6,*) ''
      WRITE(6,*) '6x6 SMA1B Output Blocks (pivot, grid):'
      DO 450 I = 1, 3
        DO 450 J = 1, 3
          WRITE(6,455) I, J
          DO 445 II = 1, 6
            WRITE(6,460) (KSTORE(II, JJ, I, J), JJ=1,6)
  445     CONTINUE
  450 CONTINUE
  455 FORMAT('  SMA1B Block (pvt=',I1,', grd=',I1,'):')
  460 FORMAT(6(1X,E13.6))
C
C     Symmetry check
      WRITE(6,*) ''
      WRITE(6,*) 'Symmetry check (max |K(i,j)-K(j,i)|):'
      MAXSYM = 0.0D0
      DO 500 I = 1, 9
        DO 500 J = I+1, 9
          SYMMERR = DABS(K9(I,J) - K9(J,I))
          IF (SYMMERR .GT. MAXSYM) MAXSYM = SYMMERR
  500 CONTINUE
      WRITE(6,510) MAXSYM
  510 FORMAT('  Max symmetry error = ', E16.9)
C
C     Row sum check (rigid body mode: constant w)
      WRITE(6,*) ''
      WRITE(6,*) 'Row sums (should be ~0 for rigid body w mode):'
      DO 600 I = 1, 9
        ROWSUM = 0.0D0
        DO 550 J = 1, 9
          ROWSUM = ROWSUM + K9(I,J)
  550   CONTINUE
        WRITE(6,610) I, ROWSUM
  600 CONTINUE
  610 FORMAT('  Row ',I1,' sum = ',E16.9)
C
C     Trace (sum of diagonal)
      TRACE = 0.0D0
      DO 700 I = 1, 9
        TRACE = TRACE + K9(I,I)
  700 CONTINUE
      WRITE(6,*) ''
      WRITE(6,710) TRACE
  710 FORMAT('  Matrix trace = ',E16.9)
C
C     Diagonal values
      WRITE(6,*) ''
      WRITE(6,*) 'Diagonal entries of 9x9 K:'
      DO 750 I = 1, 9
        WRITE(6,760) I, I, K9(I,I)
  750 CONTINUE
  760 FORMAT('  K(',I1,',',I1,') = ',E16.9)
C
      WRITE(6,*) ''
      IF (NOGO .NE. 0) THEN
        WRITE(6,*) 'WARNING: NOGO flag was set - errors occurred!'
      ELSE
        WRITE(6,*) 'Computation completed successfully.'
      END IF
C
      STOP
      END
