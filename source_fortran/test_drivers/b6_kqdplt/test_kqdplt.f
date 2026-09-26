C     TEST DRIVER FOR KQDPLT - QUADRILATERAL PLATE BENDING ELEMENT
C     ==============================================================
C     This program calls KQDPLT 4 times (once per grid point as NPVT)
C     to assemble the full 24x24 stiffness matrix (6 DOF per node in
C     SMA1B format, of which only 3 are used for bending: w, thx, thy).
C
C     Test geometry (unit square):
C       Node 1: (0, 0, 0)
C       Node 2: (1, 0, 0)
C       Node 3: (1, 1, 0)
C       Node 4: (0, 1, 0)
C
C     Material: E=200e9, NU=0.3
C     Plate properties: thickness t=0.01, I = t^3/12
C
      PROGRAM TESTKQDPLT
      IMPLICIT NONE
C
C     COMMON blocks matching KQDPLT/KTRBSC expectations
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
C     SMA1DP - double precision work area (large for KQDPLT)
      DOUBLE PRECISION DPWORK(1000)
      COMMON /SMA1DP/ DPWORK
C
C     MACHIN - machine type
      INTEGER          MACH
      COMMON /MACHIN/ MACH
C
C     Capture buffer for SMA1B
      DOUBLE PRECISION KSTORE(6,6,4,4)
      INTEGER          CURPVT, NPVTS
      COMMON /CAPTUREP/ KSTORE, CURPVT, NPVTS
C
C     Local variables
      DOUBLE PRECISION KFULL(24,24)
      INTEGER          IPVT, I, J, II, JJ
      DOUBLE PRECISION K12(12,12)
      DOUBLE PRECISION ROWSUM, TRACE, SYMMERR, MAXSYM
      REAL             ECPT_SAVE(100)
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
      NPVTS = 4
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
C     Fill ECPT for QDPLT (quadrilateral plate bending)
      DO 10 I = 1, 100
        ECPT(I) = 0.0
   10 CONTINUE
C
C     ECPT(1) = Element ID
      NECPT(1) = 1
C     ECPT(2) = Grid Point A (SIL ID)
      NECPT(2) = 101
C     ECPT(3) = Grid Point B
      NECPT(3) = 102
C     ECPT(4) = Grid Point C
      NECPT(4) = 103
C     ECPT(5) = Grid Point D
      NECPT(5) = 104
C     ECPT(6) = THETA (material angle, degrees)
      ECPT(6) = 0.0
C     ECPT(7) = Material ID 1
      NECPT(7) = 1
C     ECPT(8) = I (moment of inertia per unit width = t^3/12)
      ECPT(8) = 8.3333333E-8
C     ECPT(9) = Material ID 2 (for transverse shear, 0 = none)
      NECPT(9) = 0
C     ECPT(10) = T2 (transverse shear thickness, 0 = none)
      ECPT(10) = 0.0
C     ECPT(11) = Non-structural mass
      ECPT(11) = 0.0
C     ECPT(12) = Z1 (fiber distance)
      ECPT(12) = -0.005
C     ECPT(13) = Z2 (fiber distance)
      ECPT(13) = 0.005
C
C     Grid point coordinates (in basic system, CSID=0):
C     Node 1: (0, 0, 0)
      NECPT(14) = 0
      ECPT(15) = 0.0
      ECPT(16) = 0.0
      ECPT(17) = 0.0
C     Node 2: (1, 0, 0)
      NECPT(18) = 0
      ECPT(19) = 1.0
      ECPT(20) = 0.0
      ECPT(21) = 0.0
C     Node 3: (1, 1, 0)
      NECPT(22) = 0
      ECPT(23) = 1.0
      ECPT(24) = 1.0
      ECPT(25) = 0.0
C     Node 4: (0, 1, 0)
      NECPT(26) = 0
      ECPT(27) = 0.0
      ECPT(28) = 1.0
      ECPT(29) = 0.0
C     ECPT(30) = Element temperature
      ECPT(30) = 20.0
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
      DO 20 I = 1, 24
        DO 20 J = 1, 24
          KFULL(I, J) = 0.0D0
   20 CONTINUE
C
      WRITE(6,*) 'KQDPLT Test Driver - Unit Square Plate'
      WRITE(6,*) '======================================='
      WRITE(6,*) 'E = 200e9, NU = 0.3, t = 0.01'
      WRITE(6,*) 'I = t^3/12 = ', ECPT(8)
      WRITE(6,*) 'Node 1: (0, 0, 0)'
      WRITE(6,*) 'Node 2: (1, 0, 0)'
      WRITE(6,*) 'Node 3: (1, 1, 0)'
      WRITE(6,*) 'Node 4: (0, 1, 0)'
      WRITE(6,*) ''
C
C     Save ECPT (KQDPLT modifies it internally with shift operation)
      DO 25 I = 1, 100
        ECPT_SAVE(I) = ECPT(I)
   25 CONTINUE
C
C     Call KQDPLT 4 times, once for each grid point as pivot
      DO 100 IPVT = 1, 4
C
C       Restore ECPT before each call (KQDPLT shifts elements)
        DO 28 I = 1, 100
          ECPT(I) = ECPT_SAVE(I)
   28   CONTINUE
C
        NPVT = NECPT(IPVT + 1)
        CURPVT = IPVT
C
C       Zero the workspace
        DO 30 I = 1, 1000
          DPWORK(I) = 0.0D0
   30   CONTINUE
C
        WRITE(6,*) '  Calling KQDPLT with pivot =', IPVT,
     1             ' (SIL=', NPVT, ')'
C
        CALL KQDPLT
C
        IF (NOGO .NE. 0) THEN
          WRITE(6,*) '  ERROR: NOGO flag set after pivot', IPVT
          STOP
        END IF
C
  100 CONTINUE
C
C     Retrieve the assembled matrix
      CALL GET_KMATRIX_PLATE(KFULL, 24)
C
C     Extract the 12x12 plate bending matrix
C     From the 24x24 (6 DOF per node), plate bending DOFs are
C     positions 3,4,5 of each 6-DOF block (w, theta_x, theta_y)
      DO 200 I = 1, 4
        DO 200 J = 1, 4
          DO 200 II = 1, 3
            DO 200 JJ = 1, 3
              K12(3*(I-1)+II, 3*(J-1)+JJ) =
     1            KFULL(6*(I-1)+II+2, 6*(J-1)+JJ+2)
  200 CONTINUE
C
      WRITE(6,*) ''
      WRITE(6,*) '12x12 Plate Bending Stiffness Matrix:'
      WRITE(6,*) 'DOF: (w1,thx1,thy1, w2,..., w4,thx4,thy4)'
      WRITE(6,*) '================================================'
      DO 300 I = 1, 12
        WRITE(6,310) (K12(I,J), J=1,12)
  300 CONTINUE
  310 FORMAT(12(1X,E12.5))
C
C     Print by 3x3 blocks
      WRITE(6,*) ''
      WRITE(6,*) '3x3 Sub-blocks:'
      DO 400 I = 1, 4
        DO 400 J = 1, 4
          WRITE(6,410) I, J
          DO 350 II = 1, 3
            WRITE(6,420) (K12(3*(I-1)+II, 3*(J-1)+JJ), JJ=1,3)
  350     CONTINUE
  400 CONTINUE
  410 FORMAT('  K(',I1,',',I1,'):')
  420 FORMAT(3(2X,E16.9))
C
C     Symmetry check
      WRITE(6,*) ''
      WRITE(6,*) 'Symmetry check (max |K(i,j)-K(j,i)|):'
      MAXSYM = 0.0D0
      DO 500 I = 1, 12
        DO 500 J = I+1, 12
          SYMMERR = DABS(K12(I,J) - K12(J,I))
          IF (SYMMERR .GT. MAXSYM) MAXSYM = SYMMERR
  500 CONTINUE
      WRITE(6,510) MAXSYM
  510 FORMAT('  Max symmetry error = ', E16.9)
C
C     Row sum check (rigid body mode)
      WRITE(6,*) ''
      WRITE(6,*) 'Row sums (rigid body mode check):'
      DO 600 I = 1, 12
        ROWSUM = 0.0D0
        DO 550 J = 1, 12
          ROWSUM = ROWSUM + K12(I,J)
  550   CONTINUE
        WRITE(6,610) I, ROWSUM
  600 CONTINUE
  610 FORMAT('  Row ',I2,' sum = ',E16.9)
C
C     Trace
      TRACE = 0.0D0
      DO 700 I = 1, 12
        TRACE = TRACE + K12(I,I)
  700 CONTINUE
      WRITE(6,*) ''
      WRITE(6,710) TRACE
  710 FORMAT('  Matrix trace = ',E16.9)
C
      WRITE(6,*) ''
      IF (NOGO .NE. 0) THEN
        WRITE(6,*) 'WARNING: NOGO flag was set!'
      ELSE
        WRITE(6,*) 'Computation completed successfully.'
      END IF
C
      STOP
      END
