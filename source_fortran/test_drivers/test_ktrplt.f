C     TEST DRIVER FOR KTRPLT - TRIANGULAR PLATE BENDING ELEMENT
C     ===========================================================
C     This program calls KTRPLT 3 times (once per grid point as NPVT)
C     to assemble the full 18x18 stiffness matrix (6 DOF per node in 
C     SMA1B format, of which only 3 are used for bending: w, thx, thy).
C
C     Test geometry (equilateral triangle, side=1.0):
C       Node 1: (0, 0, 0)
C       Node 2: (1, 0, 0)
C       Node 3: (0.5, 0.866025, 0)
C
C     Material: E=200e9, NU=0.3
C     Plate properties: thickness t=0.01, I = t^3/12
C
C     DOF layout in SMA1B 6x6 blocks for plate bending:
C       Position 1,2 = zero (membrane u, v)
C       Position 3 = w (transverse displacement)
C       Position 4 = theta_x (rotation about x)
C       Position 5 = theta_y (rotation about y)
C       Position 6 = zero (theta_z)
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
C     SMA1DP - double precision work area
      DOUBLE PRECISION DPWORK(400)
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
      DOUBLE PRECISION KFULL(18,18)
      INTEGER          IPVT, I, J, II, JJ
      DOUBLE PRECISION K9(9,9)
      REAL             SQRT3O2
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
      SQRT3O2 = 0.8660254
C
C     ECPT layout for TRPLT:
      NECPT(1) = 1        ! Element ID
      NECPT(2) = 101      ! Grid Point A
      NECPT(3) = 102      ! Grid Point B
      NECPT(4) = 103      ! Grid Point C
      ECPT(5) = 0.0       ! THETA (material angle, degrees)
      NECPT(6) = 1        ! Material ID 1
      ECPT(7) = 8.3333333E-8  ! I = t^3/12, t=0.01
      NECPT(8) = 0        ! Material ID 2 (0=none)
      ECPT(9) = 0.0       ! T2 (transverse shear, 0=none)
      ECPT(10) = 0.0      ! Non-structural mass
      ECPT(11) = -0.005   ! Z1
      ECPT(12) = 0.005    ! Z2
C     Node 1: (0, 0, 0), CSID=0
      NECPT(13) = 0
      ECPT(14) = 0.0
      ECPT(15) = 0.0
      ECPT(16) = 0.0
C     Node 2: (1, 0, 0), CSID=0
      NECPT(17) = 0
      ECPT(18) = 1.0
      ECPT(19) = 0.0
      ECPT(20) = 0.0
C     Node 3: (0.5, sqrt(3)/2, 0), CSID=0
      NECPT(21) = 0
      ECPT(22) = 0.5
      ECPT(23) = SQRT3O2
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
      WRITE(6,*) 'KTRPLT Test Driver - Equilateral Triangle Plate'
      WRITE(6,*) '================================================'
      WRITE(6,*) 'E = 200e9, NU = 0.3, t = 0.01'
      WRITE(6,*) 'I = t^3/12 = ', ECPT(7)
      WRITE(6,*) 'Node 1: (0, 0, 0)'
      WRITE(6,*) 'Node 2: (1, 0, 0)'
      WRITE(6,*) 'Node 3: (0.5, 0.866025, 0)'
      WRITE(6,*) ''
C
C     Call KTRPLT 3 times, once for each grid point as pivot
      DO 100 IPVT = 1, 3
        NPVT = NECPT(IPVT + 1)
        CURPVT = IPVT
C
C       Zero the workspace
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
C     Retrieve the assembled matrix
      CALL GET_KMATRIX_PLATE(KFULL, 18)
C
C     Extract the 9x9 plate bending matrix
C     Bending DOFs are at positions 3,4,5 of each 6-DOF block
C     (1-indexed: pos 3=w, pos 4=theta_x, pos 5=theta_y)
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
      WRITE(6,*) 'DOF: (w1,thx1,thy1, w2,thx2,thy2, w3,thx3,thy3)'
      WRITE(6,*) '================================================'
      DO 300 I = 1, 9
        WRITE(6,310) (K9(I,J), J=1,9)
  300 CONTINUE
  310 FORMAT(9(1X,E14.7))
C
C     Print by 3x3 blocks
      WRITE(6,*) ''
      WRITE(6,*) '3x3 Sub-blocks:'
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
      WRITE(6,*) ''
      IF (NOGO .NE. 0) THEN
        WRITE(6,*) 'WARNING: NOGO flag was set!'
      ELSE
        WRITE(6,*) 'Computation completed successfully.'
      END IF
C
      STOP
      END
