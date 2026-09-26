C     TEST DRIVER FOR KSOLID - HEXAHEDRON SOLID ELEMENT (5 TETRAS)
C     =============================================================
C     This program calls KSOLID(2) which decomposes a hexahedron into
C     5 tetrahedra and calls KTETRA for each. The full 24x24 stiffness
C     matrix is assembled from the sub-block outputs.
C
C     Test geometry: unit cube (0,0,0) to (1,1,1)
C       Node 1: (0, 0, 0)   Node 5: (0, 0, 1)
C       Node 2: (1, 0, 0)   Node 6: (1, 0, 1)
C       Node 3: (1, 1, 0)   Node 7: (1, 1, 1)
C       Node 4: (0, 1, 0)   Node 8: (0, 1, 1)
C
C     Material: E=200e9, NU=0.3
C
      PROGRAM TESTKSOLID
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
C     Local variables for assembling the 24x24 matrix
C
      DOUBLE PRECISION KFULL(24,24)
      INTEGER          I, J, IROW, ICOL, IPVT
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
C     Fill ECPT for hexahedron (ITYPE=2, 8 nodes)
C     Layout from ksolid.f comments:
C     ECPT(1) = EL ID
C     ECPT(2) = MAT-ID
C     ECPT(3-10) = GRID-1 through GRID-8
C     ECPT(11) = CSID-1, ECPT(12-14) = X1,Y1,Z1
C     ECPT(15) = CSID-2, ECPT(16-18) = X2,Y2,Z2
C     ... and so on for all 8 nodes
C     ECPT(43) = EL-TEMP
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
C     Grid SIL IDs (nodes 1-8)
      NECPT(3)  = 101
      NECPT(4)  = 102
      NECPT(5)  = 103
      NECPT(6)  = 104
      NECPT(7)  = 105
      NECPT(8)  = 106
      NECPT(9)  = 107
      NECPT(10) = 108
C
C     Node 1: (0,0,0)
      NECPT(11) = 0
      ECPT(12) = 0.0
      ECPT(13) = 0.0
      ECPT(14) = 0.0
C
C     Node 2: (1,0,0)
      NECPT(15) = 0
      ECPT(16) = 1.0
      ECPT(17) = 0.0
      ECPT(18) = 0.0
C
C     Node 3: (1,1,0)
      NECPT(19) = 0
      ECPT(20) = 1.0
      ECPT(21) = 1.0
      ECPT(22) = 0.0
C
C     Node 4: (0,1,0)
      NECPT(23) = 0
      ECPT(24) = 0.0
      ECPT(25) = 1.0
      ECPT(26) = 0.0
C
C     Node 5: (0,0,1)
      NECPT(27) = 0
      ECPT(28) = 0.0
      ECPT(29) = 0.0
      ECPT(30) = 1.0
C
C     Node 6: (1,0,1)
      NECPT(31) = 0
      ECPT(32) = 1.0
      ECPT(33) = 0.0
      ECPT(34) = 1.0
C
C     Node 7: (1,1,1)
      NECPT(35) = 0
      ECPT(36) = 1.0
      ECPT(37) = 1.0
      ECPT(38) = 1.0
C
C     Node 8: (0,1,1)
      NECPT(39) = 0
      ECPT(40) = 0.0
      ECPT(41) = 1.0
      ECPT(42) = 1.0
C
C     Element temperature
      ECPT(43) = 20.0
C
C     Zero the full stiffness matrix
C
      DO 20 I = 1, 24
        DO 15 J = 1, 24
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
C     Call KSOLID 8 times, once for each grid point as NPVT
C     KSOLID decomposes the hex into tetrahedra, and for each tetra
C     calls KTETRA which processes only the current NPVT.
C
      WRITE(6,*) 'KSOLID Test Driver - Unit Cube Hexahedron'
      WRITE(6,*) '=========================================='
      WRITE(6,*) 'E = 200e9, NU = 0.3, G = ', GG
      WRITE(6,*) 'ITYPE = 2 (5 tetrahedra decomposition)'
      WRITE(6,*) ''
C
      DO 100 IPVT = 1, 8
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
        CALL KSOLID(2)
C
        IF (NOGOO) THEN
          WRITE(6,*) '  ERROR: NOGOO set after pivot', IPVT
          STOP
        END IF
  100 CONTINUE
C
C     Retrieve the assembled matrix from our SMA1B capture
C
      CALL GET_KMATRIX_HEX(KFULL)
C
C     Print summary information
C
      WRITE(6,*) ''
      WRITE(6,*) 'Full 24x24 Stiffness Matrix (diagonal entries):'
      WRITE(6,*) '================================================'
      DO 200 I = 1, 24
        WRITE(6,300) I, I, KFULL(I,I)
  200 CONTINUE
  300 FORMAT('  K(',I2,',',I2,') = ',E16.9)
C
C     Print row sums (should be ~0 for rigid body mode check)
C
      WRITE(6,*) ''
      WRITE(6,*) 'Row sums (rigid body mode check):'
      DO 400 I = 1, 24
        ROWSUM = 0.0D0
        DO 350 J = 1, 24
          ROWSUM = ROWSUM + KFULL(I,J)
  350   CONTINUE
        WRITE(6,450) I, ROWSUM
  400 CONTINUE
  450 FORMAT('  Row ',I2,' sum = ',E16.9)
C
C     Print trace
C
      TRACE = 0.0D0
      DO 500 I = 1, 24
        TRACE = TRACE + KFULL(I,I)
  500 CONTINUE
      WRITE(6,*) ''
      WRITE(6,550) TRACE
  550 FORMAT('  Matrix trace = ',E16.9)
C
C     Print full matrix (first 8x8 block for readability)
C
      WRITE(6,*) ''
      WRITE(6,*) 'Top-left 8x8 block of stiffness matrix:'
      WRITE(6,*) '========================================='
      DO 600 I = 1, 8
        WRITE(6,650) (KFULL(I,J), J=1,8)
  600 CONTINUE
  650 FORMAT(8(1X,E12.5))
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
