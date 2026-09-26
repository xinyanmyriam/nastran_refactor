C     TEST DRIVER FOR KTETRA - TETRAHEDRON ELEMENT STIFFNESS MATRIX
C     =============================================================
C     This program calls KTETRA 4 times (once per grid point as NPVT)
C     to assemble the full 12x12 stiffness matrix for a unit tetrahedron.
C
C     Test geometry (unit tetrahedron):
C       Node 1: (0, 0, 0)
C       Node 2: (1, 0, 0)
C       Node 3: (0, 1, 0)
C       Node 4: (0, 0, 1)
C
C     Material: E=200e9, NU=0.3, G=E/(2*(1+NU))
C
      PROGRAM TESTKTETRA
      IMPLICIT NONE
C
C     COMMON blocks matching KTETRA expectations
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
      INTEGER          IOPT4, K4GGSW, NPVT, ISKP(17)
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
      REAL             SKIP(16), VOLUME, SURFAC
      COMMON /BLANK / SKIP, VOLUME, SURFAC
C
      CHARACTER        UFM*23
      COMMON /XMSSG / UFM
C
C     Local variables for assembling the 12x12 matrix
C
      DOUBLE PRECISION KFULL(12,12)
      INTEGER          IPVT, I, J, IROW, ICOL
C
C     Initialize
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
C     Set material properties (will also be set by MAT stub)
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
C     Fill ECPT - Element Connection Property Table
C
C     Zero it out first
      DO 10 I = 1, 100
        ECPT(I) = 0.0
   10 CONTINUE
C
C     ECPT(1) = Element ID
      NECPT(1) = 1
C     ECPT(2) = Material ID
      NECPT(2) = 1
C     ECPT(3-6) = SIL grid point IDs
      NECPT(3) = 101
      NECPT(4) = 102
      NECPT(5) = 103
      NECPT(6) = 104
C
C     Grid 1: (0, 0, 0), CSID=0
      NECPT(7) = 0
      ECPT(8)  = 0.0
      ECPT(9)  = 0.0
      ECPT(10) = 0.0
C     Grid 2: (1, 0, 0), CSID=0
      NECPT(11) = 0
      ECPT(12)  = 1.0
      ECPT(13)  = 0.0
      ECPT(14)  = 0.0
C     Grid 3: (0, 1, 0), CSID=0
      NECPT(15) = 0
      ECPT(16)  = 0.0
      ECPT(17)  = 1.0
      ECPT(18)  = 0.0
C     Grid 4: (0, 0, 1), CSID=0
      NECPT(19) = 0
      ECPT(20)  = 0.0
      ECPT(21)  = 0.0
      ECPT(22)  = 1.0
C     Element temperature
      ECPT(23) = 20.0
C
C     Zero the full stiffness matrix
C
      DO 20 I = 1, 12
        DO 15 J = 1, 12
          KFULL(I,J) = 0.0D0
   15   CONTINUE
   20 CONTINUE
C
C     Call KTETRA 4 times, once for each grid point as pivot
C     KTETRA produces a pivot row: for each of the 4 grids it outputs
C     a 3x3 sub-block KIJ via SMA1B. Our SMA1B stub captures these.
C
      WRITE(6,*) 'KTETRA Test Driver - Unit Tetrahedron'
      WRITE(6,*) '======================================'
      WRITE(6,*) 'E = 200e9, NU = 0.3, G = ', GG
      WRITE(6,*) ''
C
      DO 100 IPVT = 1, 4
        NPVT = NECPT(IPVT + 2)
C
C       Reset element ID flag to force H-matrix recalculation
C       (KTETRA uses IDFLAG internally; we change NECPT(1) to force it)
C
        NECPT(1) = IPVT
C
        CALL KTETRA(0, 0)
C
  100 CONTINUE
C
C     Now retrieve the assembled matrix from our SMA1B capture
C     (stored in a separate common block by our stub)
C
      CALL GET_KMATRIX(KFULL)
C
C     Print the 12x12 stiffness matrix
C
      WRITE(6,*) ''
      WRITE(6,*) 'Full 12x12 Stiffness Matrix'
      WRITE(6,*) 'Row/Col: G1x,G1y,G1z,G2x,G2y,G2z,...'
      WRITE(6,*) '=========================================='
      DO 200 I = 1, 12
        WRITE(6,300) (KFULL(I,J), J=1,12)
  200 CONTINUE
  300 FORMAT(12(1X,E14.7))
C
C     Also print in a more readable format
      WRITE(6,*) ''
      WRITE(6,*) 'Stiffness matrix (by 3x3 blocks):'
      WRITE(6,*) ''
      DO 400 IROW = 0, 3
        DO 350 ICOL = 0, 3
          WRITE(6,450) IROW+1, ICOL+1
          DO 340 I = 1, 3
            WRITE(6,460) (KFULL(3*IROW+I, 3*ICOL+J), J=1,3)
  340     CONTINUE
  350   CONTINUE
  400 CONTINUE
  450 FORMAT('  Block K(',I1,',',I1,'):')
  460 FORMAT(3(2X,E16.9))
C
      IF (NOGOO) THEN
        WRITE(6,*) 'WARNING: NOGOO flag was set - errors occurred!'
      ELSE
        WRITE(6,*) 'Computation completed successfully.'
      END IF
C
      STOP
      END
