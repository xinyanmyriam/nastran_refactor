C     STUB SUBROUTINES FOR KSOLID TEST DRIVER
C     =========================================
C     KSOLID calls KTETRA which calls SMA1B. We need to capture
C     the outputs from all 5 tetrahedra across 8 grid points.
C
C     ===================================================================
C     CAPTURE BUFFER: stores all 3x3 sub-blocks output by KTETRA/KSOLID
C     via SMA1B calls. For hexahedron: up to 8 nodes, 8 columns each.
C     ===================================================================
C
      BLOCK DATA CAPTUREDATA
      IMPLICIT NONE
      DOUBLE PRECISION KSTORE(3,3,8,8)
      INTEGER          NCAPTURE
      COMMON /CAPTURE/ KSTORE, NCAPTURE
      DATA NCAPTURE /0/
      DATA KSTORE /576*0.0D0/
      END
C
C     ===================================================================
C     MAT - Material property routine stub
C     ===================================================================
C
      SUBROUTINE MAT (IELID)
      IMPLICIT NONE
      INTEGER IELID
C
      INTEGER   MATID, INFLAG
      REAL      ELTEMP
      COMMON /MATIN / MATID, INFLAG, ELTEMP
C
      REAL      E, GG, NU, RHO, ALPHA, TSUB0, GSUBE, SIGT, SIGC, SIGS
      COMMON /MATOUT/ E, GG, NU, RHO, ALPHA, TSUB0, GSUBE,
     1                SIGT, SIGC, SIGS
C
      E = 200.0E9
      NU = 0.3
      GG = E / (2.0*(1.0+NU))
      RHO = 7800.0
      ALPHA = 12.0E-6
      TSUB0 = 20.0
      GSUBE = 0.0
      SIGT = 0.0
      SIGC = 0.0
      SIGS = 0.0
C
      RETURN
      END
C
C     ===================================================================
C     SMA1B - Stiffness matrix assembly stub
C     Captures sub-blocks and accumulates into global 24x24 matrix
C     ===================================================================
C
      SUBROUTINE SMA1B (KIJ, IGRID, ICODE, IFFILE, FACTOR)
      IMPLICIT NONE
      DOUBLE PRECISION KIJ(6,6)
      INTEGER          IGRID, ICODE, IFFILE
      DOUBLE PRECISION FACTOR
C
      DOUBLE PRECISION KSTORE(3,3,8,8)
      INTEGER          NCAPTURE
      COMMON /CAPTURE/ KSTORE, NCAPTURE
C
      REAL             ECPT(100)
      INTEGER          NECPT(100)
      EQUIVALENCE     (NECPT(1), ECPT(1))
      COMMON /SMA1ET/ ECPT
C
      INTEGER          IOPT4, K4GGSW, NPVT, ISKP(19)
      LOGICAL          NOGOO
      COMMON /SMA1CL/ IOPT4, K4GGSW, NPVT, ISKP, NOGOO
C
      INTEGER I, J, IPVT, IGRP
      INTEGER GRIDS(8)
C
C     Only capture the primary stiffness (FACTOR = 0)
C
      IF (FACTOR .NE. 0.0D0) RETURN
C
C     The grid SIL IDs for the hexahedron are stored in ECPT(53-60)
C     (KSOLID copies ECPT to ECPT+50 and then fills sub-tetra ECPT)
C     But NPVT and IGRID reference the original hex grid IDs (101-108)
C
      GRIDS(1) = 101
      GRIDS(2) = 102
      GRIDS(3) = 103
      GRIDS(4) = 104
      GRIDS(5) = 105
      GRIDS(6) = 106
      GRIDS(7) = 107
      GRIDS(8) = 108
C
C     Determine which pivot point (row block) this is
C
      DO 10 IPVT = 1, 8
        IF (NPVT .EQ. GRIDS(IPVT)) GO TO 20
   10 CONTINUE
      RETURN
C
C     Determine which grid point (column block) this is
C
   20 DO 30 IGRP = 1, 8
        IF (IGRID .EQ. GRIDS(IGRP)) GO TO 40
   30 CONTINUE
      RETURN
C
C     Accumulate the 3x3 sub-block (add since multiple tetras contribute)
C
   40 NCAPTURE = NCAPTURE + 1
      DO 50 I = 1, 3
        DO 50 J = 1, 3
          KSTORE(I, J, IPVT, IGRP) = KSTORE(I, J, IPVT, IGRP)
     1                                + KIJ(I, J)
   50 CONTINUE
C
      RETURN
      END
C
C     ===================================================================
C     GET_KMATRIX_HEX - Retrieve the assembled 24x24 stiffness matrix
C     ===================================================================
C
      SUBROUTINE GET_KMATRIX_HEX (KFULL)
      IMPLICIT NONE
      DOUBLE PRECISION KFULL(24,24)
C
      DOUBLE PRECISION KSTORE(3,3,8,8)
      INTEGER          NCAPTURE
      COMMON /CAPTURE/ KSTORE, NCAPTURE
C
      INTEGER I, J, IPVT, IGRP, IROW, ICOL
C
      DO 100 IPVT = 1, 8
        DO 90 IGRP = 1, 8
          DO 80 I = 1, 3
            DO 80 J = 1, 3
              IROW = (IPVT-1)*3 + I
              ICOL = (IGRP-1)*3 + J
              KFULL(IROW, ICOL) = KSTORE(I, J, IPVT, IGRP)
   80       CONTINUE
   90     CONTINUE
  100 CONTINUE
C
      RETURN
      END
C
C     ===================================================================
C     SAXB - Cross product A x B = C
C     ===================================================================
C
      SUBROUTINE SAXB (A, B, C)
      IMPLICIT NONE
      REAL A(3), B(3), C(3)
C
      C(1) = A(2)*B(3) - A(3)*B(2)
      C(2) = A(3)*B(1) - A(1)*B(3)
      C(3) = A(1)*B(2) - A(2)*B(1)
C
      RETURN
      END
C
C     ===================================================================
C     SADOTB - Dot product function A . B
C     ===================================================================
C
      FUNCTION SADOTB (A, B)
      IMPLICIT NONE
      REAL SADOTB
      REAL A(3), B(3)
C
      SADOTB = A(1)*B(1) + A(2)*B(2) + A(3)*B(3)
C
      RETURN
      END
C
C     ===================================================================
C     KPLTST - Planarity test stub (checks if 4 points are coplanar)
C     For a unit cube all faces are planar, so just return
C     ===================================================================
C
      SUBROUTINE KPLTST (P1, P2, P3, P4)
      IMPLICIT NONE
      REAL P1(3), P2(3), P3(3), P4(3)
C     All faces of a cube are planar - do nothing
      RETURN
      END
C
C     ===================================================================
C     TRANSD - Coordinate transformation stub
C     ===================================================================
C
      SUBROUTINE TRANSD (ICSID, T)
      IMPLICIT NONE
      INTEGER ICSID
      DOUBLE PRECISION T(9)
C
      INTEGER OUT, SYSBUF
      LOGICAL NOGO
      COMMON /SYSTEM/ SYSBUF, OUT, NOGO
C
      WRITE(OUT,*) 'WARNING: TRANSD called with CSID=', ICSID
C     Return identity matrix
      T(1) = 1.0D0
      T(2) = 0.0D0
      T(3) = 0.0D0
      T(4) = 0.0D0
      T(5) = 1.0D0
      T(6) = 0.0D0
      T(7) = 0.0D0
      T(8) = 0.0D0
      T(9) = 1.0D0
      RETURN
      END
C
C     ===================================================================
C     MESAGE - Error message stub
C     ===================================================================
C
      SUBROUTINE MESAGE (ICODE, IMSG, IPARM)
      IMPLICIT NONE
      INTEGER ICODE, IMSG, IPARM(2)
C
      INTEGER OUT, SYSBUF
      LOGICAL NOGO
      COMMON /SYSTEM/ SYSBUF, OUT, NOGO
C
      WRITE(OUT,*) 'MESAGE called: ICODE=', ICODE, ' IMSG=', IMSG
      IF (ICODE .LT. 0) THEN
        WRITE(OUT,*) 'FATAL ERROR - stopping.'
        STOP
      END IF
      RETURN
      END
C
C     ===================================================================
C     WRITE - File write stub (not needed for testing)
C     ===================================================================
C
      SUBROUTINE WRITE (IFILE, BUF, NWORDS, IFLG)
      IMPLICIT NONE
      INTEGER IFILE, NWORDS, IFLG
      REAL    BUF(1)
      RETURN
      END
C
C     ===================================================================
C     HMAT - Heat material stub
C     ===================================================================
C
      SUBROUTINE HMAT (NECPT)
      IMPLICIT NONE
      INTEGER NECPT(4)
      REAL MATBUF(7)
      COMMON /HMTOUT/ MATBUF
      MATBUF(1) = 50.0
      RETURN
      END
