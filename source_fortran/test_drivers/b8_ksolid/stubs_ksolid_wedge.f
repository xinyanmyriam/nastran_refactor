C     STUB SUBROUTINES FOR KSOLID WEDGE TEST DRIVER
C     =================================================
C     KSOLID(1) decomposes a wedge into 3 tetrahedra.
C     We need to capture outputs from KTETRA via SMA1B for 6 nodes.
C
C     ===================================================================
C     CAPTURE BUFFER: stores all 3x3 sub-blocks for 6x6 node pairs
C     ===================================================================
C
      BLOCK DATA CAPTUREDATAW
      IMPLICIT NONE
      DOUBLE PRECISION KSTORE(3,3,6,6)
      INTEGER          NCAPTURE
      COMMON /CAPTURW/ KSTORE, NCAPTURE
      DATA NCAPTURE /0/
      DATA KSTORE /324*0.0D0/
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
C     SMA1B - Stiffness matrix assembly stub for wedge (6 nodes)
C     Captures sub-blocks and accumulates into global 18x18 matrix
C     ===================================================================
C
      SUBROUTINE SMA1B (KIJ, IGRID, ICODE, IFFILE, FACTOR)
      IMPLICIT NONE
      DOUBLE PRECISION KIJ(6,6)
      INTEGER          IGRID, ICODE, IFFILE
      DOUBLE PRECISION FACTOR
C
      DOUBLE PRECISION KSTORE(3,3,6,6)
      INTEGER          NCAPTURE
      COMMON /CAPTURW/ KSTORE, NCAPTURE
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
      INTEGER GRIDS(6)
C
C     Only capture the primary stiffness (FACTOR = 0)
C
      IF (FACTOR .NE. 0.0D0) RETURN
C
C     The original wedge grid SIL IDs
C
      GRIDS(1) = 101
      GRIDS(2) = 102
      GRIDS(3) = 103
      GRIDS(4) = 104
      GRIDS(5) = 105
      GRIDS(6) = 106
C
C     Determine which pivot point (row block) this is
C
      DO 10 IPVT = 1, 6
        IF (NPVT .EQ. GRIDS(IPVT)) GO TO 20
   10 CONTINUE
      RETURN
C
C     Determine which grid point (column block) this is
C
   20 DO 30 IGRP = 1, 6
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
C     GET_KMATRIX_WEDGE - Retrieve the assembled 18x18 stiffness matrix
C     ===================================================================
C
      SUBROUTINE GET_KMATRIX_WEDGE (KFULL)
      IMPLICIT NONE
      DOUBLE PRECISION KFULL(18,18)
C
      DOUBLE PRECISION KSTORE(3,3,6,6)
      INTEGER          NCAPTURE
      COMMON /CAPTURW/ KSTORE, NCAPTURE
C
      INTEGER I, J, IPVT, IGRP, IROW, ICOL
C
      DO 100 IPVT = 1, 6
        DO 90 IGRP = 1, 6
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
C     KPLTST - Planarity test stub
C     For our test geometry all faces are planar, so just return
C     ===================================================================
C
      SUBROUTINE KPLTST (P1, P2, P3, P4)
      IMPLICIT NONE
      REAL P1(3), P2(3), P3(3), P4(3)
C     All faces of the test wedge are planar - do nothing
      RETURN
      END
C
C     ===================================================================
C     TRANSD - Coordinate transformation stub (all CSID=0)
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
