C     STUB SUBROUTINES FOR KTETRA TEST DRIVER
C     =========================================
C
C     These stubs replace NASTRAN infrastructure subroutines that are
C     not needed for standalone testing of the element stiffness matrix.
C
C     ===================================================================
C     CAPTURE BUFFER: stores all 3x3 sub-blocks output by KTETRA
C     via SMA1B calls. We store up to 16 blocks (4 pivots x 4 grids).
C     ===================================================================
C
      BLOCK DATA CAPTUREDATA
      IMPLICIT NONE
      DOUBLE PRECISION KSTORE(3,3,4,4)
      INTEGER          NBLOCKS
      COMMON /CAPTURE/ KSTORE, NBLOCKS
      DATA NBLOCKS /0/
      DATA KSTORE /144*0.0D0/
      END
C
C     ===================================================================
C     MAT - Material property routine stub
C     Sets material properties in COMMON /MATOUT/
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
C     Set material properties for steel-like material
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
C     Captures the 3x3 (in 6x6 format) KIJ sub-block for each grid point
C     ===================================================================
C
      SUBROUTINE SMA1B (KIJ, IGRID, ICODE, IFFILE, FACTOR)
      IMPLICIT NONE
      DOUBLE PRECISION KIJ(6,6)
      INTEGER          IGRID, ICODE, IFFILE
      DOUBLE PRECISION FACTOR
C
      DOUBLE PRECISION KSTORE(3,3,4,4)
      INTEGER          NBLOCKS
      COMMON /CAPTURE/ KSTORE, NBLOCKS
C
      REAL             ECPT(100)
      INTEGER          NECPT(100)
      EQUIVALENCE     (NECPT(1), ECPT(1))
      COMMON /SMA1ET/ ECPT
C
      INTEGER          IOPT4, K4GGSW, NPVT, ISKP(17)
      LOGICAL          NOGOO
      COMMON /SMA1CL/ IOPT4, K4GGSW, NPVT, ISKP, NOGOO
C
      INTEGER I, J, IPVT, IGRP
C
C     Only capture the primary stiffness (FACTOR = 0)
C     Skip K4GG contributions
C
      IF (FACTOR .NE. 0.0D0) RETURN
C
C     Determine which pivot point (row block) this is
C
      DO 10 IPVT = 1, 4
        IF (NPVT .EQ. NECPT(IPVT+2)) GO TO 20
   10 CONTINUE
      RETURN
C
C     Determine which grid point (column block) this is
C
   20 DO 30 IGRP = 1, 4
        IF (IGRID .EQ. NECPT(IGRP+2)) GO TO 40
   30 CONTINUE
      RETURN
C
C     Store the 3x3 sub-block (KIJ stores 3x3 in upper-left of 6x6)
C
   40 DO 50 I = 1, 3
        DO 50 J = 1, 3
          KSTORE(I, J, IPVT, IGRP) = KIJ(I, J)
   50 CONTINUE
C
      RETURN
      END
C
C     ===================================================================
C     GET_KMATRIX - Retrieve the assembled 12x12 stiffness matrix
C     ===================================================================
C
      SUBROUTINE GET_KMATRIX (KFULL)
      IMPLICIT NONE
      DOUBLE PRECISION KFULL(12,12)
C
      DOUBLE PRECISION KSTORE(3,3,4,4)
      INTEGER          NBLOCKS
      COMMON /CAPTURE/ KSTORE, NBLOCKS
C
      INTEGER I, J, IPVT, IGRP, IROW, ICOL
C
      DO 100 IPVT = 1, 4
        DO 90 IGRP = 1, 4
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
C     TRANSD - Coordinate transformation stub (should not be called
C     since all CSID = 0)
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
      WRITE(OUT,*) 'ERROR: TRANSD called with CSID=', ICSID
      WRITE(OUT,*) '       All CSID should be 0 for this test!'
C     Return identity matrix anyway
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
      WRITE(OUT,*) '  IPARM=', IPARM(1), IPARM(2)
      IF (ICODE .LT. 0) THEN
        WRITE(OUT,*) 'FATAL ERROR - stopping.'
        STOP
      END IF
      RETURN
      END
C
C     ===================================================================
C     WRITE - File write stub (for volume/surface output, not needed)
C     ===================================================================
C
      SUBROUTINE WRITE (IFILE, BUF, NWORDS, IFLG)
      IMPLICIT NONE
      INTEGER IFILE, NWORDS, IFLG
      REAL    BUF(1)
C     Do nothing - volume/surface output not needed for this test
      RETURN
      END
C
C     ===================================================================
C     HMAT - Heat material stub (should not be called since HEAT=.FALSE.)
C     ===================================================================
C
      SUBROUTINE HMAT (NECPT)
      IMPLICIT NONE
      INTEGER NECPT(4)
C
      REAL MATBUF(7)
      COMMON /HMTOUT/ MATBUF
C
      INTEGER OUT, SYSBUF
      LOGICAL NOGO
      COMMON /SYSTEM/ SYSBUF, OUT, NOGO
C
      WRITE(OUT,*) 'WARNING: HMAT called but HEAT should be .FALSE.'
C     Fill with thermal conductivity anyway
      MATBUF(1) = 50.0
      MATBUF(2) = 0.0
      MATBUF(3) = 0.0
      MATBUF(4) = 50.0
      MATBUF(5) = 0.0
      MATBUF(6) = 50.0
      MATBUF(7) = 0.0
      RETURN
      END
