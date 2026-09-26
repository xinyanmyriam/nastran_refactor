C     STUB SUBROUTINES FOR KELBOW TEST DRIVER
C     ==========================================
C
C     These stubs replace NASTRAN infrastructure subroutines that are
C     not needed for standalone testing of the elbow element stiffness.
C
C     ===================================================================
C     CAPTURE BUFFER: stores the 6x6 sub-blocks output by KELBOW
C     via SMA1B calls. KELBOW outputs 2 blocks per pivot call:
C       - K(NPVT,NPVT): self-coupling (6x6)
C       - K(NPVT,J): coupling to other grid (6x6)
C     Called twice (pivot A, pivot B) = 4 total 6x6 blocks.
C     KSTORE(:,:, pivot_idx, block_idx)
C       pivot_idx=1 for A, 2 for B
C       block_idx=1 for self, 2 for coupling
C     ===================================================================
C
      BLOCK DATA CAPTUREDATA
      IMPLICIT NONE
      DOUBLE PRECISION KSTORE(6,6,2,2)
      INTEGER          NCALLS, CPIVOT
      COMMON /CAPTURE/ KSTORE, NCALLS, CPIVOT
      DATA NCALLS /0/
      DATA CPIVOT /0/
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
      INTEGER   MATIDC, MATFLG
      REAL      ELTEMP, STRESS, SINTH, COSTH
      COMMON /MATIN / MATIDC, MATFLG, ELTEMP, STRESS, SINTH, COSTH
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
      SIGT = 250.0E6
      SIGC = 250.0E6
      SIGS = 145.0E6
C
      RETURN
      END
C
C     ===================================================================
C     SMA1B - Stiffness matrix assembly stub
C     Captures 6x6 sub-blocks from KELBOW.
C     KELBOW calls SMA1B twice per pivot:
C       1st call: K(NPVT,NPVT) self-coupling
C       2nd call: K(NPVT,J) coupling
C     INDEX = ISILNO(pivot), ICODE = -1
C     ===================================================================
C
      SUBROUTINE SMA1B (KIJ, INDEX, ICODE, IFFILE, FACTOR)
      IMPLICIT NONE
      DOUBLE PRECISION KIJ(6,6)
      INTEGER          INDEX, ICODE, IFFILE
      DOUBLE PRECISION FACTOR
C
      DOUBLE PRECISION KSTORE(6,6,2,2)
      INTEGER          NCALLS, CPIVOT
      COMMON /CAPTURE/ KSTORE, NCALLS, CPIVOT
C
      REAL             ECPT(100)
      INTEGER          NECPT(100)
      EQUIVALENCE     (NECPT(1), ECPT(1))
      COMMON /SMA1ET/ ECPT
C
      INTEGER          IOPT4, K4GGSW, NPVT, ISKP(11)
      COMMON /SMA1CL/ IOPT4, K4GGSW, NPVT, ISKP
C
      INTEGER I, J, IPVT, IBLK
C
C     Only capture the primary stiffness (FACTOR = 0)
C     Skip K4GG contributions (damping)
C
      IF (FACTOR .NE. 0.0D0) RETURN
C
C     Determine which pivot (1=A, 2=B)
      IPVT = 0
      IF (NPVT .EQ. NECPT(2)) IPVT = 1
      IF (NPVT .EQ. NECPT(3)) IPVT = 2
      IF (IPVT .EQ. 0) RETURN
C
C     Track calls per pivot. KELBOW calls SMA1B twice per pivot:
C       1st call = self-coupling K(NPVT,NPVT)
C       2nd call = coupling K(NPVT,J)
C     We detect new pivot by comparing to CPIVOT
C
      IF (IPVT .NE. CPIVOT) THEN
        CPIVOT = IPVT
        NCALLS = 0
      END IF
      NCALLS = NCALLS + 1
      IBLK = NCALLS
      IF (IBLK .GT. 2) IBLK = 2
C
C     Store the 6x6 sub-block
      DO 50 I = 1, 6
        DO 50 J = 1, 6
          KSTORE(I, J, IPVT, IBLK) = KIJ(I, J)
   50 CONTINUE
C
      RETURN
      END
C
C     ===================================================================
C     GET_KMATRIX_ELBOW - Retrieve assembled 12x12 stiffness matrix
C     Layout:
C       KSTORE(:,:,1,1) = K(A,A) -> rows 1-6, cols 1-6
C       KSTORE(:,:,1,2) = K(A,B) -> rows 1-6, cols 7-12
C       KSTORE(:,:,2,1) = K(B,B) -> rows 7-12, cols 7-12
C       KSTORE(:,:,2,2) = K(B,A) -> rows 7-12, cols 1-6
C     ===================================================================
C
      SUBROUTINE GET_KMATRIX_ELBOW (KFULL)
      IMPLICIT NONE
      DOUBLE PRECISION KFULL(12,12)
C
      DOUBLE PRECISION KSTORE(6,6,2,2)
      INTEGER          NCALLS, CPIVOT
      COMMON /CAPTURE/ KSTORE, NCALLS, CPIVOT
C
      INTEGER I, J
C
C     K(A,A): pivot=A(1), self(1) -> rows 1-6, cols 1-6
      DO 100 I = 1, 6
        DO 100 J = 1, 6
          KFULL(I, J) = KSTORE(I, J, 1, 1)
  100 CONTINUE
C
C     K(A,B): pivot=A(1), coupling(2) -> rows 1-6, cols 7-12
      DO 200 I = 1, 6
        DO 200 J = 1, 6
          KFULL(I, J+6) = KSTORE(I, J, 1, 2)
  200 CONTINUE
C
C     K(B,B): pivot=B(2), self(1) -> rows 7-12, cols 7-12
      DO 300 I = 1, 6
        DO 300 J = 1, 6
          KFULL(I+6, J+6) = KSTORE(I, J, 2, 1)
  300 CONTINUE
C
C     K(B,A): pivot=B(2), coupling(2) -> rows 7-12, cols 1-6
      DO 400 I = 1, 6
        DO 400 J = 1, 6
          KFULL(I+6, J) = KSTORE(I, J, 2, 2)
  400 CONTINUE
C
      RETURN
      END
C
C     ===================================================================
C     TRANSD - Coordinate transformation stub
C     Since all CSID = 0 (basic coords), this should not be called.
C     Returns identity if called anyway.
C     ===================================================================
C
      SUBROUTINE TRANSD (ICSID, T)
      IMPLICIT NONE
      INTEGER ICSID
      DOUBLE PRECISION T(9)
C
      INTEGER SYSBUF, NOUT
      COMMON /SYSTEM/ SYSBUF, NOUT
C
      WRITE(NOUT,*) 'NOTE: TRANSD called with CSID=', ICSID
      WRITE(NOUT,*) '      Returning identity (basic coords assumed)'
C     Return identity matrix (row-major: T11,T12,T13,T21,T22,T23,...)
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
      INTEGER ICODE, IMSG, IPARM
C
      INTEGER SYSBUF, NOUT
      COMMON /SYSTEM/ SYSBUF, NOUT
C
      WRITE(NOUT,*) 'MESAGE called: ICODE=', ICODE, ' IMSG=', IMSG
      WRITE(NOUT,*) '  IPARM=', IPARM
      IF (ICODE .LT. 0) THEN
        WRITE(NOUT,*) 'FATAL ERROR in MESAGE - stopping.'
        STOP
      END IF
      RETURN
      END
C
C     ===================================================================
C     HMAT - Heat material stub (should not be called since HEAT=.FALSE.)
C     ===================================================================
C
      SUBROUTINE HMAT (IELID)
      IMPLICIT NONE
      INTEGER IELID
C
      REAL FK
      COMMON /HMTOUT/ FK
C
      INTEGER SYSBUF, NOUT
      COMMON /SYSTEM/ SYSBUF, NOUT
C
      WRITE(NOUT,*) 'WARNING: HMAT called but HEAT should be .FALSE.'
      FK = 50.0
      RETURN
      END
