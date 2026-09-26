C     STUB SUBROUTINES FOR SROD1 TEST DRIVER
C     =========================================
C     These stubs replace NASTRAN infrastructure subroutines.
C     SROD1 calls: MAT, TRANSS, GMMATS
C     For basic coordinates (CSID=0), TRANSS/GMMATS are not called.
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
      REAL      E, G, NU, RHO, ALPHA, TSUB0, GSUBE, SIGT, SIGC, SIGS
      COMMON /MATOUT/ E, G, NU, RHO, ALPHA, TSUB0,
     1                GSUBE, SIGT, SIGC, SIGS
C
C     Steel-like material: E=2.1e11, NU=0.3
C
      E = 2.1E11
      NU = 0.3
      G = E / (2.0*(1.0+NU))
      RHO = 7800.0
      ALPHA = 12.0E-6
      TSUB0 = 20.0
      GSUBE = 0.0
      SIGT = 250.0E6
      SIGC = 250.0E6
      SIGS = 150.0E6
C
      RETURN
      END
C
C     ===================================================================
C     TRANSS - Coordinate transformation (single precision version)
C     Should not be called since both nodes are in basic (CSID=0)
C     ===================================================================
C
      SUBROUTINE TRANSS (ICSID, TI)
      IMPLICIT NONE
      INTEGER ICSID
      REAL    TI(9)
C
      WRITE(6,*) 'WARNING: TRANSS called with CSID=', ICSID
      WRITE(6,*) '  All CSID should be 0 for this test!'
C     Return identity
      TI(1) = 1.0
      TI(2) = 0.0
      TI(3) = 0.0
      TI(4) = 0.0
      TI(5) = 1.0
      TI(6) = 0.0
      TI(7) = 0.0
      TI(8) = 0.0
      TI(9) = 1.0
      RETURN
      END
C
C     ===================================================================
C     GMMATS - General matrix multiply (single precision)
C     C = A * B  where dimensions are specified
C     GMMATS(A, NROWA, NCOLA, TRANSA, B, NROWB, NCOLB, TRANSB, C)
C     ===================================================================
C
      SUBROUTINE GMMATS (A, NROWA, NCOLA, ITRANSA,
     1                   B, NROWB, NCOLB, ITRANSB, C)
      IMPLICIT NONE
      INTEGER NROWA, NCOLA, ITRANSA, NROWB, NCOLB, ITRANSB
      REAL    A(NROWA, NCOLA), B(NROWB, NCOLB), C(*)
C
      INTEGER I, J, K, NRA, NCA, NRB, NCB, IDX
      REAL    SUM
C
C     Handle transpose flags
C     ITRANSA=0: use A as-is; ITRANSA=1: use A^T
C     ITRANSB=0: use B as-is; ITRANSB=1: use B^T
C
      IF (ITRANSA .EQ. 0) THEN
        NRA = NROWA
        NCA = NCOLA
      ELSE
        NRA = NCOLA
        NCA = NROWA
      END IF
C
      IF (ITRANSB .EQ. 0) THEN
        NRB = NROWB
        NCB = NCOLB
      ELSE
        NRB = NCOLB
        NCB = NROWB
      END IF
C
C     Result is NRA x NCB
      IDX = 0
      DO 100 J = 1, NCB
        DO 90 I = 1, NRA
          SUM = 0.0
          DO 80 K = 1, NCA
            IF (ITRANSA .EQ. 0) THEN
              IF (ITRANSB .EQ. 0) THEN
                SUM = SUM + A(I,K) * B(K,J)
              ELSE
                SUM = SUM + A(I,K) * B(J,K)
              END IF
            ELSE
              IF (ITRANSB .EQ. 0) THEN
                SUM = SUM + A(K,I) * B(K,J)
              ELSE
                SUM = SUM + A(K,I) * B(J,K)
              END IF
            END IF
   80     CONTINUE
          IDX = IDX + 1
          C(IDX) = SUM
   90   CONTINUE
  100 CONTINUE
C
      RETURN
      END
