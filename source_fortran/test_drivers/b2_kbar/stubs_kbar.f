C=======================================================================
C     Stubs for the KBAR (B2, CBAR) standalone test driver.
C
C     Lets the UNMODIFIED nastran/NASTRAN-95/mis/kbar.f run outside
C     NASTRAN so its output can serve as benchmark B2's reference,
C     replacing the analytical Euler-Bernoulli matrix computed in Python.
C
C     /MATOUT/ layout here is KBAR's:
C       (E, G, NU, RHO, ALPHA, TSUBO, GSUBE, SIGT, SIGC, SIGS)
C     which is NOT the (G11,G12,...) layout KTRMEM uses.
C=======================================================================

      SUBROUTINE MAT (IELID)
      IMPLICIT NONE
      INTEGER IELID
C
      INTEGER   MATIDC, MATFLG
      REAL      ELTEMP, STRESS
      COMMON /MATIN / MATIDC, MATFLG, ELTEMP, STRESS
C
      REAL      E, G, NU, RHO, ALPHA, TSUBO, GSUBE, SIGT, SIGC, SIGS
      COMMON /MATOUT/ E, G, NU, RHO, ALPHA, TSUBO, GSUBE,
     1                SIGT, SIGC, SIGS
C
C     Same constants the analytical reference used.
      E     = 200.0E9
      G     = 76923000000.0
      NU    = 0.3
      RHO   = 0.0
      ALPHA = 0.0
      TSUBO = 0.0
      GSUBE = 0.0
      SIGT  = 0.0
      SIGC  = 0.0
      SIGS  = 0.0
      RETURN
      END


C=======================================================================
C     SMA1B captures each 6x6 block KBAR produces.
C
C     KBAR runs two passes per pivot:
C       pass 1  -> (pivot, pivot)
C       pass 2  -> (pivot, far node)
C     each ending in CALL SMA1B (KEP(109), INDEX, -1, IFKGG, 0.0D0),
C     where KEP(109..144) is a 6x6 stored BY ROWS and INDEX is the SIL
C     of the column node. Accumulated into a 12x12 (2 nodes x 6 DOF).
C=======================================================================

      SUBROUTINE SMA1B (KEP, IGRID, ICODE, IFFILE, FACTOR)
      IMPLICIT NONE
      DOUBLE PRECISION KEP(36), FACTOR
      INTEGER IGRID, ICODE, IFFILE
C
      DOUBLE PRECISION KSTORE(12,12)
      INTEGER          NPVTG
      COMMON /KCAP  /  KSTORE, NPVTG
C
      INTEGER I, J, IR, IC, IROW, ICOL
C
C     SILs are 101,102 -> node index 1,2
      IROW = NPVTG - 100
      ICOL = IGRID - 100
      IF (IROW .LT. 1 .OR. IROW .GT. 2) RETURN
      IF (ICOL .LT. 1 .OR. ICOL .GT. 2) RETURN
C
      DO 20 I = 1, 6
        DO 10 J = 1, 6
          IR = 6*(IROW-1) + I
          IC = 6*(ICOL-1) + J
          KSTORE(IR,IC) = KSTORE(IR,IC) + KEP(6*(I-1) + J)
   10   CONTINUE
   20 CONTINUE
      RETURN
      END


      SUBROUTINE TRANSD (RCSID, T)
      IMPLICIT NONE
      REAL RCSID
      DOUBLE PRECISION T(9)
      INTEGER I
C     Identity. Both grid points are in basic coordinates in this test
C     case (MCSIDA = MCSIDB = 0) so KBAR never calls this; present so
C     the link resolves.
      DO 10 I = 1, 9
        T(I) = 0.0D0
   10 CONTINUE
      T(1) = 1.0D0
      T(5) = 1.0D0
      T(9) = 1.0D0
      RETURN
      END


      SUBROUTINE MESAGE (ICODE, IMSG, IPARM)
      IMPLICIT NONE
      INTEGER ICODE, IMSG
      REAL    IPARM
      WRITE (6,900) ICODE, IMSG
  900 FORMAT (' MESAGE called: code=', I6, ' msg=', I6)
      RETURN
      END


C=======================================================================
C     HMAT is only referenced from KBAR's HEAT branch (statement 2000),
C     never entered here because HEAT = .FALSE.
C=======================================================================
      SUBROUTINE HMAT (IELID)
      IMPLICIT NONE
      INTEGER IELID
      REAL    FK
      COMMON /HMTOUT/ FK
      WRITE (6,*) 'WARNING: HMAT reached but HEAT should be .FALSE.'
      FK = 50.0
      RETURN
      END
