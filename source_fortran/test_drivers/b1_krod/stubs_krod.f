C=======================================================================
C     Stubs for the KROD (B1, CROD) standalone test driver.
C
C     Lets the UNMODIFIED nastran/NASTRAN-95/mis/krod.f run outside
C     NASTRAN so its output can serve as benchmark B1's reference,
C     replacing the analytical EA/L + GJ/L matrix computed in Python.
C
C     NOTE on /MATOUT/: NASTRAN reuses this COMMON block with different
C     layouts depending on INFLAG/MATFLG. KROD reads it as
C     (E, G, NU, RHO, ALPHA, TSUBO, GSUBE, SIGT, SIGC, SIGS), which is
C     NOT the (G11,G12,...) layout KTRMEM uses. The layout below must
C     match KROD's declaration, hence a per-driver stub.
C=======================================================================

      SUBROUTINE MAT (IELID)
      IMPLICIT NONE
      INTEGER IELID
C
      INTEGER   MATIDC, MATFLG
      REAL      ELTEMP, STRESS, SINTH, COSTH
      COMMON /MATIN / MATIDC, MATFLG, ELTEMP, STRESS, SINTH, COSTH
C
      REAL      E, G, NU, RHO, ALPHA, TSUBO, GSUBE, SIGT, SIGC, SIGS
      COMMON /MATOUT/ E, G, NU, RHO, ALPHA, TSUBO, GSUBE,
     1                SIGT, SIGC, SIGS
C
C     Same constants the analytical reference used, so that any residual
C     difference is attributable to the formulation rather than to the
C     material data.
C       E = 200e9 Pa, G = 76.923e9 Pa
C
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
C     SMA1B captures each 6x6 block KROD produces.
C
C     KROD calls SMA1B(KE, ECPT(NONPVT), -1, IFKGG, 0.0D0) twice per
C     pivot: first with NONPVT=2 (the pivot itself, positive constants),
C     then with NONPVT=3 (the far node, constants negated).
C
C     KE is 6x6 stored BY ROWS. Only two diagonal 3x3 blocks are filled:
C       translational  KE(1..3), KE(7..9), KE(13..15)     scaled by A*E/L
C       rotational     KE(22..24), KE(28..30), KE(34..36) scaled by J*G/L
C     Accumulated into a 12x12 (2 nodes x 6 DOF).
C=======================================================================

      SUBROUTINE SMA1B (KE, IGRID, ICODE, IFFILE, FACTOR)
      IMPLICIT NONE
      DOUBLE PRECISION KE(36), FACTOR
      INTEGER IGRID, ICODE, IFFILE
C
      DOUBLE PRECISION KSTORE(12,12)
      INTEGER          NPVTG
      COMMON /KCAP  /  KSTORE, NPVTG
C
      INTEGER I, J, IR, IC, IROW, ICOL
C
C     Grid ids are 101,102 -> node index 1,2
      IROW = NPVTG - 100
      ICOL = IGRID - 100
      IF (IROW .LT. 1 .OR. IROW .GT. 2) RETURN
      IF (ICOL .LT. 1 .OR. ICOL .GT. 2) RETURN
C
      DO 20 I = 1, 6
        DO 10 J = 1, 6
          IR = 6*(IROW-1) + I
          IC = 6*(ICOL-1) + J
          KSTORE(IR,IC) = KSTORE(IR,IC) + KE(6*(I-1) + J)
   10   CONTINUE
   20 CONTINUE
      RETURN
      END


      SUBROUTINE TRANSD (RCSID, T)
      IMPLICIT NONE
      REAL RCSID
      DOUBLE PRECISION T(9)
      INTEGER I
C     Identity. Every CSID in this test case is 0 so KROD skips the
C     transform entirely; present so the link resolves.
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
C     HMAT is only referenced from KROD's HEAT branch (statement 400),
C     which this driver never enters because HEAT = .FALSE.
C=======================================================================
      SUBROUTINE HMAT (NECPT)
      IMPLICIT NONE
      INTEGER NECPT(4)
      REAL    FK
      COMMON /HMTOUT/ FK
      WRITE (6,*) 'WARNING: HMAT reached but HEAT should be .FALSE.'
      FK = 50.0
      RETURN
      END
