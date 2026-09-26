C=======================================================================
C     Stubs for the KTRMEM (B3, CTRIA3 membrane) standalone test driver.
C
C     Purpose: let the UNMODIFIED nastran/NASTRAN-95/mis/ktrmem.f run
C     outside NASTRAN so that its output can serve as the reference for
C     benchmark B3. Previously B3 was scored against an analytical
C     textbook CST matrix computed in Python, which is a different
C     criterion from "reproduces the legacy implementation".
C
C     Provided here:
C       MAT     - material data; returns the plane-stress G matrix
C       SMA1B   - insertion routine; captures the 3x3 blocks KTRMEM emits
C       TRANSD  - coordinate transform; identity (all CSIDs are 0 anyway)
C       MESAGE  - error reporter
C     GMMATD is NOT stubbed: the real one is linked, because the storage
C     convention it implements is part of what is under test.
C=======================================================================

      SUBROUTINE MAT (IELID)
      IMPLICIT NONE
      INTEGER IELID
C
      INTEGER   MATID, INFLAG
      REAL      ELTEMP, STRESS, SINTH, COSTH
      COMMON /MATIN / MATID, INFLAG, ELTEMP, STRESS, SINTH, COSTH
C
      REAL      G11, G12, G13, G22, G23, G33, RHO
      REAL      ALPHA1, ALPHA2, ALP12, TSUB0, GSUBE
      REAL      SIGTEN, SIGCOM, SIGSHE
      REAL      G2X211, G2X212, G2X222
      COMMON /MATOUT/ G11, G12, G13, G22, G23, G33, RHO,
     1                ALPHA1, ALPHA2, ALP12, TSUB0, GSUBE,
     2                SIGTEN, SIGCOM, SIGSHE,
     3                G2X211, G2X212, G2X222
C
      REAL      E, NU, FAC
C
C     Same material the analytical reference used, so that any difference
C     between the two references is attributable to the formulation and
C     not to the material constants.
C       E = 2.1e11 Pa, nu = 0.3, plane stress
C
      E  = 2.1E11
      NU = 0.3
      FAC = E / (1.0 - NU*NU)
C
      G11 = FAC
      G12 = NU * FAC
      G13 = 0.0
      G22 = FAC
      G23 = 0.0
      G33 = FAC * (1.0 - NU) / 2.0
C
      RHO    = 0.0
      ALPHA1 = 0.0
      ALPHA2 = 0.0
      ALP12  = 0.0
      TSUB0  = 0.0
      GSUBE  = 0.0
      SIGTEN = 0.0
      SIGCOM = 0.0
      SIGSHE = 0.0
      G2X211 = 0.0
      G2X212 = 0.0
      G2X222 = 0.0
      RETURN
      END


C=======================================================================
C     SMA1B captures each 3x3 block KTRMEM produces.
C
C     KTRMEM calls SMA1B(KIJ, NECPT(I+1), -1, IFKGG, 0.0D0) once per
C     column-node I while NPVT selects the row node. KIJ is declared
C     6x6 but only the leading 3x3 is filled, stored BY ROWS:
C         KIJ(1),(2),(3)   = row 1, cols 1..3
C         KIJ(7),(8),(9)   = row 2
C         KIJ(13),(14),(15)= row 3
C     We accumulate into a 9x9 (3 nodes x 3 translational DOF).
C=======================================================================

      SUBROUTINE SMA1B (KIJ, IGRID, ICODE, IFFILE, FACTOR)
      IMPLICIT NONE
      DOUBLE PRECISION KIJ(36), FACTOR
      INTEGER IGRID, ICODE, IFFILE
C
      DOUBLE PRECISION KSTORE(9,9)
      INTEGER          NPVTG
      COMMON /KCAP  /  KSTORE, NPVTG
C
      INTEGER I, J, IR, IC, IROW, ICOL
C
C     Grid ids are 101,102,103 -> node index 1,2,3
      IROW = NPVTG - 100
      ICOL = IGRID - 100
      IF (IROW .LT. 1 .OR. IROW .GT. 3) RETURN
      IF (ICOL .LT. 1 .OR. ICOL .GT. 3) RETURN
C
      DO 20 I = 1, 3
        DO 10 J = 1, 3
          IR = 3*(IROW-1) + I
          IC = 3*(ICOL-1) + J
          KSTORE(IR,IC) = KSTORE(IR,IC) + KIJ(6*(I-1) + J)
   10   CONTINUE
   20 CONTINUE
      RETURN
      END


      SUBROUTINE TRANSD (ICSID, T)
      IMPLICIT NONE
      INTEGER ICSID
      DOUBLE PRECISION T(9)
      INTEGER I
C     Identity: every CSID in this test case is 0, so KTRMEM should not
C     even reach here. Present so the link succeeds.
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
C     HMAT supplies material data for HEAT problems. KTRMEM references it
C     from the HEAT branch (line ~1010), which this test case never takes
C     because HEAT = .FALSE. Present so the link resolves; if it is ever
C     reached that is a bug in the driver setup, so it says so loudly.
C=======================================================================

      SUBROUTINE HMAT (NECPT)
      IMPLICIT NONE
      INTEGER NECPT(4)
      REAL    MATBUF(7)
      COMMON /HMTOUT/ MATBUF
      WRITE (6,*) 'WARNING: HMAT reached but HEAT should be .FALSE.'
      MATBUF(1) = 50.0
      RETURN
      END
