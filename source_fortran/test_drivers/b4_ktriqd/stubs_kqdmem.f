C=======================================================================
C     Stubs for the KQDMEM/KTRMEM (B4, CQUAD4 membrane) test driver.
C
C     Lets the UNMODIFIED kqdmem.f + ktrmem.f run outside NASTRAN so
C     their output can serve as benchmark B4's reference, replacing the
C     analytical two-CST assembly computed in Python.
C
C     /MATOUT/ uses KTRMEM's view (G11..G33, GSUBE at word 12). KQDMEM
C     declares the same block as DUM99(11), GSUBE, DUM88(6), which puts
C     GSUBE at word 12 as well, so the two views agree.
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
C     Same material the analytical reference used for B4.
C       E = 200e9 Pa, nu = 0.3, plane stress
C
      E  = 200.0E9
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
C     SMA1B captures the four 3x3 blocks KQDMEM emits per pivot.
C
C     KQDMEM sums the sub-triangle contributions into KSUM, copies each
C     3x3 into the leading rows of a 6x6 KIJ (words 1,2,3 / 7,8,9 /
C     13,14,15, i.e. stored BY ROWS) and calls
C         SMA1B (KIJ(1), NGRID(J), -1, IFKGG, 0.0D0)   J = 1..4
C     Accumulated into a 12x12 (4 nodes x 3 translational DOF).
C=======================================================================

      SUBROUTINE SMA1B (KIJ, IGRID, ICODE, IFFILE, FACTOR)
      IMPLICIT NONE
      DOUBLE PRECISION KIJ(36), FACTOR
      INTEGER IGRID, ICODE, IFFILE
C
      DOUBLE PRECISION KSTORE(12,12)
      INTEGER          NPVTG
      COMMON /KCAP  /  KSTORE, NPVTG
C
      INTEGER I, J, IR, IC, IROW, ICOL
C
C     Grid ids are 101..104 -> node index 1..4
      IROW = NPVTG - 100
      ICOL = IGRID - 100
      IF (IROW .LT. 1 .OR. IROW .GT. 4) RETURN
      IF (ICOL .LT. 1 .OR. ICOL .GT. 4) RETURN
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
C     Identity. All four CSIDs are 0 in this test case.
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
C     HMAT is referenced from KTRMEM's HEAT branch, never entered here.
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
