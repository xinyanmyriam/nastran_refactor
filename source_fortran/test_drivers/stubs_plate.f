C     STUB SUBROUTINES FOR PLATE BENDING TEST DRIVERS (KTRPLT/KQDPLT)
C     ==================================================================
C
C     These stubs replace NASTRAN infrastructure subroutines that are
C     not needed for standalone testing of the plate bending elements.
C     The plate elements call KTRBSC which in turn calls:
C       MAT, TRANSD, SMA1B, INVERD, GMMATD, MESAGE
C
C     SMA1B for plates uses 6x6 blocks (3 per pivot for TRPLT, 
C     4 per pivot for QDPLT).
C
C     ===================================================================
C     CAPTURE BUFFER: stores all 6x6 sub-blocks output by the plate
C     elements via SMA1B calls. 
C     For TRPLT: 3 pivots x 3 grids = 9 blocks of 6x6
C     For QDPLT: 4 pivots x 4 grids = 16 blocks of 6x6
C     We use the larger buffer (16 blocks).
C     ===================================================================
C
      BLOCK DATA CAPTUREDATAP
      IMPLICIT NONE
      DOUBLE PRECISION KSTORE(6,6,4,4)
      INTEGER          CURPVT, NPVTS
      COMMON /CAPTUREP/ KSTORE, CURPVT, NPVTS
      DATA CURPVT /0/
      DATA NPVTS /0/
      DATA KSTORE /576*0.0D0/
      END
C
C     ===================================================================
C     MAT - Material property routine stub for plate bending
C     Sets material properties in COMMON /MATOUT/
C     For INFLAG=2: returns G-matrix for plate bending
C     ===================================================================
C
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
      REAL      E, NU
C
C     Material: E = 200e9 Pa, nu = 0.3
      E = 200.0E9
      NU = 0.3
C
      IF (INFLAG .EQ. 2) THEN
C       Plate bending G-matrix (for use with moment of inertia I)
C       G11 = E/(1-nu^2), G12 = nu*G11, G22 = G11
C       G33 = G11*(1-nu)/2, G13=G23=0
        G11 = E / (1.0 - NU*NU)
        G12 = NU * G11
        G13 = 0.0
        G22 = G11
        G23 = 0.0
        G33 = G11 * (1.0 - NU) / 2.0
      ELSE IF (INFLAG .EQ. 3) THEN
C       Transverse shear G2x2 matrix (for T2 coupling)
C       Return zeros to indicate no transverse shear
        G2X211 = 0.0
        G2X212 = 0.0
        G2X222 = 0.0
      ELSE
C       Default: isotropic 3D
        G11 = E / (1.0 - NU*NU)
        G12 = NU * G11
        G13 = 0.0
        G22 = G11
        G23 = 0.0
        G33 = G11 * (1.0 - NU) / 2.0
      END IF
C
      RHO = 7800.0
      ALPHA1 = 12.0E-6
      ALPHA2 = 12.0E-6
      ALP12 = 0.0
      TSUB0 = 20.0
      GSUBE = 0.0
      SIGTEN = 0.0
      SIGCOM = 0.0
      SIGSHE = 0.0
C
      RETURN
      END
C
C     ===================================================================
C     SMA1B - Stiffness matrix assembly stub for plate bending
C     Captures the 6x6 KIJ sub-block for each grid point
C     ===================================================================
C
      SUBROUTINE SMA1B (KIJ, IGRID, ICODE, IFFILE, FACTOR)
      IMPLICIT NONE
      DOUBLE PRECISION KIJ(6,6)
      INTEGER          IGRID, ICODE, IFFILE
      DOUBLE PRECISION FACTOR
C
      DOUBLE PRECISION KSTORE(6,6,4,4)
      INTEGER          CURPVT, NPVTS
      COMMON /CAPTUREP/ KSTORE, CURPVT, NPVTS
C
      REAL             ECPT(100)
      INTEGER          NECPT(100)
      EQUIVALENCE     (NECPT(1), ECPT(1))
      COMMON /SMA1ET/ ECPT
C
      INTEGER          IOPT4, K4GGSW, NPVT
      REAL             DUMCL(7)
      INTEGER          LINK(10), IDETCK
      REAL             DODET
      INTEGER          NOGO
      COMMON /SMA1CL/ IOPT4, K4GGSW, NPVT, DUMCL, LINK, IDETCK,
     1                DODET, NOGO
C
      INTEGER I, J, IPVT, IGRP
C
C     Only capture the primary stiffness (FACTOR = 0)
C
      IF (FACTOR .NE. 0.0D0) RETURN
C
C     Debug output
      WRITE(6,*) '  SMA1B: NPVT=',NPVT,' IGRID=',IGRID,' NPVTS=',NPVTS
C
C     Determine which pivot point (row block) this is
C     For TRPLT: grid points are at NECPT(2), NECPT(3), NECPT(4)
C     For QDPLT: grid points are at NECPT(2), NECPT(3), NECPT(4), NECPT(5)
C
      DO 10 IPVT = 1, NPVTS
        IF (NPVT .EQ. NECPT(IPVT+1)) GO TO 20
   10 CONTINUE
      WRITE(6,*) 'SMA1B: Could not find pivot NPVT=', NPVT
      WRITE(6,*) '  NECPT(2..5)=',NECPT(2),NECPT(3),NECPT(4),NECPT(5)
      RETURN
C
C     Determine which grid point (column block) this is
C
   20 DO 30 IGRP = 1, NPVTS
        IF (IGRID .EQ. NECPT(IGRP+1)) GO TO 40
   30 CONTINUE
      WRITE(6,*) 'SMA1B: Could not find grid IGRID=', IGRID
      RETURN
C
C     Store the 6x6 sub-block
C
   40 DO 50 I = 1, 6
        DO 50 J = 1, 6
          KSTORE(I, J, IPVT, IGRP) = KIJ(I, J)
   50 CONTINUE
C
      RETURN
      END
C
C     ===================================================================
C     GET_KMATRIX_PLATE - Retrieve the assembled stiffness matrix
C     For N nodes: returns N*6 x N*6 matrix (but only 3 DOF per node
C     are used for plate bending: w, theta_x, theta_y in positions 1-3
C     of each 6-DOF block)
C     ===================================================================
C
      SUBROUTINE GET_KMATRIX_PLATE (KFULL, NDOF)
      IMPLICIT NONE
      INTEGER NDOF
      DOUBLE PRECISION KFULL(NDOF, NDOF)
C
      DOUBLE PRECISION KSTORE(6,6,4,4)
      INTEGER          CURPVT, NPVTS
      COMMON /CAPTUREP/ KSTORE, CURPVT, NPVTS
C
      INTEGER I, J, IPVT, IGRP, IROW, ICOL
C
C     Zero out
      DO 5 I = 1, NDOF
        DO 5 J = 1, NDOF
          KFULL(I, J) = 0.0D0
    5 CONTINUE
C
C     Assemble from stored 6x6 blocks
C     Each 6x6 block goes into the (IPVT, IGRP) position
C     Only first 3x3 of each 6x6 contains plate bending DOFs
C     (w, theta_x, theta_y) since NASTRAN stores in displacement form.
C     Actually for plate bending elements, SMA1B receives a full 6x6
C     where the plate bending DOFs are in positions (1,2,3) = (w,thx,thy).
C     Positions 4,5,6 are membrane DOFs (u,v,theta_z) which are zero
C     for pure bending elements.
C
      DO 100 IPVT = 1, NPVTS
        DO 90 IGRP = 1, NPVTS
          DO 80 I = 1, 6
            DO 80 J = 1, 6
              IROW = (IPVT-1)*6 + I
              ICOL = (IGRP-1)*6 + J
              IF (IROW .LE. NDOF .AND. ICOL .LE. NDOF) THEN
                KFULL(IROW, ICOL) = KSTORE(I, J, IPVT, IGRP)
              END IF
   80       CONTINUE
   90     CONTINUE
  100 CONTINUE
C
      RETURN
      END
C
C     ===================================================================
C     TRANSD - Coordinate transformation stub
C     Returns identity since all CSID = 0 for our test case
C     ===================================================================
C
      SUBROUTINE TRANSD (ICSID, T)
      IMPLICIT NONE
      INTEGER ICSID
      DOUBLE PRECISION T(9)
C
C     Return identity matrix (all nodes in basic coordinate system)
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
      WRITE(6,*) 'MESAGE called: ICODE=', ICODE, ' IMSG=', IMSG
      WRITE(6,*) '  IPARM=', IPARM(1), IPARM(2)
      IF (ICODE .LT. 0) THEN
        WRITE(6,*) 'FATAL ERROR - stopping.'
        STOP
      END IF
      RETURN
      END
