C     TEST DRIVER FOR KELBOW - CURVED BAR / ELBOW ELEMENT
C     =====================================================
C     This program calls KELBOW twice (once per grid point as NPVT)
C     to assemble the full 12x12 stiffness matrix for a 90-degree elbow.
C
C     Test geometry: 90-degree elbow (quarter circle) in the XZ plane
C       Grid A at (1, 0, 0) - start of elbow
C       Grid B at (0, 0, 1) - end of elbow (quarter circle in XZ)
C       Reference vector V = (0, 1, 0) - perpendicular to plane of curve
C       Radius of curvature R = 1.0
C       Angle BETAR = 90 degrees
C
C     Material: E = 200e9, NU = 0.3, G = E/(2*(1+NU))
C     Cross section: circular pipe
C       Outer radius r_o = 0.05, wall thickness t = 0.005
C       A = pi*(r_o^2 - r_i^2) ~ 0.001414
C       I1 = I2 = pi/4*(r_o^4 - r_i^4) ~ 1.688e-6
C       FJ = pi/2*(r_o^4 - r_i^4) ~ 3.376e-6
C
      PROGRAM TESTKELBOW
      IMPLICIT NONE
C
C     COMMON blocks matching KELBOW expectations
C
C     /SMA1ET/ - Element Connection Property Table
      REAL             ECPT(100)
      INTEGER          NECPT(100)
      EQUIVALENCE     (NECPT(1), ECPT(1))
      COMMON /SMA1ET/ ECPT
C
C     /SMA1CL/ - Control parameters
      INTEGER          IOPT4, K4GGSW, NPVT, LEFT, FROWIC, LROWIC
      INTEGER          NROWSC, TNROWS, JMAX, NLINKS, LINK(10)
      INTEGER          IDETCK, DODET, NOGO_CL
      COMMON /SMA1CL/ IOPT4, K4GGSW, NPVT, LEFT, FROWIC, LROWIC,
     1                 NROWSC, TNROWS, JMAX, NLINKS, LINK, IDETCK,
     2                 DODET, NOGO_CL
C
C     /SMA1HT/ - Heat flag
      LOGICAL          HEAT
      COMMON /SMA1HT/ HEAT
C
C     /MATIN/ - Material input
      INTEGER          MATIDC, MATFLG
      REAL             ELTEMP, STRESS, SINTH, COSTH
      COMMON /MATIN / MATIDC, MATFLG, ELTEMP, STRESS, SINTH, COSTH
C
C     /MATOUT/ - Material output
      REAL             E, GG, NU, RHO, ALPHA, TSUB0, GSUBE
      REAL             SIGT, SIGC, SIGS
      COMMON /MATOUT/ E, GG, NU, RHO, ALPHA, TSUB0, GSUBE,
     1                SIGT, SIGC, SIGS
C
C     /HMTOUT/ - Heat material output
      REAL             FK
      COMMON /HMTOUT/ FK
C
C     /SMA1IO/ - I/O parameters
      INTEGER          IFCSTM, IFMPT, IFDIT, IDUM1, IFECPT, IGECPT
      INTEGER          IFGPCT, IGGPCT, IFGEI, IGGEI, IFKGG, IGKGG
      INTEGER          IF4GG, IG4GG
      INTEGER          IFGPST, IGGPST, INRW, OUTRW, CLSNRW, CLSRW
      INTEGER          NEOR, EOR
      INTEGER          MCBKGG(7), MCB4GG(7)
      COMMON /SMA1IO/ IFCSTM, IFMPT, IFDIT, IDUM1, IFECPT, IGECPT,
     1                 IFGPCT, IGGPCT, IFGEI, IGGEI, IFKGG, IGKGG,
     2                 IF4GG, IG4GG,
     3                 IFGPST, IGGPST, INRW, OUTRW, CLSNRW, CLSRW,
     4                 NEOR, EOR, MCBKGG, MCB4GG
C
C     /SMA1DP/ - Double precision work arrays for KELBOW
C     KE(144), KEP(144), DELA(6), DELB(6) = 300 DP words
      DOUBLE PRECISION KE(144), KEP(144), DELA(6), DELB(6)
      COMMON /SMA1DP/ KE, KEP, DELA, DELB
C
C     /SMA1BK/ - Block pointers (not used but must exist)
      INTEGER          ICSTM, NCSTM, IGPCT, NGPCT, IPOINT, NPOINT
      INTEGER          I6X6K, N6X6K, I6X64, N6X64
      COMMON /SMA1BK/ ICSTM, NCSTM, IGPCT, NGPCT, IPOINT, NPOINT,
     1                 I6X6K, N6X6K, I6X64, N6X64
C
C     /SYSTEM/ - System parameters
      INTEGER          SYSBUF, NOUT
      COMMON /SYSTEM/ SYSBUF, NOUT
C
C     /MACHIN/ - Machine type
      INTEGER          MACH
      COMMON /MACHIN/ MACH
C
C     /ZZZZZZ/ - Open core (provide some workspace)
      DOUBLE PRECISION Z(1000)
      COMMON /ZZZZZZ/ Z
C
C     Local variables
      DOUBLE PRECISION KFULL(12,12)
      DOUBLE PRECISION ROWSUM, TRACE, MAXVAL, OFFDIAG
      INTEGER          IPVT, I, J, IROW, ICOL
      REAL             PI, RO, RI, TWALL, AREA, AI1, AI2, AFJ
C
C     === INITIALIZATION ===
C
      MACH = 1
      NOUT = 6
      SYSBUF = 1
      HEAT = .FALSE.
      IOPT4 = 0
      K4GGSW = 0
      NOGO_CL = 0
C
C     Zero out SMA1IO
      IFCSTM = 0
      IFMPT = 0
      IFDIT = 0
      IDUM1 = 0
      IFECPT = 0
      IGECPT = 0
      IFGPCT = 0
      IGGPCT = 0
      IFGEI = 0
      IGGEI = 0
      IFKGG = 0
      IGKGG = 0
      IF4GG = 0
      IG4GG = 0
      IFGPST = 0
      IGGPST = 0
      INRW = 0
      OUTRW = 0
      CLSNRW = 0
      CLSRW = 0
      NEOR = 0
      EOR = 0
      DO 3 I = 1, 7
        MCBKGG(I) = 0
        MCB4GG(I) = 0
    3 CONTINUE
C
C     Zero SMA1BK
      ICSTM = 0
      NCSTM = 0
      IGPCT = 0
      NGPCT = 0
      IPOINT = 0
      NPOINT = 0
      I6X6K = 0
      N6X6K = 0
      I6X64 = 0
      N6X64 = 0
C
C     Zero workspace
      DO 4 I = 1, 1000
        Z(I) = 0.0D0
    4 CONTINUE
C
C     Set material properties
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
C     === COMPUTE CROSS-SECTION PROPERTIES ===
      PI = 3.14159265358979
      RO = 0.05
      TWALL = 0.005
      RI = RO - TWALL
      AREA = PI * (RO**2 - RI**2)
      AI1 = PI / 4.0 * (RO**4 - RI**4)
      AI2 = AI1
      AFJ = PI / 2.0 * (RO**4 - RI**4)
C
C     === FILL ECPT ===
      DO 10 I = 1, 100
        ECPT(I) = 0.0
   10 CONTINUE
C
C     ECPT(1) = Element ID
      NECPT(1) = 101
C     ECPT(2) = SIL grid A
      NECPT(2) = 1
C     ECPT(3) = SIL grid B
      NECPT(3) = 2
C     ECPT(4-6) = Reference vector V = (0, 1, 0)
      ECPT(4) = 0.0
      ECPT(5) = 1.0
      ECPT(6) = 0.0
C     ECPT(7) = Coord sys ID for V (basic = 0)
      NECPT(7) = 0
C     ECPT(8) = Material ID
      NECPT(8) = 1
C     ECPT(9) = Area
      ECPT(9) = AREA
C     ECPT(10) = I1
      ECPT(10) = AI1
C     ECPT(11) = I2
      ECPT(11) = AI2
C     ECPT(12) = FJ (torsional constant)
      ECPT(12) = AFJ
C     ECPT(13) = NSM (non-structural mass)
      ECPT(13) = 0.0
C     ECPT(14) = FE (force element)
      ECPT(14) = 0.0
C     ECPT(15-22) = Stress recovery points (set to 0)
C     (already zeroed)
C     ECPT(23) = K1 (shear area factor 1)
      ECPT(23) = 0.9
C     ECPT(24) = K2 (shear area factor 2)
      ECPT(24) = 0.9
C     ECPT(25) = C (stress intensification factor)
      ECPT(25) = 1.0
C     ECPT(26) = KX (flexibility correction)
      ECPT(26) = 1.0
C     ECPT(27) = KY (flexibility correction)
      ECPT(27) = 1.0
C     ECPT(28) = KZ (flexibility correction)
      ECPT(28) = 1.0
C     ECPT(29) = R (radius of curvature)
      ECPT(29) = 1.0
C     ECPT(30) = BETAR (angle from GA to GB in DEGREES for trig)
C     NOTE: KELBOW uses DTR(BETAR) = BETAR * 0.017453292 for radians
C     and SID/COD compute sin/cos using degrees. So BETAR = 90.0
      ECPT(30) = 90.0
C     ECPT(31) = MCSIDA (coord sys ID for A, basic = 0)
      NECPT(31) = 0
C     ECPT(32-34) = Grid A coordinates (1, 0, 0)
      ECPT(32) = 1.0
      ECPT(33) = 0.0
      ECPT(34) = 0.0
C     ECPT(35) = MCSIDB (coord sys ID for B, basic = 0)
      NECPT(35) = 0
C     ECPT(36-38) = Grid B coordinates (0, 0, 1)
      ECPT(36) = 0.0
      ECPT(37) = 0.0
      ECPT(38) = 1.0
C     ECPT(39) = Element temperature
      ECPT(39) = 20.0
C
C     === ZERO THE FULL STIFFNESS MATRIX ===
      DO 20 I = 1, 12
        DO 15 J = 1, 12
          KFULL(I,J) = 0.0D0
   15   CONTINUE
   20 CONTINUE
C
C     === PRINT HEADER ===
      WRITE(6,*) 'KELBOW Test Driver - 90-degree Elbow Element'
      WRITE(6,*) '============================================='
      WRITE(6,*) ''
      WRITE(6,*) 'Geometry:'
      WRITE(6,*) '  Grid A: (1.0, 0.0, 0.0)'
      WRITE(6,*) '  Grid B: (0.0, 0.0, 1.0)'
      WRITE(6,*) '  Ref vector V: (0.0, 1.0, 0.0)'
      WRITE(6,*) '  Radius R = 1.0'
      WRITE(6,*) '  Angle = 90 degrees'
      WRITE(6,*) ''
      WRITE(6,*) 'Material:'
      WRITE(6,900) E
      WRITE(6,901) NU
      WRITE(6,902) GG
  900 FORMAT('  E  = ', E12.5)
  901 FORMAT('  NU = ', F6.3)
  902 FORMAT('  G  = ', E12.5)
      WRITE(6,*) ''
      WRITE(6,*) 'Cross section (circular pipe):'
      WRITE(6,903) AREA
      WRITE(6,904) AI1
      WRITE(6,905) AI2
      WRITE(6,906) AFJ
  903 FORMAT('  A   = ', E12.5)
  904 FORMAT('  I1  = ', E12.5)
  905 FORMAT('  I2  = ', E12.5)
  906 FORMAT('  FJ  = ', E12.5)
      WRITE(6,*) ''
C
C     === CALL KELBOW TWICE: PIVOT=A, THEN PIVOT=B ===
C
      WRITE(6,*) 'Calling KELBOW with NPVT = Grid A (SIL=1)...'
      NPVT = NECPT(2)
      CALL KELBOW
      WRITE(6,*) '  Done (pivot A).'
C
      WRITE(6,*) 'Calling KELBOW with NPVT = Grid B (SIL=2)...'
      NPVT = NECPT(3)
      CALL KELBOW
      WRITE(6,*) '  Done (pivot B).'
      WRITE(6,*) ''
C
C     === RETRIEVE ASSEMBLED 12x12 MATRIX ===
      CALL GET_KMATRIX_ELBOW(KFULL)
C
C     === PRINT THE FULL 12x12 STIFFNESS MATRIX ===
      WRITE(6,*) 'Full 12x12 Stiffness Matrix (K_global)'
      WRITE(6,*) 'Rows 1-6: Grid A (Fx,Fy,Fz,Mx,My,Mz)'
      WRITE(6,*) 'Rows 7-12: Grid B (Fx,Fy,Fz,Mx,My,Mz)'
      WRITE(6,*) '=========================================='
      DO 200 I = 1, 12
        WRITE(6,300) I, (KFULL(I,J), J=1,6)
        WRITE(6,301) (KFULL(I,J), J=7,12)
  200 CONTINUE
  300 FORMAT('Row',I3,':',6(1X,E13.6))
  301 FORMAT('      ',6(1X,E13.6))
      WRITE(6,*) ''
C
C     === PRINT BY 6x6 BLOCKS ===
      WRITE(6,*) 'Stiffness matrix by 6x6 sub-blocks:'
      WRITE(6,*) ''
      WRITE(6,*) 'K(A,A) - Self-coupling of Grid A:'
      DO 410 I = 1, 6
        WRITE(6,460) (KFULL(I,J), J=1,6)
  410 CONTINUE
      WRITE(6,*) ''
      WRITE(6,*) 'K(A,B) - Coupling A to B:'
      DO 420 I = 1, 6
        WRITE(6,460) (KFULL(I,J), J=7,12)
  420 CONTINUE
      WRITE(6,*) ''
      WRITE(6,*) 'K(B,A) - Coupling B to A:'
      DO 430 I = 7, 12
        WRITE(6,460) (KFULL(I,J), J=1,6)
  430 CONTINUE
      WRITE(6,*) ''
      WRITE(6,*) 'K(B,B) - Self-coupling of Grid B:'
      DO 440 I = 7, 12
        WRITE(6,460) (KFULL(I,J), J=7,12)
  440 CONTINUE
  460 FORMAT(6(1X,E16.9))
      WRITE(6,*) ''
C
C     === SYMMETRY CHECK ===
      WRITE(6,*) '=== VERIFICATION CHECKS ==='
      WRITE(6,*) ''
      MAXVAL = 0.0D0
      DO 500 I = 1, 12
        DO 500 J = 1, 12
          IF (DABS(KFULL(I,J)) .GT. MAXVAL) MAXVAL = DABS(KFULL(I,J))
  500 CONTINUE
C
      OFFDIAG = 0.0D0
      DO 510 I = 1, 12
        DO 510 J = I+1, 12
          IF (DABS(KFULL(I,J)-KFULL(J,I)) .GT. OFFDIAG) THEN
            OFFDIAG = DABS(KFULL(I,J)-KFULL(J,I))
          END IF
  510 CONTINUE
      WRITE(6,520) MAXVAL
      WRITE(6,521) OFFDIAG
      IF (MAXVAL .GT. 0.0D0) THEN
        WRITE(6,522) OFFDIAG/MAXVAL
      END IF
  520 FORMAT('1. Symmetry check:  max |K_ij| = ', E12.5)
  521 FORMAT('   max |K_ij - K_ji|           = ', E12.5)
  522 FORMAT('   Relative asymmetry          = ', E12.5)
      WRITE(6,*) ''
C
C     === EQUILIBRIUM CHECK (row sums for rigid body modes) ===
C     For a proper stiffness matrix, rigid body translation should
C     give zero forces: sum of each row in translation DOFs = 0
      WRITE(6,*) '2. Equilibrium check (row sums for rigid body):'
      DO 600 I = 1, 12
        ROWSUM = 0.0D0
        DO 590 J = 1, 12
          ROWSUM = ROWSUM + KFULL(I,J)
  590   CONTINUE
        WRITE(6,610) I, ROWSUM
  600 CONTINUE
  610 FORMAT('   Row ',I2,' sum = ', E12.5)
      WRITE(6,*) ''
C
C     === POSITIVE DIAGONAL CHECK ===
      WRITE(6,*) '3. Diagonal values (should be positive):'
      TRACE = 0.0D0
      DO 700 I = 1, 12
        WRITE(6,710) I, I, KFULL(I,I)
        TRACE = TRACE + KFULL(I,I)
  700 CONTINUE
  710 FORMAT('   K(',I2,',',I2,') = ', E16.9)
      WRITE(6,720) TRACE
  720 FORMAT('   Trace = ', E16.9)
      WRITE(6,*) ''
C
C     === K(A,B) = K(B,A)^T CHECK ===
      WRITE(6,*) '4. Off-diagonal block transpose check:'
      WRITE(6,*) '   K(A,B) should equal K(B,A)^T'
      OFFDIAG = 0.0D0
      DO 800 I = 1, 6
        DO 800 J = 1, 6
          IF (DABS(KFULL(I,J+6)-KFULL(J+6,I)) .GT. OFFDIAG) THEN
            OFFDIAG = DABS(KFULL(I,J+6)-KFULL(J+6,I))
          END IF
  800 CONTINUE
      WRITE(6,810) OFFDIAG
  810 FORMAT('   max |K_AB(i,j) - K_BA(j,i)| = ', E12.5)
      WRITE(6,*) ''
C
      IF (NOGO_CL .NE. 0) THEN
        WRITE(6,*) 'WARNING: NOGO flag was set - errors occurred!'
      ELSE
        WRITE(6,*) 'Computation completed successfully.'
      END IF
C
      STOP
      END
