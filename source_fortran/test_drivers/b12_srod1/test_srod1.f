C     TEST DRIVER FOR SROD1 - Rod Stress Recovery (B12)
C     ===================================================
C     This program sets up a rod element with known displacements
C     and calls SROD1 to compute stress recovery coefficients.
C     It then manually applies those coefficients to get stresses.
C
C     Test case:
C       Rod along x-axis, L=2.0, A=0.01, J=5.0e-6
C       Material: E=2.1e11, G=E/(2*(1+NU)), NU=0.3
C       Displacement: u_a = [0, 0, 0], u_b = [0.001, 0, 0]
C       Expected: axial_force = E*A*(u_b1-u_a1)/L
C                 axial_stress = force/A = E*(u_b1-u_a1)/L
C       Torsion: theta_a=0, theta_b=0.01 rad
C       Expected: torque = G*J*(theta_b-theta_a)/L
C
C     SROD1 computes the stress-displacement transformation.
C     The actual stress is: sigma = SAT . u_a + SBT . u_b (axial)
C                           tau   = SAR . u_a + SBR . u_b (torsion)
C
      PROGRAM TESTSROD1
      IMPLICIT NONE
C
C     =================================================================
C     COMMON blocks matching SROD1 expectations
C     =================================================================
C
C     SDR2X5 - Input/Output block for stress data recovery
C     SROD1 reads ECPT from here and writes stress coefficients
C
      REAL             ECPT(17), DUMMY1(83)
      INTEGER          IELID, ISILNO(2)
      REAL             SAT(3), SBT(3), SAR(3), SBR(3)
      REAL             ST, SDELTA, AREA, FJOVRC
      REAL             TSUBC0, SIGMAT, SIGMAC, SIGMAS
      REAL             SIGVEC(77), FORVEC(25)
      COMMON /SDR2X5/ ECPT, DUMMY1,
     1                IELID, ISILNO,
     2                SAT, SBT, SAR, SBR,
     3                ST, SDELTA, AREA, FJOVRC,
     4                TSUBC0, SIGMAT, SIGMAC, SIGMAS,
     5                SIGVEC, FORVEC
C
C     SDR2X6 - Scratch block
C
      REAL             XN(6), TI(9)
      REAL             XL, EOVERL
      INTEGER          IBASE
      COMMON /SDR2X6/ XN, TI, XL, EOVERL, IBASE
C
C     Material I/O
C
      INTEGER          MATIDC, MATFLG
      REAL             ELTEMP, STRESS, SINTH, COSTH
      COMMON /MATIN / MATIDC, MATFLG, ELTEMP, STRESS, SINTH, COSTH
C
      REAL             E, G, NU, RHO, ALPHA, TSUB0
      REAL             GSUBE, SIGT, SIGC, SIGS
      COMMON /MATOUT/ E, G, NU, RHO, ALPHA, TSUB0,
     1                GSUBE, SIGT, SIGC, SIGS
C
C     Local variables
C
      INTEGER          IECPT(17)
      EQUIVALENCE     (IECPT(1), ECPT(1))
C
      REAL             UA(3), UB(3), THETAA(3), THETAB(3)
      REAL             AXFORCE, AXSTRESS, TORQUE
      REAL             EXPECTED_STRESS, EXPECTED_TORQUE
      REAL             ROD_L, ROD_A, ROD_J, ROD_C
      INTEGER          I
C
C     =================================================================
C     Set up test parameters
C     =================================================================
C
      ROD_L = 2.0
      ROD_A = 0.01
      ROD_J = 5.0E-6
      ROD_C = 0.005
C
      WRITE(6,*) 'SROD1 Test Driver - Rod Stress Recovery'
      WRITE(6,*) '========================================'
      WRITE(6,*) ''
      WRITE(6,*) 'Rod parameters:'
      WRITE(6,*) '  Length L = 2.0 m'
      WRITE(6,*) '  Area   A = 0.01 m^2'
      WRITE(6,*) '  J        = 5.0e-6 m^4'
      WRITE(6,*) '  C (J/c)  = 0.005 m'
      WRITE(6,*) '  E        = 2.1e11 Pa'
      WRITE(6,*) '  NU       = 0.3'
      WRITE(6,*) ''
C
C     =================================================================
C     Fill ECPT for rod element
C     =================================================================
C     From srod1.f comments, the ECPT layout for ROD is:
C       IECPT(1)  = Element ID
C       IECPT(2)  = SIL node A
C       IECPT(3)  = SIL node B
C       IECPT(4)  = Material ID
C       ECPT(5)   = Cross-sectional area A
C       ECPT(6)   = Torsional constant J
C       ECPT(7)   = Torsional stress recovery coeff C (actually J/C)
C       ECPT(8)   = Non-structural mass per length
C       IECPT(9)  = Coord system ID for point A
C       ECPT(10)  = XA  (x-coordinate of point A in basic)
C       ECPT(11)  = YA
C       ECPT(12)  = ZA
C       IECPT(13) = Coord system ID for point B
C       ECPT(14)  = XB  (x-coordinate of point B in basic)
C       ECPT(15)  = YB
C       ECPT(16)  = ZB
C       ECPT(17)  = Element temperature
C
      DO 10 I = 1, 17
        ECPT(I) = 0.0
   10 CONTINUE
C
      IECPT(1) = 100
      IECPT(2) = 1
      IECPT(3) = 2
      IECPT(4) = 1
      ECPT(5)  = ROD_A
      ECPT(6)  = ROD_J
      ECPT(7)  = ROD_J / ROD_C
      ECPT(8)  = 0.0
      IECPT(9) = 0
      ECPT(10) = 0.0
      ECPT(11) = 0.0
      ECPT(12) = 0.0
      IECPT(13) = 0
      ECPT(14) = ROD_L
      ECPT(15) = 0.0
      ECPT(16) = 0.0
      ECPT(17) = 20.0
C
C     Zero the dummy1 area
      DO 15 I = 1, 83
        DUMMY1(I) = 0.0
   15 CONTINUE
C
C     Zero output arrays
      DO 16 I = 1, 77
        SIGVEC(I) = 0.0
   16 CONTINUE
      DO 17 I = 1, 25
        FORVEC(I) = 0.0
   17 CONTINUE
C
C     =================================================================
C     Set material properties (will also be set by MAT stub)
C     =================================================================
C
      E = 2.1E11
      NU = 0.3
      G = E / (2.0*(1.0+NU))
      RHO = 7800.0
      ALPHA = 12.0E-6
      TSUB0 = 20.0
      GSUBE = 0.0
      SIGT = 0.0
      SIGC = 0.0
      SIGS = 0.0
C
C     =================================================================
C     Call SROD1
C     =================================================================
C
      WRITE(6,*) 'Calling SROD1...'
      WRITE(6,*) ''
C
      CALL SROD1
C
C     =================================================================
C     Print SROD1 outputs (stress-displacement coefficients)
C     =================================================================
C
      WRITE(6,*) 'SROD1 Output - Stress Recovery Coefficients:'
      WRITE(6,*) '---------------------------------------------'
      WRITE(6,100) 'SAT', SAT(1), SAT(2), SAT(3)
      WRITE(6,100) 'SBT', SBT(1), SBT(2), SBT(3)
      WRITE(6,100) 'SAR', SAR(1), SAR(2), SAR(3)
      WRITE(6,100) 'SBR', SBR(1), SBR(2), SBR(3)
  100 FORMAT(2X,A3,' = [',E14.7,', ',E14.7,', ',E14.7,']')
C
      WRITE(6,*) ''
      WRITE(6,110) 'ST     ', ST
      WRITE(6,110) 'SDELTA ', SDELTA
      WRITE(6,110) 'AREA   ', AREA
      WRITE(6,110) 'FJOVRC ', FJOVRC
      WRITE(6,110) 'TSUBC0 ', TSUBC0
  110 FORMAT(2X,A7,' = ',E14.7)
C
      WRITE(6,*) ''
      WRITE(6,*) 'Element ID: ', IELID
      WRITE(6,*) 'SIL nodes:  ', ISILNO(1), ISILNO(2)
C
C     =================================================================
C     Apply coefficients to known displacements
C     =================================================================
C
      WRITE(6,*) ''
      WRITE(6,*) '================================================'
      WRITE(6,*) 'Applying to test displacements:'
      WRITE(6,*) ''
C
C     Test 1: Pure axial extension
C     u_a = (0, 0, 0), u_b = (0.001, 0, 0)
C
      UA(1) = 0.0
      UA(2) = 0.0
      UA(3) = 0.0
      UB(1) = 0.001
      UB(2) = 0.0
      UB(3) = 0.0
C
      WRITE(6,*) 'Test 1: Pure axial extension'
      WRITE(6,*) '  u_a = (0, 0, 0)'
      WRITE(6,*) '  u_b = (0.001, 0, 0)'
C
C     Axial stress = SAT.u_a + SBT.u_b
      AXSTRESS = SAT(1)*UA(1) + SAT(2)*UA(2) + SAT(3)*UA(3)
     1         + SBT(1)*UB(1) + SBT(2)*UB(2) + SBT(3)*UB(3)
      AXFORCE = AXSTRESS * ROD_A
C
C     Expected: stress = E * (u_b - u_a) / L = 2.1e11 * 0.001 / 2.0
      EXPECTED_STRESS = E * (UB(1) - UA(1)) / ROD_L
C
      WRITE(6,*) ''
      WRITE(6,200) AXSTRESS
      WRITE(6,210) EXPECTED_STRESS
      WRITE(6,220) AXFORCE
  200 FORMAT('  Computed axial stress  = ',E14.7,' Pa')
  210 FORMAT('  Expected axial stress  = ',E14.7,' Pa')
  220 FORMAT('  Axial force            = ',E14.7,' N')
C
      IF (ABS(AXSTRESS - EXPECTED_STRESS) .LT.
     1    ABS(EXPECTED_STRESS)*1.0E-5) THEN
        WRITE(6,*) '  >> PASS: Axial stress matches.'
      ELSE
        WRITE(6,*) '  >> FAIL: Axial stress mismatch!'
      END IF
C
C     =================================================================
C     Test 2: Pure torsion
C     theta_a = (0.01, 0, 0), theta_b = (0.02, 0, 0)
C     For rod along x, torsion is rotation about x-axis
C     Torque = G*J*(theta_b - theta_a)/L
C
      WRITE(6,*) ''
      WRITE(6,*) 'Test 2: Pure torsion (rod along x-axis)'
      WRITE(6,*) '  theta_a = (0.01, 0, 0)  [rot about x at A]'
      WRITE(6,*) '  theta_b = (0.02, 0, 0)  [rot about x at B]'
C
      THETAA(1) = 0.01
      THETAA(2) = 0.0
      THETAA(3) = 0.0
      THETAB(1) = 0.02
      THETAB(2) = 0.0
      THETAB(3) = 0.0
C
C     Torsional stress = SAR.theta_a + SBR.theta_b
      TORQUE = SAR(1)*THETAA(1) + SAR(2)*THETAA(2) + SAR(3)*THETAA(3)
     1       + SBR(1)*THETAB(1) + SBR(2)*THETAB(2) + SBR(3)*THETAB(3)
C
C     Expected torsional stress:
C     The rod formulation gives torsional stress as G*C*(theta_b-theta_a)/L
C     where C is the torsional stress coefficient
C     Actually SAR/SBR are stress-per-displacement, and for torsion:
C       torsional_stress = G*J/(C*L) * (delta_theta along rod axis)
C     SROD1: GCOVRL = G * J / L, SAR = XN * GCOVRL
C     So: torsion_result = SAR.theta_a + SBR.theta_b
C     For rod along x: SAR = (G*J/L, 0, 0), SBR = (-G*J/L, 0, 0)
C     result = G*J/L * theta_a_x + (-G*J/L) * theta_b_x
C            = G*J/L * (theta_a_x - theta_b_x)
C     The sign convention: SROD1 computes from A to B direction
C     Actually looking at srod1.f: XN = (XA-XB) direction
C     So for XA=(0,0,0), XB=(2,0,0): XN = (-1,0,0) after normalization
C     SAR = XN * GCOVRL = (-G*J/L, 0, 0)
C     SBR = -XN * GCOVRL = (G*J/L, 0, 0)
C     result = -G*J/L * theta_a_x + G*J/L * theta_b_x
C            = G*J/L * (theta_b_x - theta_a_x)
C
      EXPECTED_TORQUE = G * ROD_J / ROD_L * (THETAB(1) - THETAA(1))
C
      WRITE(6,*) ''
      WRITE(6,300) TORQUE
      WRITE(6,310) EXPECTED_TORQUE
  300 FORMAT('  Computed torsion result = ',E14.7)
  310 FORMAT('  Expected G*J/L*dtheta  = ',E14.7)
C
      IF (ABS(TORQUE - EXPECTED_TORQUE) .LT.
     1    ABS(EXPECTED_TORQUE)*1.0E-5) THEN
        WRITE(6,*) '  >> PASS: Torsion matches.'
      ELSE
        WRITE(6,*) '  >> FAIL: Torsion mismatch!'
        WRITE(6,*) '  (Note: sign depends on SROD1 XN direction)'
      END IF
C
C     =================================================================
C     Summary
C     =================================================================
C
      WRITE(6,*) ''
      WRITE(6,*) '================================================'
      WRITE(6,*) 'Additional verification:'
      WRITE(6,*) '  E/L (EOVERL) should be:', E/ROD_L
      WRITE(6,400) EOVERL
  400 FORMAT('  SROD1 computed EOVERL = ',E14.7)
      WRITE(6,*) '  Rod length should be:', ROD_L
      WRITE(6,410) XL
  410 FORMAT('  SROD1 computed XL    = ',F10.6)
      WRITE(6,*) ''
      WRITE(6,*) 'Computation completed successfully.'
C
      STOP
      END
