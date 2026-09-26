C======================================================================
C     TEST_NEWMARK.F - Standalone Newmark-Beta Time Integration
C
C     This program implements the Newmark-Beta (average acceleration)
C     algorithm extracted from NASTRAN-95's TRD1C.F/STEP.F subroutines.
C     Instead of GINO tape I/O, it uses explicit arrays.
C
C     Uses NASTRAN's INVERD subroutine for the linear system solve.
C
C     Problem: 3-DOF spring-mass system (both ends fixed), free vibration
C       M = I (identity), C = 0 (no damping)
C       K = [2 -1 0; -1 2 -1; 0 -1 2]  (standard tridiagonal)
C       u0 = [1, 0, 0], v0 = [0, 0, 0]
C       dt = 0.01, 100 steps
C
C     Eigenvalues: lambda_k = 2 - 2*cos(k*pi/4), k=1,2,3
C       lambda_1 = 2 - sqrt(2) = 0.5858
C       lambda_2 = 2.0
C       lambda_3 = 2 + sqrt(2) = 3.4142
C
C     Eigenvectors (M-orthonormal):
C       phi_k(i) = sqrt(2/4) * sin(i*k*pi/4)
C
C     Newmark parameters: beta=1/4, gamma=1/2 (average acceleration)
C     This is the same scheme used by NASTRAN-95's TRD1C subroutine.
C======================================================================
      PROGRAM TEST_NEWMARK
      IMPLICIT NONE
C
C     Parameters
      INTEGER NDOF, NSTEP
      PARAMETER (NDOF=3, NSTEP=100)
      DOUBLE PRECISION PI
      PARAMETER (PI=3.141592653589793D0)
C
C     Arrays
      DOUBLE PRECISION K(NDOF,NDOF), M(NDOF,NDOF), KEFF(NDOF,NDOF)
      DOUBLE PRECISION U(NDOF), V(NDOF), A(NDOF)
      DOUBLE PRECISION UNEW(NDOF), VNEW(NDOF), ANEW(NDOF)
      DOUBLE PRECISION RHS(NDOF,1)
      DOUBLE PRECISION KEFF_SAVE(NDOF,NDOF)
      INTEGER INDEX(NDOF,3)
C
C     Newmark parameters
      DOUBLE PRECISION BETA, GAMMA, DT
      DOUBLE PRECISION A0, A2, A3, A6, A7
C
C     Analytical solution arrays
      DOUBLE PRECISION OMEGA(NDOF), PHI(NDOF,NDOF)
      DOUBLE PRECISION UANAL(NDOF), MODAL_AMP(NDOF)
C
C     Energy check
      DOUBLE PRECISION KE, PE, TE, TE0
C
C     Working variables
      DOUBLE PRECISION DETERM, T, ERR, MAXERR, RMSERR
      DOUBLE PRECISION TEMP(NDOF)
      INTEGER I, J, KK, ISTEP, ISING, IOUT
      INTEGER MACH_VAL
C
C     Common block required by INVERD
      COMMON /MACHIN/ MACH_VAL
C
C     Initialize machine type (1=generic)
      MACH_VAL = 1
C
C     Set Newmark parameters
      BETA  = 0.25D0
      GAMMA = 0.50D0
      DT    = 0.01D0
C
C     Integration constants (Newmark-beta average acceleration)
      A0 = 1.0D0 / (BETA * DT * DT)
      A2 = 1.0D0 / (BETA * DT)
      A3 = 1.0D0 / (2.0D0 * BETA) - 1.0D0
      A6 = DT * (1.0D0 - GAMMA)
      A7 = GAMMA * DT
C
C     ================================================================
C     Define the 3-DOF system (both ends fixed)
C     ================================================================
C
C     Mass matrix M = I
      DO I = 1, NDOF
        DO J = 1, NDOF
          M(I,J) = 0.0D0
        END DO
        M(I,I) = 1.0D0
      END DO
C
C     Stiffness matrix K = [2 -1 0; -1 2 -1; 0 -1 2]
      DO I = 1, NDOF
        DO J = 1, NDOF
          K(I,J) = 0.0D0
        END DO
      END DO
      K(1,1) =  2.0D0
      K(1,2) = -1.0D0
      K(2,1) = -1.0D0
      K(2,2) =  2.0D0
      K(2,3) = -1.0D0
      K(3,2) = -1.0D0
      K(3,3) =  2.0D0
C
C     ================================================================
C     Form effective stiffness: K_eff = K + a0*M (C=0)
C     ================================================================
      DO I = 1, NDOF
        DO J = 1, NDOF
          KEFF(I,J) = K(I,J) + A0 * M(I,J)
        END DO
      END DO
C
C     Save a copy of K_eff (INVERD modifies the matrix in-place)
      DO I = 1, NDOF
        DO J = 1, NDOF
          KEFF_SAVE(I,J) = KEFF(I,J)
        END DO
      END DO
C
C     ================================================================
C     Initial conditions
C     ================================================================
      U(1) = 1.0D0
      U(2) = 0.0D0
      U(3) = 0.0D0
      V(1) = 0.0D0
      V(2) = 0.0D0
      V(3) = 0.0D0
C
C     Initial acceleration: M*a0 = F0 - K*u0 - C*v0
C     F0 = 0 (free vibration), C=0, M=I, so a0 = -K*u0
      DO I = 1, NDOF
        A(I) = 0.0D0
        DO J = 1, NDOF
          A(I) = A(I) - K(I,J) * U(J)
        END DO
      END DO
C
C     ================================================================
C     Set up analytical solution (modal superposition)
C     ================================================================
C
C     For tridiagonal K with fixed ends (N=3, N+1=4 segments):
C     Eigenvalues: lambda_k = 2 - 2*cos(k*pi/(N+1))
C     Eigenvectors: phi_k(i) = sqrt(2/(N+1)) * sin(i*k*pi/(N+1))
C
      DO KK = 1, NDOF
        OMEGA(KK) = DSQRT(2.0D0 - 2.0D0*DCOS(DBLE(KK)*PI/4.0D0))
        DO I = 1, NDOF
          PHI(I,KK) = DSQRT(2.0D0/4.0D0)
     &              * DSIN(DBLE(I)*DBLE(KK)*PI/4.0D0)
        END DO
      END DO
C
C     Modal amplitudes: q_k(0) = phi_k^T * M * u0 = phi_k^T * u0
      DO KK = 1, NDOF
        MODAL_AMP(KK) = 0.0D0
        DO J = 1, NDOF
          MODAL_AMP(KK) = MODAL_AMP(KK) + PHI(J,KK) * U(J)
        END DO
      END DO
C
C     ================================================================
C     Compute initial total energy (KE + PE)
C     ================================================================
      KE = 0.0D0
      DO I = 1, NDOF
        KE = KE + 0.5D0 * V(I) * V(I)
      END DO
      PE = 0.0D0
      DO I = 1, NDOF
        DO J = 1, NDOF
          PE = PE + 0.5D0 * U(I) * K(I,J) * U(J)
        END DO
      END DO
      TE0 = KE + PE
C
C     ================================================================
C     Print header
C     ================================================================
      WRITE(*,900)
  900 FORMAT('=============================================='/
     &       ' NASTRAN-95 Newmark-Beta Standalone Test'/
     &       ' (Algorithm from TRD1C.F + STEP.F)'/
     &       ' Linear solve via INVERD.F'/
     &       '=============================================='/)
      WRITE(*,901) NDOF, NSTEP, DT, BETA, GAMMA
  901 FORMAT(' Problem: ',I1,'-DOF spring-mass (fixed-fixed),',
     &       ' free vibration'/
     &       ' Steps:   ',I4/
     &       ' dt:      ',F8.5/
     &       ' beta:    ',F6.4/
     &       ' gamma:   ',F6.4/)
      WRITE(*,902) OMEGA(1), OMEGA(2), OMEGA(3)
  902 FORMAT(' Natural frequencies:'/
     &       '   omega_1 = ',F10.6,' rad/s'/
     &       '   omega_2 = ',F10.6,' rad/s'/
     &       '   omega_3 = ',F10.6,' rad/s'/)
      WRITE(*,903) (MODAL_AMP(I), I=1,NDOF)
  903 FORMAT(' Modal participation factors:'/
     &       '   q1(0) = ',F10.6/
     &       '   q2(0) = ',F10.6/
     &       '   q3(0) = ',F10.6/)
      WRITE(*,904) TE0
  904 FORMAT(' Initial total energy: ',F12.6/)
      WRITE(*,905)
  905 FORMAT(' Step    Time      u1_num      u1_anal  ',
     &       '   u2_num      u2_anal     u3_num   ',
     &       '   u3_anal     Error       Energy'/
     &       ' ----  ------  ----------  ----------',
     &       '  ----------  ----------  ----------',
     &       '  ----------  ----------  ----------')
C
C     ================================================================
C     Print initial state
C     ================================================================
      T = 0.0D0
      DO I = 1, NDOF
        UANAL(I) = 0.0D0
        DO KK = 1, NDOF
          UANAL(I) = UANAL(I)
     &             + PHI(I,KK)*MODAL_AMP(KK)*DCOS(OMEGA(KK)*T)
        END DO
      END DO
      ERR = 0.0D0
      DO I = 1, NDOF
        ERR = ERR + (U(I) - UANAL(I))**2
      END DO
      ERR = DSQRT(ERR)
      WRITE(*,910) 0, T, U(1), UANAL(1), U(2), UANAL(2),
     &             U(3), UANAL(3), ERR, TE0
  910 FORMAT(I5,F8.4,3(F12.7,F12.7),2E12.4)
C
C     ================================================================
C     Time integration loop (Newmark-beta average acceleration)
C     ================================================================
      MAXERR = 0.0D0
      RMSERR = 0.0D0
C
      DO 500 ISTEP = 1, NSTEP
        T = DBLE(ISTEP) * DT
C
C       Form RHS = F(t+dt) + M*(a0*u + a2*v + a3*a)
C       F(t+dt) = 0 (free vibration), C=0, so:
C       RHS = M * (a0*u + a2*v + a3*a) = a0*u + a2*v + a3*a (M=I)
C
        DO I = 1, NDOF
          RHS(I,1) = A0*U(I) + A2*V(I) + A3*A(I)
        END DO
C
C       Solve K_eff * u_{n+1} = RHS using INVERD
C       Restore K_eff (INVERD destroys it)
C
        DO I = 1, NDOF
          DO J = 1, NDOF
            KEFF(I,J) = KEFF_SAVE(I,J)
          END DO
        END DO
C
        ISING = -1
        CALL INVERD(NDOF, KEFF, NDOF, RHS, 1, DETERM, ISING, INDEX)
C
        IF (ISING .EQ. 2) THEN
          WRITE(*,*) 'ERROR: Singular matrix at step ', ISTEP
          STOP
        END IF
C
C       RHS now contains u_{n+1}
        DO I = 1, NDOF
          UNEW(I) = RHS(I,1)
        END DO
C
C       Update acceleration: a_{n+1} = a0*(u_{n+1}-u_n) - a2*v_n - a3*a_n
        DO I = 1, NDOF
          ANEW(I) = A0*(UNEW(I) - U(I)) - A2*V(I) - A3*A(I)
        END DO
C
C       Update velocity: v_{n+1} = v_n + a6*a_n + a7*a_{n+1}
        DO I = 1, NDOF
          VNEW(I) = V(I) + A6*A(I) + A7*ANEW(I)
        END DO
C
C       Shift state
        DO I = 1, NDOF
          U(I) = UNEW(I)
          V(I) = VNEW(I)
          A(I) = ANEW(I)
        END DO
C
C       Compute analytical solution at this time
        DO I = 1, NDOF
          UANAL(I) = 0.0D0
          DO KK = 1, NDOF
            UANAL(I) = UANAL(I)
     &               + PHI(I,KK)*MODAL_AMP(KK)*DCOS(OMEGA(KK)*T)
          END DO
        END DO
C
C       Compute error (L2 norm)
        ERR = 0.0D0
        DO I = 1, NDOF
          ERR = ERR + (U(I) - UANAL(I))**2
        END DO
        ERR = DSQRT(ERR)
        IF (ERR .GT. MAXERR) MAXERR = ERR
        RMSERR = RMSERR + ERR*ERR
C
C       Compute total energy
        KE = 0.0D0
        DO I = 1, NDOF
          KE = KE + 0.5D0 * V(I) * V(I)
        END DO
        PE = 0.0D0
        DO I = 1, NDOF
          DO J = 1, NDOF
            PE = PE + 0.5D0 * U(I) * K(I,J) * U(J)
          END DO
        END DO
        TE = KE + PE
C
C       Output at selected steps
        IOUT = 0
        IF (ISTEP .LE. 5) IOUT = 1
        IF (MOD(ISTEP,10) .EQ. 0) IOUT = 1
        IF (ISTEP .EQ. NSTEP) IOUT = 1
C
        IF (IOUT .EQ. 1) THEN
          WRITE(*,910) ISTEP, T, U(1), UANAL(1), U(2), UANAL(2),
     &                 U(3), UANAL(3), ERR, TE
        END IF
C
  500 CONTINUE
C
C     ================================================================
C     Summary
C     ================================================================
      RMSERR = DSQRT(RMSERR / DBLE(NSTEP))
      WRITE(*,920) MAXERR, RMSERR, TE0, TE, DABS(TE-TE0)/TE0
  920 FORMAT(//' ================================================'/
     &       ' Summary Results'/
     &       ' ================================================'/
     &       '   Max L2 error vs analytical: ',E12.5/
     &       '   RMS L2 error vs analytical: ',E12.5/
     &       '   Initial energy:             ',F12.8/
     &       '   Final energy:               ',F12.8/
     &       '   Energy drift (relative):    ',E12.5/)
C
      IF (MAXERR .LT. 1.0D-4) THEN
        WRITE(*,*) 'RESULT: PASS - Newmark integration matches'
        WRITE(*,*) '  analytical solution within 1e-4'
      ELSE IF (MAXERR .LT. 1.0D-2) THEN
        WRITE(*,*) 'RESULT: ACCEPTABLE - Max error ',MAXERR
        WRITE(*,*) '  (within 1e-2; reduce dt for better accuracy)'
      ELSE
        WRITE(*,*) 'RESULT: FAIL - Error exceeds 1e-2'
      END IF
C
      IF (DABS(TE-TE0)/TE0 .LT. 1.0D-10) THEN
        WRITE(*,*) 'ENERGY: CONSERVED (drift < 1e-10)'
      ELSE
        WRITE(*,*) 'ENERGY: drift = ', DABS(TE-TE0)/TE0
      END IF
C
      WRITE(*,930)
  930 FORMAT(/
     &  ' Algorithm correspondence to NASTRAN-95 TRD1C:'/
     &  '   - Newmark-beta (beta=1/4, gamma=1/2) = average accel.'/
     &  '   - K_eff = K + (4/dt^2)*M  [matches SCR1 formation]'/
     &  '   - At each step: solve K_eff*u_{n+1} = effective load'/
     &  '   - This is equivalent to STEP.F: MATVEC+MATVEC+INTFBS'/
     &  '   - Linear solve via INVERD (Gauss-Jordan w/ pivoting)'/
     &  '   - INVERD is the actual NASTRAN-95 equation solver'/)
C
      STOP
      END
