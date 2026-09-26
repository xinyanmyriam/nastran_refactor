C     TEST DRIVER FOR TRD1C - Newmark-Beta Time Integration (B11)
C     =============================================================
C     Standalone reference implementation of Newmark-beta method
C     for transient dynamic response.
C
C     trd1c.f in NASTRAN is deeply embedded in the transient response
C     module with many COMMON blocks and file I/O. This standalone
C     implements the same Newmark-beta algorithm.
C
C     Test system 1: Single DOF damped forced vibration
C       m*u'' + c*u' + k*u = F(t)
C       m = 1.0, c = 0.1, k = 100.0
C       F(t) = sin(10*t)
C       u(0) = 0, u'(0) = 0
C       beta = 0.25 (average acceleration), gamma = 0.5
C       dt = 0.01, 100 time steps
C
C     Test system 2: 2-DOF spring-mass chain
C       Verifies matrix assembly in Newmark integration
C
      PROGRAM TESTTRD1C
      IMPLICIT NONE
C
      INTEGER NSTEPS
      PARAMETER (NSTEPS = 100)
C
      DOUBLE PRECISION MASS, DAMP, STIFF
      DOUBLE PRECISION BETA, GAMMA, DT
      DOUBLE PRECISION U, V, ACC, UNEW, VNEW, ANEW
      DOUBLE PRECISION FORCE, FNEW, TIME
      DOUBLE PRECISION OMEGA, OMEGA_D, ZETA, OMEGA_N
      DOUBLE PRECISION A0, A1, A2, A3, A4, A5, KEFF, FEFF
      DOUBLE PRECISION PI, ENERGY, KE, PE
      INTEGER I
C
C     2-DOF system variables
      INTEGER N2
      PARAMETER (N2 = 2)
      DOUBLE PRECISION M2(N2,N2), C2B(N2,N2), K2B(N2,N2)
      DOUBLE PRECISION U2(N2), V2(N2), A2(N2)
      DOUBLE PRECISION U2N(N2), V2N(N2), A2N(N2)
      DOUBLE PRECISION F2(N2), KEFF2(N2,N2), FEFF2(N2)
      DOUBLE PRECISION L2(N2,N2), TEMP2(N2), X2(N2)
      DOUBLE PRECISION SUM1
      INTEGER J, KK
C
      PI = 4.0D0*DATAN(1.0D0)
C
C     System parameters
      MASS  = 1.0D0
      DAMP  = 0.1D0
      STIFF = 100.0D0
      BETA  = 0.25D0
      GAMMA = 0.5D0
      DT    = 0.01D0
C
      WRITE(6,*) 'TRD1C Test Driver - Newmark-Beta Time Integration'
      WRITE(6,*) '=================================================='
      WRITE(6,*) ''
      WRITE(6,*) '=== Test 1: Single DOF System ==='
      WRITE(6,*) '  m = 1.0, c = 0.1, k = 100.0'
      WRITE(6,*) '  F(t) = sin(10t)'
      WRITE(6,*) '  beta = 0.25, gamma = 0.5, dt = 0.01'
      WRITE(6,*) ''
C
C     Derived quantities
      OMEGA_N = DSQRT(STIFF/MASS)
      ZETA = DAMP / (2.0D0*MASS*OMEGA_N)
      OMEGA_D = OMEGA_N * DSQRT(1.0D0 - ZETA**2)
      OMEGA = 10.0D0
C
      WRITE(6,200) OMEGA_N, ZETA, OMEGA_D
  200 FORMAT('  omega_n = ',F8.4,', zeta = ',F8.6,
     1       ', omega_d = ',F8.4)
      WRITE(6,*) ''
C
C     Newmark integration constants
      A0 = 1.0D0/(BETA*DT**2)
      A1 = GAMMA/(BETA*DT)
      A2 = 1.0D0/(BETA*DT)
      A3 = 1.0D0/(2.0D0*BETA) - 1.0D0
      A4 = GAMMA/BETA - 1.0D0
      A5 = DT*(GAMMA/(2.0D0*BETA) - 1.0D0)
C
C     Effective stiffness
      KEFF = STIFF + A0*MASS + A1*DAMP
C
C     Initial conditions
      U = 0.0D0
      V = 0.0D0
      FORCE = 0.0D0
C     Initial acceleration: a(0) = (F(0) - c*v(0) - k*u(0))/m
      ACC = (FORCE - DAMP*V - STIFF*U) / MASS
C
C     Time integration loop
      WRITE(6,*) '  Step    Time      Displacement    Velocity'
     1          //'        Acceleration'
      WRITE(6,*) '  ----    ----      ------------    --------'
     1          //'        ------------'
C
      DO 100 I = 1, NSTEPS
        TIME = DBLE(I) * DT
C
C       Force at new time
        FNEW = DSIN(OMEGA*TIME)
C
C       Effective force
        FEFF = FNEW + MASS*(A0*U + A2*V + A3*ACC)
     1              + DAMP*(A1*U + A4*V + A5*ACC)
C
C       Solve for new displacement
        UNEW = FEFF / KEFF
C
C       New acceleration and velocity
        ANEW = A0*(UNEW - U) - A2*V - A3*ACC
        VNEW = V + DT*((1.0D0-GAMMA)*ACC + GAMMA*ANEW)
C
C       Print every 10 steps
        IF (MOD(I,10) .EQ. 0) THEN
          WRITE(6,300) I, TIME, UNEW, VNEW, ANEW
        END IF
  300   FORMAT(2X,I4,4X,F6.3,4X,E14.7,4X,E14.7,4X,E14.7)
C
C       Update for next step
        U = UNEW
        V = VNEW
        ACC = ANEW
C
  100 CONTINUE
C
C     Energy balance check
      KE = 0.5D0*MASS*V**2
      PE = 0.5D0*STIFF*U**2
      ENERGY = KE + PE
C
      WRITE(6,*) ''
      WRITE(6,400) U
  400 FORMAT('  Final displacement at t=1.0: u = ',E16.9)
      WRITE(6,410) V
  410 FORMAT('  Final velocity at t=1.0:     v = ',E16.9)
      WRITE(6,420) KE, PE, ENERGY
  420 FORMAT('  Energy: KE=',E12.5,', PE=',E12.5,
     1       ', Total=',E12.5)
C
C     =========================================================
C     Test 2: 2-DOF Newmark (matrix operations verification)
C     =========================================================
C
      WRITE(6,*) ''
      WRITE(6,*) '=== Test 2: 2-DOF Spring-Mass Chain ==='
      WRITE(6,*) '  m1=m2=1, k1=k2=100, c1=c2=0.2'
      WRITE(6,*) '  F1(t)=sin(10t), F2(t)=0'
      WRITE(6,*) ''
C
C     Mass matrix (lumped)
      M2(1,1) = 1.0D0
      M2(1,2) = 0.0D0
      M2(2,1) = 0.0D0
      M2(2,2) = 1.0D0
C
C     Stiffness matrix (2-spring chain, fixed-free)
      K2B(1,1) = 200.0D0
      K2B(1,2) = -100.0D0
      K2B(2,1) = -100.0D0
      K2B(2,2) = 100.0D0
C
C     Damping matrix (Rayleigh: C = 0.001*K)
      C2B(1,1) = 0.2D0
      C2B(1,2) = -0.1D0
      C2B(2,1) = -0.1D0
      C2B(2,2) = 0.1D0
C
C     Effective stiffness: Keff = K + a0*M + a1*C
      DO 510 I = 1, N2
        DO 510 J = 1, N2
          KEFF2(I,J) = K2B(I,J) + A0*M2(I,J) + A1*C2B(I,J)
  510 CONTINUE
C
C     Factor Keff (2x2 Cholesky)
      L2(1,1) = DSQRT(KEFF2(1,1))
      L2(1,2) = 0.0D0
      L2(2,1) = KEFF2(2,1) / L2(1,1)
      L2(2,2) = DSQRT(KEFF2(2,2) - L2(2,1)**2)
C
C     Initial conditions
      DO 520 I = 1, N2
        U2(I) = 0.0D0
        V2(I) = 0.0D0
        A2(I) = 0.0D0
  520 CONTINUE
C
C     Time integration for 2-DOF
      WRITE(6,*) '  Step    Time       U1             U2'
      WRITE(6,*) '  ----    ----       --             --'
C
      DO 600 I = 1, NSTEPS
        TIME = DBLE(I) * DT
C
C       External force
        F2(1) = DSIN(OMEGA*TIME)
        F2(2) = 0.0D0
C
C       Effective force: F + M*(a0*u + a2*v + a3*a) + C*(a1*u + a4*v + a5*a)
        DO 530 J = 1, N2
          FEFF2(J) = F2(J)
          DO 530 KK = 1, N2
            FEFF2(J) = FEFF2(J)
     1        + M2(J,KK)*(A0*U2(KK) + A2*V2(KK) + A3*A2(KK))
     2        + C2B(J,KK)*(A1*U2(KK) + A4*V2(KK) + A5*A2(KK))
  530   CONTINUE
C
C       Solve Keff * u_new = Feff (forward/backward sub)
C       Forward: L*y = Feff
        TEMP2(1) = FEFF2(1) / L2(1,1)
        TEMP2(2) = (FEFF2(2) - L2(2,1)*TEMP2(1)) / L2(2,2)
C       Backward: L^T*x = y
        X2(2) = TEMP2(2) / L2(2,2)
        X2(1) = (TEMP2(1) - L2(2,1)*X2(2)) / L2(1,1)
C
        U2N(1) = X2(1)
        U2N(2) = X2(2)
C
C       New accelerations and velocities
        DO 540 J = 1, N2
          A2N(J) = A0*(U2N(J) - U2(J)) - A2*V2(J) - A3*A2(J)
          V2N(J) = V2(J) + DT*((1.0D0-GAMMA)*A2(J) + GAMMA*A2N(J))
  540   CONTINUE
C
C       Print every 10 steps
        IF (MOD(I,10) .EQ. 0) THEN
          WRITE(6,610) I, TIME, U2N(1), U2N(2)
        END IF
  610   FORMAT(2X,I4,4X,F6.3,4X,E14.7,4X,E14.7)
C
C       Update
        DO 550 J = 1, N2
          U2(J) = U2N(J)
          V2(J) = V2N(J)
          A2(J) = A2N(J)
  550   CONTINUE
C
  600 CONTINUE
C
      WRITE(6,*) ''
      WRITE(6,700) U2(1), U2(2)
  700 FORMAT('  Final: U1=',E14.7,', U2=',E14.7)
      WRITE(6,710) V2(1), V2(2)
  710 FORMAT('         V1=',E14.7,', V2=',E14.7)
C
      WRITE(6,*) ''
      WRITE(6,*) 'Computation completed successfully.'
      STOP
      END
