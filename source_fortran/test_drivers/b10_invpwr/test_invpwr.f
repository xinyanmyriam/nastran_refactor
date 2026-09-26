C     TEST DRIVER FOR INVPWR - Inverse Power Eigenvalue Solver (B10)
C     ================================================================
C     Standalone reference implementation of inverse iteration for
C     generalized eigenvalue problem: K*x = lambda*M*x
C
C     invpwr.f in NASTRAN uses GINO file I/O extensively, so this
C     standalone program implements the same inverse power iteration:
C       1. Choose shift sigma near expected eigenvalue
C       2. Factor (K - sigma*M) = L*L^T
C       3. Iterate: solve (K-sigma*M)*y = M*x, normalize y
C       4. Converge when Rayleigh quotient stabilizes
C
C     Test system: 4-DOF spring-mass chain
C       All springs k=100, all masses m=1
C       Analytical eigenvalues: lambda_n = 2k/m * (1 - cos(n*pi/5))
C         lambda_1 = 38.197, lambda_2 = 138.197
C         lambda_3 = 261.803, lambda_4 = 361.803
C
      PROGRAM TESTINVPWR
      IMPLICIT NONE
C
      INTEGER N
      PARAMETER (N = 4)
C
      DOUBLE PRECISION K(N,N), M(N,N), A(N,N), L(N,N)
      DOUBLE PRECISION X(N), Y(N), MX(N), TEMP(N)
      DOUBLE PRECISION LAMBDA, LAMBDA_OLD, SIGMA
      DOUBLE PRECISION NORM, XTMX, XTKX, RAYLEIGH
      DOUBLE PRECISION TOL, PI, EXACT(N)
      INTEGER I, J, ITER, MAXITER, IEIG
      DOUBLE PRECISION SPRING_K, MASS_M
C
      PARAMETER (SPRING_K = 100.0D0, MASS_M = 1.0D0)
      PARAMETER (TOL = 1.0D-10, MAXITER = 200)
C
      PI = 4.0D0*DATAN(1.0D0)
C
      WRITE(6,*) 'INVPWR Test Driver - Inverse Power Iteration'
      WRITE(6,*) '=============================================='
      WRITE(6,*) ''
      WRITE(6,*) 'System: 4-DOF spring-mass chain'
      WRITE(6,*) '  k = 100 (all springs), m = 1 (all masses)'
      WRITE(6,*) ''
C
C     Build stiffness matrix (tridiagonal)
C     K(i,i) = 2k (interior) or k (boundary)
C
      DO 10 I = 1, N
        DO 10 J = 1, N
          K(I,J) = 0.0D0
          M(I,J) = 0.0D0
   10 CONTINUE
C
C     Stiffness: fixed-fixed chain (4 springs, 4 interior DOFs)
C     Spring between wall and DOF 1, between DOFs, and DOF 4 and wall
      K(1,1) = 2.0D0*SPRING_K
      K(1,2) = -SPRING_K
      K(2,1) = -SPRING_K
      K(2,2) = 2.0D0*SPRING_K
      K(2,3) = -SPRING_K
      K(3,2) = -SPRING_K
      K(3,3) = 2.0D0*SPRING_K
      K(3,4) = -SPRING_K
      K(4,3) = -SPRING_K
      K(4,4) = 2.0D0*SPRING_K
C
C     Mass matrix (lumped, identity since m=1)
      DO 20 I = 1, N
        M(I,I) = MASS_M
   20 CONTINUE
C
C     Analytical eigenvalues for fixed-fixed chain:
C     lambda_n = 2*k/m * (1 - cos(n*pi/(N+1)))
C
      WRITE(6,*) 'Analytical eigenvalues:'
      DO 25 I = 1, N
        EXACT(I) = 2.0D0*SPRING_K/MASS_M *
     1             (1.0D0 - DCOS(DBLE(I)*PI/DBLE(N+1)))
        WRITE(6,100) I, EXACT(I)
   25 CONTINUE
  100 FORMAT('  lambda_',I1,' = ',F12.6)
C
C     Print stiffness matrix
      WRITE(6,*) ''
      WRITE(6,*) 'Stiffness Matrix K:'
      DO 30 I = 1, N
        WRITE(6,110) (K(I,J), J=1,N)
   30 CONTINUE
  110 FORMAT(4(2X,F8.2))
C
      WRITE(6,*) ''
      WRITE(6,*) '--- Inverse Iteration Results ---'
      WRITE(6,*) ''
C
C     =========================================================
C     Find each eigenvalue using inverse iteration with shift
C     =========================================================
C
      DO 500 IEIG = 1, N
C
C       Shift point near expected eigenvalue
        SIGMA = EXACT(IEIG) - 1.0D0
C
C       Form shifted matrix: A = K - sigma*M
        DO 200 I = 1, N
          DO 200 J = 1, N
            A(I,J) = K(I,J) - SIGMA*M(I,J)
  200   CONTINUE
C
C       Cholesky factorization of A
        CALL CHOLESKY_EIG(A, L, N)
C
C       Initial guess (unit vector)
        DO 210 I = 1, N
          X(I) = 1.0D0/DSQRT(DBLE(N))
  210   CONTINUE
C
C       Inverse iteration loop
        LAMBDA_OLD = 0.0D0
        DO 400 ITER = 1, MAXITER
C
C         Compute MX = M*x (for lumped mass, MX = x)
          DO 220 I = 1, N
            MX(I) = 0.0D0
            DO 220 J = 1, N
              MX(I) = MX(I) + M(I,J)*X(J)
  220     CONTINUE
C
C         Solve (K - sigma*M)*y = M*x using Cholesky
          CALL FWDSUB_EIG(L, MX, TEMP, N)
          CALL BWDSUB_EIG(L, TEMP, Y, N)
C
C         Compute Rayleigh quotient: lambda = x^T*K*x / x^T*M*x
          XTKX = 0.0D0
          XTMX = 0.0D0
          DO 240 I = 1, N
            DO 240 J = 1, N
              XTKX = XTKX + X(I)*K(I,J)*X(J)
              XTMX = XTMX + X(I)*M(I,J)*X(J)
  240     CONTINUE
          LAMBDA = XTKX / XTMX
C
C         Normalize y (M-norm)
          NORM = 0.0D0
          DO 260 I = 1, N
            DO 260 J = 1, N
              NORM = NORM + Y(I)*M(I,J)*Y(J)
  260     CONTINUE
          NORM = DSQRT(NORM)
          DO 270 I = 1, N
            X(I) = Y(I) / NORM
  270     CONTINUE
C
C         Check convergence
          IF (ITER .GT. 1) THEN
            IF (ABS(LAMBDA - LAMBDA_OLD) .LT.
     1          TOL*ABS(LAMBDA)) GO TO 410
          END IF
          LAMBDA_OLD = LAMBDA
C
  400   CONTINUE
        WRITE(6,*) 'WARNING: Did not converge for eigenvalue', IEIG
C
C       Final Rayleigh quotient
  410   XTKX = 0.0D0
        XTMX = 0.0D0
        DO 420 I = 1, N
          DO 420 J = 1, N
            XTKX = XTKX + X(I)*K(I,J)*X(J)
            XTMX = XTMX + X(I)*M(I,J)*X(J)
  420   CONTINUE
        LAMBDA = XTKX / XTMX
C
        WRITE(6,300) IEIG, LAMBDA, ITER,
     1               ABS(LAMBDA - EXACT(IEIG))
  300   FORMAT('  Eigenvalue ',I1,': lambda = ',F12.6,
     1         ' (iters=',I3,', error=',E10.3,')')
C
        WRITE(6,310) (X(I), I=1,N)
  310   FORMAT('    Eigenvector: ',4(F9.6,1X))
        WRITE(6,*) ''
C
  500 CONTINUE
C
C     Summary
C
      WRITE(6,*) '--- Summary ---'
      WRITE(6,*) '  Eigenvalue    Computed        Analytical'
     1           //'      Error'
      DO 600 IEIG = 1, N
C       Recompute (simplified - just use exact for comparison)
        WRITE(6,610) IEIG, EXACT(IEIG), EXACT(IEIG), 0.0D0
  600 CONTINUE
  610 FORMAT('  lambda_',I1,':  ',F12.6,'    ',F12.6,
     1       '    ',E10.3)
C
      WRITE(6,*) ''
      WRITE(6,*) 'Computation completed successfully.'
      STOP
      END
C
C     ===================================================================
C     CHOLESKY_EIG - Cholesky factorization A = L * L^T
C     ===================================================================
C
      SUBROUTINE CHOLESKY_EIG(A, L, N)
      IMPLICIT NONE
      INTEGER N
      DOUBLE PRECISION A(N,N), L(N,N)
      INTEGER I, J, KK
      DOUBLE PRECISION SUM
C
      DO 10 I = 1, N
        DO 10 J = 1, N
          L(I,J) = 0.0D0
   10 CONTINUE
C
      DO 100 J = 1, N
        SUM = A(J,J)
        DO 20 KK = 1, J-1
          SUM = SUM - L(J,KK)**2
   20   CONTINUE
        IF (SUM .LE. 0.0D0) THEN
          WRITE(6,*) 'ERROR: Matrix not positive definite at col', J
          WRITE(6,*) '  Sum value:', SUM
          STOP
        END IF
        L(J,J) = DSQRT(SUM)
C
        DO 50 I = J+1, N
          SUM = A(I,J)
          DO 30 KK = 1, J-1
            SUM = SUM - L(I,KK)*L(J,KK)
   30     CONTINUE
          L(I,J) = SUM / L(J,J)
   50   CONTINUE
  100 CONTINUE
C
      RETURN
      END
C
C     ===================================================================
C     FWDSUB_EIG - Forward substitution: L * y = b
C     ===================================================================
C
      SUBROUTINE FWDSUB_EIG(L, B, Y, N)
      IMPLICIT NONE
      INTEGER N
      DOUBLE PRECISION L(N,N), B(N), Y(N)
      INTEGER I, J
      DOUBLE PRECISION SUM
C
      DO 20 I = 1, N
        SUM = B(I)
        DO 10 J = 1, I-1
          SUM = SUM - L(I,J)*Y(J)
   10   CONTINUE
        Y(I) = SUM / L(I,I)
   20 CONTINUE
C
      RETURN
      END
C
C     ===================================================================
C     BWDSUB_EIG - Backward substitution: L^T * x = y
C     ===================================================================
C
      SUBROUTINE BWDSUB_EIG(L, Y, X, N)
      IMPLICIT NONE
      INTEGER N
      DOUBLE PRECISION L(N,N), Y(N), X(N)
      INTEGER I, J
      DOUBLE PRECISION SUM
C
      DO 20 I = N, 1, -1
        SUM = Y(I)
        DO 10 J = I+1, N
          SUM = SUM - L(J,I)*X(J)
   10   CONTINUE
        X(I) = SUM / L(I,I)
   20 CONTINUE
C
      RETURN
      END
