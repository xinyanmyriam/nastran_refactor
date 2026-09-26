C     TEST DRIVER FOR SOLVE - Linear Equation Solver (B9)
C     =====================================================
C     Standalone reference implementation of symmetric LU (Cholesky)
C     decomposition + forward/backward substitution.
C
C     solve.f in NASTRAN is a DMAP-level driver using GINO file I/O,
C     so this standalone program implements the same core algorithm:
C       1. Cholesky factorization: A = L * L^T
C       2. Forward substitution: L * y = b
C       3. Backward substitution: L^T * x = y
C
C     Test system: 4x4 SPD matrix from 3-spring FEA model
C       Springs: k1=100, k2=200, k3=150 connecting 4 DOFs
C       Known solution: x = [1, 2, 3, 4]
C       RHS: b = K * x_known
C
      PROGRAM TESTSOLVE
      IMPLICIT NONE
C
      INTEGER N
      PARAMETER (N = 4)
C
      DOUBLE PRECISION K(N,N), L(N,N), B(N), X(N), Y(N)
      DOUBLE PRECISION XKNOWN(N), ERR, MAXERR
      INTEGER I, J
C
      WRITE(6,*) 'SOLVE Test Driver - Cholesky Decomposition'
      WRITE(6,*) '============================================'
      WRITE(6,*) ''
C
C     Build stiffness matrix from 3-spring system
C     Spring 1 (k=100): connects DOF 1-2
C     Spring 2 (k=200): connects DOF 2-3
C     Spring 3 (k=150): connects DOF 3-4
C     Plus grounding springs: k_g=50 at each DOF (to make SPD)
C
C     K(i,i) = sum of springs at DOF i + grounding
C     K(i,j) = -k_spring between i and j
C
      DO 10 I = 1, N
        DO 10 J = 1, N
          K(I,J) = 0.0D0
   10 CONTINUE
C
C     Grounding springs (50 at each DOF)
      K(1,1) = 50.0D0
      K(2,2) = 50.0D0
      K(3,3) = 50.0D0
      K(4,4) = 50.0D0
C
C     Spring 1: k=100, DOFs 1-2
      K(1,1) = K(1,1) + 100.0D0
      K(2,2) = K(2,2) + 100.0D0
      K(1,2) = K(1,2) - 100.0D0
      K(2,1) = K(2,1) - 100.0D0
C
C     Spring 2: k=200, DOFs 2-3
      K(2,2) = K(2,2) + 200.0D0
      K(3,3) = K(3,3) + 200.0D0
      K(2,3) = K(2,3) - 200.0D0
      K(3,2) = K(3,2) - 200.0D0
C
C     Spring 3: k=150, DOFs 3-4
      K(3,3) = K(3,3) + 150.0D0
      K(4,4) = K(4,4) + 150.0D0
      K(3,4) = K(3,4) - 150.0D0
      K(4,3) = K(4,3) - 150.0D0
C
C     Print stiffness matrix
C
      WRITE(6,*) 'Stiffness Matrix K:'
      DO 20 I = 1, N
        WRITE(6,100) (K(I,J), J=1,N)
   20 CONTINUE
  100 FORMAT(4(2X,F10.2))
C
C     Known solution
      XKNOWN(1) = 1.0D0
      XKNOWN(2) = 2.0D0
      XKNOWN(3) = 3.0D0
      XKNOWN(4) = 4.0D0
C
C     Compute RHS: b = K * x_known
      DO 30 I = 1, N
        B(I) = 0.0D0
        DO 30 J = 1, N
          B(I) = B(I) + K(I,J) * XKNOWN(J)
   30 CONTINUE
C
      WRITE(6,*) ''
      WRITE(6,*) 'Known solution x = [1, 2, 3, 4]'
      WRITE(6,*) 'RHS vector b = K*x:'
      WRITE(6,110) (B(I), I=1,N)
  110 FORMAT(4(2X,F10.2))
C
C     =========================================================
C     CHOLESKY DECOMPOSITION: K = L * L^T
C     =========================================================
C
      CALL CHOLESKY(K, L, N)
C
      WRITE(6,*) ''
      WRITE(6,*) 'Lower triangular factor L:'
      DO 40 I = 1, N
        WRITE(6,120) (L(I,J), J=1,N)
   40 CONTINUE
  120 FORMAT(4(2X,F12.6))
C
C     =========================================================
C     FORWARD SUBSTITUTION: L * y = b
C     =========================================================
C
      CALL FWDSUB(L, B, Y, N)
C
      WRITE(6,*) ''
      WRITE(6,*) 'Intermediate vector y (L*y=b):'
      WRITE(6,130) (Y(I), I=1,N)
  130 FORMAT(4(2X,F12.6))
C
C     =========================================================
C     BACKWARD SUBSTITUTION: L^T * x = y
C     =========================================================
C
      CALL BWDSUB(L, Y, X, N)
C
      WRITE(6,*) ''
      WRITE(6,*) 'Computed solution x:'
      WRITE(6,140) (X(I), I=1,N)
  140 FORMAT(4(2X,F12.8))
C
C     =========================================================
C     VERIFICATION
C     =========================================================
C
      WRITE(6,*) ''
      WRITE(6,*) 'Verification:'
      WRITE(6,*) '  Known x:    Computed x:   Error:'
      MAXERR = 0.0D0
      DO 50 I = 1, N
        ERR = ABS(X(I) - XKNOWN(I))
        IF (ERR .GT. MAXERR) MAXERR = ERR
        WRITE(6,150) XKNOWN(I), X(I), ERR
   50 CONTINUE
  150 FORMAT(2X,F10.6,4X,F10.6,4X,E12.5)
C
      WRITE(6,*) ''
      WRITE(6,160) MAXERR
  160 FORMAT('  Maximum error: ',E12.5)
C
      IF (MAXERR .LT. 1.0D-10) THEN
        WRITE(6,*) '  PASS: Solution matches to machine precision.'
      ELSE IF (MAXERR .LT. 1.0D-6) THEN
        WRITE(6,*) '  PASS: Solution matches within tolerance.'
      ELSE
        WRITE(6,*) '  FAIL: Solution does not match!'
      END IF
C
      WRITE(6,*) ''
      WRITE(6,*) 'Computation completed successfully.'
      STOP
      END
C
C     ===================================================================
C     CHOLESKY - Cholesky factorization A = L * L^T
C     ===================================================================
C
      SUBROUTINE CHOLESKY(A, L, N)
      IMPLICIT NONE
      INTEGER N
      DOUBLE PRECISION A(N,N), L(N,N)
      INTEGER I, J, K
      DOUBLE PRECISION SUM
C
C     Zero L
      DO 10 I = 1, N
        DO 10 J = 1, N
          L(I,J) = 0.0D0
   10 CONTINUE
C
C     Compute L column by column
      DO 100 J = 1, N
C       Diagonal element
        SUM = A(J,J)
        DO 20 K = 1, J-1
          SUM = SUM - L(J,K)**2
   20   CONTINUE
        IF (SUM .LE. 0.0D0) THEN
          WRITE(6,*) 'ERROR: Matrix not positive definite at col', J
          STOP
        END IF
        L(J,J) = DSQRT(SUM)
C
C       Off-diagonal elements
        DO 50 I = J+1, N
          SUM = A(I,J)
          DO 30 K = 1, J-1
            SUM = SUM - L(I,K)*L(J,K)
   30     CONTINUE
          L(I,J) = SUM / L(J,J)
   50   CONTINUE
  100 CONTINUE
C
      RETURN
      END
C
C     ===================================================================
C     FWDSUB - Forward substitution: L * y = b
C     ===================================================================
C
      SUBROUTINE FWDSUB(L, B, Y, N)
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
C     BWDSUB - Backward substitution: L^T * x = y
C     ===================================================================
C
      SUBROUTINE BWDSUB(L, Y, X, N)
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
