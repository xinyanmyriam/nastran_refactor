C     TEST DRIVER FOR SOLVE - 6x6 Linear Equation Solver (B9)
C     ===========================================================
C     Standalone reference implementation of symmetric LU (Cholesky)
C     decomposition + forward/backward substitution.
C
C     This is the full 6x6 version matching the user specification:
C       - 6x6 SPD matrix from a spring system
C       - Known solution x = [1, 2, 3, 4, 5, 6]
C       - Computes b = A*x, then solves for x
C       - Prints L/U factors for reference
C
      PROGRAM SOLVE6X6
      IMPLICIT NONE
C
      INTEGER N
      PARAMETER (N = 6)
C
      DOUBLE PRECISION K(N,N), L(N,N), B(N), X(N), Y(N)
      DOUBLE PRECISION XKNOWN(N), ERR, MAXERR, BCHECK(N)
      INTEGER I, J, IDUM
C
      WRITE(6,*) 'SOLVE Test Driver - 6x6 Cholesky Decomposition'
      WRITE(6,*) '================================================'
      WRITE(6,*) ''
C
C     =========================================================
C     Build 6x6 stiffness matrix from a 5-spring chain system
C     with grounding springs to ensure SPD.
C
C     System: 6 DOFs, springs connecting adjacent DOFs:
C       k1=100 (DOF 1-2), k2=200 (DOF 2-3), k3=150 (DOF 3-4)
C       k4=175 (DOF 4-5), k5=125 (DOF 5-6)
C     Grounding springs: kg=80 at each DOF (anchored to ground)
C     =========================================================
C
      DO 10 I = 1, N
        DO 10 J = 1, N
          K(I,J) = 0.0D0
   10 CONTINUE
C
C     Grounding springs (80 at each DOF)
      DO 12 I = 1, N
        K(I,I) = 80.0D0
   12 CONTINUE
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
C     Spring 4: k=175, DOFs 4-5
      K(4,4) = K(4,4) + 175.0D0
      K(5,5) = K(5,5) + 175.0D0
      K(4,5) = K(4,5) - 175.0D0
      K(5,4) = K(5,4) - 175.0D0
C
C     Spring 5: k=125, DOFs 5-6
      K(5,5) = K(5,5) + 125.0D0
      K(6,6) = K(6,6) + 125.0D0
      K(5,6) = K(5,6) - 125.0D0
      K(6,5) = K(6,5) - 125.0D0
C
C     Print stiffness matrix
C
      WRITE(6,*) 'Stiffness Matrix K (6x6):'
      DO 20 I = 1, N
        WRITE(6,100) (K(I,J), J=1,N)
   20 CONTINUE
  100 FORMAT(6(2X,F8.1))
C
C     Known solution: x = [1, 2, 3, 4, 5, 6]
      XKNOWN(1) = 1.0D0
      XKNOWN(2) = 2.0D0
      XKNOWN(3) = 3.0D0
      XKNOWN(4) = 4.0D0
      XKNOWN(5) = 5.0D0
      XKNOWN(6) = 6.0D0
C
C     Compute RHS: b = K * x_known
      DO 30 I = 1, N
        B(I) = 0.0D0
        DO 30 J = 1, N
          B(I) = B(I) + K(I,J) * XKNOWN(J)
   30 CONTINUE
C
      WRITE(6,*) ''
      WRITE(6,*) 'Known solution x = [1, 2, 3, 4, 5, 6]'
      WRITE(6,*) 'RHS vector b = K*x:'
      WRITE(6,110) (B(I), I=1,N)
  110 FORMAT(6(2X,F10.2))
C
C     =========================================================
C     CHOLESKY DECOMPOSITION: K = L * L^T
C     =========================================================
C
      CALL CHOLESKY6(K, L, N)
C
      WRITE(6,*) ''
      WRITE(6,*) 'Lower triangular factor L (K = L * L^T):'
      DO 40 I = 1, N
        WRITE(6,120) (L(I,J), J=1,N)
   40 CONTINUE
  120 FORMAT(6(2X,F11.6))
C
C     Verify L*L^T = K
      WRITE(6,*) ''
      WRITE(6,*) 'Verification: L*L^T (should equal K):'
      DO 45 I = 1, N
        DO 44 J = 1, N
          BCHECK(J) = 0.0D0
          DO 43 IDUM = 1, N
            BCHECK(J) = BCHECK(J) + L(I,IDUM)*L(J,IDUM)
   43     CONTINUE
   44   CONTINUE
        WRITE(6,100) (BCHECK(J), J=1,N)
   45 CONTINUE
C
C     =========================================================
C     FORWARD SUBSTITUTION: L * y = b
C     =========================================================
C
      CALL FWDSUB6(L, B, Y, N)
C
      WRITE(6,*) ''
      WRITE(6,*) 'Intermediate vector y (L*y = b):'
      WRITE(6,130) (Y(I), I=1,N)
  130 FORMAT(6(2X,F12.6))
C
C     =========================================================
C     BACKWARD SUBSTITUTION: L^T * x = y
C     =========================================================
C
      CALL BWDSUB6(L, Y, X, N)
C
      WRITE(6,*) ''
      WRITE(6,*) 'Computed solution x:'
      WRITE(6,140) (X(I), I=1,N)
  140 FORMAT(6(2X,F12.8))
C
C     =========================================================
C     VERIFICATION
C     =========================================================
C
      WRITE(6,*) ''
      WRITE(6,*) 'Verification (Known vs Computed):'
      WRITE(6,*) '  DOF   Known x   Computed x      Error'
      WRITE(6,*) '  ---   -------   ----------      -----'
      MAXERR = 0.0D0
      DO 50 I = 1, N
        ERR = DABS(X(I) - XKNOWN(I))
        IF (ERR .GT. MAXERR) MAXERR = ERR
        WRITE(6,150) I, XKNOWN(I), X(I), ERR
   50 CONTINUE
  150 FORMAT(2X,I3,3X,F8.4,5X,F12.8,4X,E12.5)
C
      WRITE(6,*) ''
      WRITE(6,160) MAXERR
  160 FORMAT('  Maximum absolute error: ',E12.5)
C
C     Also check residual: r = b - K*x_computed
      WRITE(6,*) ''
      WRITE(6,*) 'Residual check (r = b - K*x_computed):'
      MAXERR = 0.0D0
      DO 60 I = 1, N
        ERR = B(I)
        DO 55 J = 1, N
          ERR = ERR - K(I,J) * X(J)
   55   CONTINUE
        IF (DABS(ERR) .GT. MAXERR) MAXERR = DABS(ERR)
        WRITE(6,170) I, ERR
   60 CONTINUE
  170 FORMAT('  r(',I1,') = ',E12.5)
C
      WRITE(6,*) ''
      WRITE(6,180) MAXERR
  180 FORMAT('  Maximum residual: ',E12.5)
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
C     CHOLESKY6 - Cholesky factorization A = L * L^T
C     ===================================================================
C
      SUBROUTINE CHOLESKY6(A, L, N)
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
C     FWDSUB6 - Forward substitution: L * y = b
C     ===================================================================
C
      SUBROUTINE FWDSUB6(L, B, Y, N)
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
C     BWDSUB6 - Backward substitution: L^T * x = y
C     ===================================================================
C
      SUBROUTINE BWDSUB6(L, Y, X, N)
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
