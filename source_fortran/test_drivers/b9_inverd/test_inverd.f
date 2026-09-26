      PROGRAM TEST_INVERD
C
C     Test driver for INVERD (NASTRAN-95 B9 module)
C     Tests matrix inversion and linear equation solving
C
C     Author: Test driver for verification
C     Date: 2025
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(4,4), AORIG(4,4), B(4,3), BCHECK(4,3)
      DIMENSION XKNOWN(4,3), INDEX(4,3)
      DOUBLE PRECISION DETERM
      INTEGER ISING, NDIM, N, M
      DOUBLE PRECISION ERR, MAXERR, TOL
      DOUBLE PRECISION PROD
C
      COMMON /MACHIN/ MACH
C
      PARAMETER (TOL = 1.0D-10)
      PARAMETER (NDIM = 4, N = 4)
C
      INTEGER NPASS, NFAIL
      NPASS = 0
      NFAIL = 0
C
      WRITE(*,*) '=============================================='
      WRITE(*,*) ' INVERD Test Driver (NASTRAN-95 B9 Module)'
      WRITE(*,*) '=============================================='
      WRITE(*,*)
C
C     Set MACH = 1 (required by INVERD COMMON /MACHIN/)
      MACH = 1
C
C     ===================================================================
C     TEST 1: Matrix Inversion Only (M=0)
C     ===================================================================
C
      WRITE(*,*) '----------------------------------------------'
      WRITE(*,*) ' TEST 1: Matrix Inversion (M=0)'
      WRITE(*,*) '----------------------------------------------'
      WRITE(*,*) ' 4x4 stiffness matrix K (4-spring system):'
      WRITE(*,*) ' K = [200 -100 0 0; -100 200 -100 0;'
      WRITE(*,*) '      0 -100 200 -100; 0 0 -100 100]'
      WRITE(*,*)
C
C     Define 4x4 SPD stiffness matrix (4-spring system)
      CALL SETK(A, NDIM)
C
C     Save original for verification
      DO 10 I = 1,N
        DO 10 J = 1,N
          AORIG(I,J) = A(I,J)
   10 CONTINUE
C
C     Call INVERD with M=0 (inversion only)
      M = 0
      ISING = 0
      CALL INVERD(NDIM, A, N, B, M, DETERM, ISING, INDEX)
C
      WRITE(*,*) ' ISING  = ', ISING
      WRITE(*,*) ' DETERM = ', DETERM
      WRITE(*,*)
C
      IF (ISING .NE. 1) THEN
        WRITE(*,*) ' ERROR: INVERD reported singular matrix!'
        NFAIL = NFAIL + 1
        GO TO 100
      END IF
C
C     Print K_inv
      WRITE(*,*) ' K_inverse:'
      DO 20 I = 1,N
        WRITE(*,30) (A(I,J), J=1,N)
   20 CONTINUE
   30 FORMAT(4X, 4(F12.6))
      WRITE(*,*)
C
C     Verify: K_inv * K_original = I
      WRITE(*,*) ' Verification: K_inv * K_original = I'
      MAXERR = 0.0D0
      DO 50 I = 1,N
        DO 50 J = 1,N
          PROD = 0.0D0
          DO 40 K = 1,N
            PROD = PROD + A(I,K) * AORIG(K,J)
   40     CONTINUE
          IF (I .EQ. J) THEN
            ERR = DABS(PROD - 1.0D0)
          ELSE
            ERR = DABS(PROD)
          END IF
          IF (ERR .GT. MAXERR) MAXERR = ERR
   50 CONTINUE
C
      WRITE(*,*) ' Max error in K_inv*K - I = ', MAXERR
      IF (MAXERR .LT. TOL) THEN
        WRITE(*,*) ' TEST 1: PASSED'
        NPASS = NPASS + 1
      ELSE
        WRITE(*,*) ' TEST 1: FAILED'
        NFAIL = NFAIL + 1
      END IF
      WRITE(*,*)
C
C     ===================================================================
C     TEST 2: Linear Solve with Single RHS (M=1)
C     ===================================================================
C
  100 CONTINUE
      WRITE(*,*) '----------------------------------------------'
      WRITE(*,*) ' TEST 2: Linear Solve (M=1)'
      WRITE(*,*) '----------------------------------------------'
      WRITE(*,*) ' Solve K*x = b where x_known = [1, 2, 3, 4]^T'
      WRITE(*,*)
C
C     Reset K matrix
      CALL SETK(A, NDIM)
C
C     Known solution x = [1, 2, 3, 4]
      XKNOWN(1,1) = 1.0D0
      XKNOWN(2,1) = 2.0D0
      XKNOWN(3,1) = 3.0D0
      XKNOWN(4,1) = 4.0D0
C
C     Compute RHS: b = K * x_known
      DO 110 I = 1,N
        B(I,1) = 0.0D0
        DO 110 J = 1,N
          B(I,1) = B(I,1) + A(I,J) * XKNOWN(J,1)
  110 CONTINUE
C
      WRITE(*,*) ' RHS b = K * [1,2,3,4]^T ='
      WRITE(*,120) (B(I,1), I=1,N)
  120 FORMAT(4X, 4(F10.2))
      WRITE(*,*)
C
C     Call INVERD with M=1
      M = 1
      ISING = 0
      CALL INVERD(NDIM, A, N, B, M, DETERM, ISING, INDEX)
C
      IF (ISING .NE. 1) THEN
        WRITE(*,*) ' ERROR: INVERD reported singular matrix!'
        NFAIL = NFAIL + 1
        GO TO 200
      END IF
C
C     Print solution
      WRITE(*,*) ' Solution x (should be [1, 2, 3, 4]):'
      WRITE(*,130) (B(I,1), I=1,N)
  130 FORMAT(4X, 4(F12.6))
      WRITE(*,*)
C
C     Verify solution
      MAXERR = 0.0D0
      DO 140 I = 1,N
        ERR = DABS(B(I,1) - XKNOWN(I,1))
        IF (ERR .GT. MAXERR) MAXERR = ERR
  140 CONTINUE
C
      WRITE(*,*) ' Max error = ', MAXERR
      IF (MAXERR .LT. TOL) THEN
        WRITE(*,*) ' TEST 2: PASSED'
        NPASS = NPASS + 1
      ELSE
        WRITE(*,*) ' TEST 2: FAILED'
        NFAIL = NFAIL + 1
      END IF
      WRITE(*,*)
C
C     ===================================================================
C     TEST 3: Multiple RHS (M=3)
C     ===================================================================
C
  200 CONTINUE
      WRITE(*,*) '----------------------------------------------'
      WRITE(*,*) ' TEST 3: Multiple RHS (M=3)'
      WRITE(*,*) '----------------------------------------------'
      WRITE(*,*) ' Three simultaneous RHS columns'
      WRITE(*,*)
C
C     Reset K matrix
      CALL SETK(A, NDIM)
C
C     Known solutions:
C     x1 = [1, 1, 1, 1]
C     x2 = [0, 1, 0, 1]
C     x3 = [1, 2, 3, 4]
      XKNOWN(1,1) = 1.0D0
      XKNOWN(2,1) = 1.0D0
      XKNOWN(3,1) = 1.0D0
      XKNOWN(4,1) = 1.0D0
C
      XKNOWN(1,2) = 0.0D0
      XKNOWN(2,2) = 1.0D0
      XKNOWN(3,2) = 0.0D0
      XKNOWN(4,2) = 1.0D0
C
      XKNOWN(1,3) = 1.0D0
      XKNOWN(2,3) = 2.0D0
      XKNOWN(3,3) = 3.0D0
      XKNOWN(4,3) = 4.0D0
C
C     Compute RHS columns: B(:,k) = K * x_known(:,k)
      DO 220 K = 1,3
        DO 210 I = 1,N
          B(I,K) = 0.0D0
          DO 210 J = 1,N
            B(I,K) = B(I,K) + A(I,J) * XKNOWN(J,K)
  210   CONTINUE
  220 CONTINUE
C
      WRITE(*,*) ' RHS columns b1, b2, b3:'
      DO 230 I = 1,N
        WRITE(*,240) B(I,1), B(I,2), B(I,3)
  230 CONTINUE
  240 FORMAT(4X, 3(F10.2))
      WRITE(*,*)
C
C     Call INVERD with M=3
      M = 3
      ISING = 0
      CALL INVERD(NDIM, A, N, B, M, DETERM, ISING, INDEX)
C
      IF (ISING .NE. 1) THEN
        WRITE(*,*) ' ERROR: INVERD reported singular matrix!'
        NFAIL = NFAIL + 1
        GO TO 300
      END IF
C
C     Print solutions
      WRITE(*,*) ' Solutions (columns should be [1,1,1,1],'
      WRITE(*,*) '            [0,1,0,1], [1,2,3,4]):'
      DO 250 I = 1,N
        WRITE(*,260) B(I,1), B(I,2), B(I,3)
  250 CONTINUE
  260 FORMAT(4X, 3(F12.6))
      WRITE(*,*)
C
C     Verify all three solutions
      MAXERR = 0.0D0
      DO 280 K = 1,3
        DO 270 I = 1,N
          ERR = DABS(B(I,K) - XKNOWN(I,K))
          IF (ERR .GT. MAXERR) MAXERR = ERR
  270   CONTINUE
  280 CONTINUE
C
      WRITE(*,*) ' Max error (all 3 columns) = ', MAXERR
      IF (MAXERR .LT. TOL) THEN
        WRITE(*,*) ' TEST 3: PASSED'
        NPASS = NPASS + 1
      ELSE
        WRITE(*,*) ' TEST 3: FAILED'
        NFAIL = NFAIL + 1
      END IF
      WRITE(*,*)
C
C     ===================================================================
C     SUMMARY
C     ===================================================================
C
  300 CONTINUE
      WRITE(*,*) '=============================================='
      WRITE(*,*) ' SUMMARY'
      WRITE(*,*) '=============================================='
      WRITE(*,*) ' Tests passed: ', NPASS
      WRITE(*,*) ' Tests failed: ', NFAIL
      IF (NFAIL .EQ. 0) THEN
        WRITE(*,*) ' ALL TESTS PASSED'
      ELSE
        WRITE(*,*) ' SOME TESTS FAILED'
      END IF
      WRITE(*,*) '=============================================='
C
      STOP
      END
C
C     ===================================================================
C     Utility subroutine: Set 4x4 stiffness matrix
C     ===================================================================
C
      SUBROUTINE SETK(A, NDIM)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(NDIM,NDIM)
C
C     4-spring system stiffness matrix:
C     K = [200 -100 0 0; -100 200 -100 0; 0 -100 200 -100; 0 0 -100 100]
C
      DO 10 I = 1,4
        DO 10 J = 1,4
          A(I,J) = 0.0D0
   10 CONTINUE
C
      A(1,1) = 200.0D0
      A(1,2) = -100.0D0
      A(2,1) = -100.0D0
      A(2,2) = 200.0D0
      A(2,3) = -100.0D0
      A(3,2) = -100.0D0
      A(3,3) = 200.0D0
      A(3,4) = -100.0D0
      A(4,3) = -100.0D0
      A(4,4) = 100.0D0
C
      RETURN
      END
C
C     ===================================================================
C     BLOCK DATA for COMMON /MACHIN/ initialization
C     ===================================================================
C
      BLOCK DATA MACHIN_INIT
      COMMON /MACHIN/ MACH
      DATA MACH /1/
      END
