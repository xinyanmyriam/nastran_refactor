C     ================================================================
C     TEST DRIVER B10: QRITER - QR Iteration for Tridiagonal
C                      Eigenvalue Problem
C     ================================================================
C
C     Test case: 5x5 tridiagonal matrix from spring-mass system
C
C         [2  -1   0   0   0]
C     T = [-1   2  -1   0   0]
C         [0  -1   2  -1   0]
C         [0   0  -1   2  -1]
C         [0   0   0  -1   2]
C
C     Analytical eigenvalues: lambda_k = 2 - 2*cos(k*pi/(N+1))
C     For N=5: [0.26795, 1.00000, 2.00000, 3.00000, 3.73205]
C
C     ================================================================
      PROGRAM TEST_QRITER
      IMPLICIT NONE
C
C     --- COMMON blocks matching QRITER expectations ---
C
C     /GIVN/ layout: IDUM0(100), N, LFREQ, IDUM3, IDUM4, HFREQ, LAMA,
C                    NV, NE, IDUM9, NFOUND, IDUM11, IDUM12, IDUM13,
C                    NEVER, MAX
C     Total words: 100 + 16 = 116
C
      INTEGER          IDUM0(100), NN, IDUM3, IDUM4, NLAMA, IDUM8
      INTEGER          NNV, NNE, IDUM9, NNFOUND
      INTEGER          IDUM11, IDUM12, IDUM13, NNEVER, NMAX
      REAL             RLFREQ, RHFREQ
      COMMON /GIVN  /  IDUM0, NN, RLFREQ, IDUM3, IDUM4, RHFREQ,
     1                 NLAMA, NNV, NNE, IDUM9, NNFOUND,
     2                 IDUM11, IDUM12, IDUM13, NNEVER, NMAX
C
      INTEGER          SYSBUF, NOUT
      COMMON /SYSTEM/  SYSBUF, NOUT
C
      INTEGER          IOPTN
      COMMON /REIGKR/  IOPTN
C
      DOUBLE PRECISION DLMDAS
      COMMON /MGIVXX/  DLMDAS
C
      CHARACTER        UFM*23, UWM*25, UIM*29
      COMMON /XMSSG /  UFM, UWM, UIM
C
C     --- Local variables ---
C
      INTEGER          NDIM
      PARAMETER        (NDIM = 5)
C
      DOUBLE PRECISION VAL(NDIM), O(NDIM)
      INTEGER          LOC(NDIM), QR
      DOUBLE PRECISION EXACT(NDIM), ERR, MAXERR, PI
      INTEGER          I, NPASS, NFAIL
      DOUBLE PRECISION TOL
      PARAMETER        (TOL = 1.0D-6)
C
C     ================================================================
C     SETUP
C     ================================================================
C
      WRITE(*,'(A)') '=============================================='
      WRITE(*,'(A)') ' B10 TEST DRIVER: QRITER'
      WRITE(*,'(A)') ' QR Iteration for Tridiagonal Eigenvalues'
      WRITE(*,'(A)') ' (from NASTRAN-95 MIS module)'
      WRITE(*,'(A)') '=============================================='
      WRITE(*,*)
C
C     Initialize COMMON blocks
C
      DO I = 1, 100
        IDUM0(I) = 0
      END DO
      NN      = NDIM
      RLFREQ  = 0.0
      RHFREQ  = 99.0
      NLAMA   = 0
      NNV     = NDIM
      NNE     = NDIM
      NNFOUND = 0
      NNEVER  = 0
      NMAX    = 500 * NDIM
      IDUM3   = 0
      IDUM4   = 0
      IDUM9   = 0
      IDUM11  = 0
      IDUM12  = 0
      IDUM13  = 0
C
      SYSBUF  = 100
      NOUT    = 6
      IOPTN   = 0
      DLMDAS  = 0.0D0
      UFM     = '*** USER FATAL MESSAGE'
      UWM     = '*** USER WARNING MESSAGE'
      UIM     = '*** USER INFORMATION MESSAGE'
C
C     ================================================================
C     TEST CASE: 5x5 Spring-Mass Tridiagonal Matrix
C     ================================================================
C
      WRITE(*,'(A)') ' Test: 5x5 spring-mass tridiagonal matrix'
      WRITE(*,'(A)') ' Diagonal = [2, 2, 2, 2, 2]'
      WRITE(*,'(A)') ' Off-diagonal = [-1, -1, -1, -1]'
      WRITE(*,*)
C
C     Fill diagonal (VAL) and squared off-diagonal (O)
C
      DO I = 1, NDIM
        VAL(I) = 2.0D0
      END DO
C
C     O(i) = square of off-diagonal element e(i)
C     e = [-1, -1, -1, -1], so O = [1, 1, 1, 1, 0]
C     Note: O(N) is not used as off-diagonal (only N-1 off-diags)
C
      DO I = 1, NDIM-1
        O(I) = 1.0D0
      END DO
      O(NDIM) = 0.0D0
C
C     QR = 0 means perform QR iteration (not just reorder)
C
      QR = 0
C
C     ================================================================
C     CALL QRITER
C     ================================================================
C
      WRITE(*,'(A)') ' Calling QRITER...'
      CALL QRITER(VAL, O, LOC, QR)
      WRITE(*,'(A)') ' QRITER returned successfully.'
      WRITE(*,*)
C
C     ================================================================
C     COMPUTE ANALYTICAL EIGENVALUES
C     ================================================================
C
      PI = 4.0D0 * DATAN(1.0D0)
      DO I = 1, NDIM
        EXACT(I) = 2.0D0 - 2.0D0 * DCOS(DBLE(I) * PI / DBLE(NDIM+1))
      END DO
C
C     ================================================================
C     COMPARE RESULTS
C     ================================================================
C
      WRITE(*,'(A)') ' Results:'
      WRITE(*,'(A)') ' -----------------------------------------------'
      WRITE(*,'(A5, A18, A18, A14, A6)') 
     1  ' k', 'Computed', 'Analytical', 'Error', 'Pass?'
      WRITE(*,'(A)') ' -----------------------------------------------'
C
      NPASS  = 0
      NFAIL  = 0
      MAXERR = 0.0D0
C
      DO I = 1, NDIM
        ERR = DABS(VAL(I) - EXACT(I))
        IF (ERR .GT. MAXERR) MAXERR = ERR
        IF (ERR .LT. TOL) THEN
          NPASS = NPASS + 1
          WRITE(*,'(I4, F18.10, F18.10, ES14.3, A6)') 
     1      I, VAL(I), EXACT(I), ERR, '  PASS'
        ELSE
          NFAIL = NFAIL + 1
          WRITE(*,'(I4, F18.10, F18.10, ES14.3, A6)') 
     1      I, VAL(I), EXACT(I), ERR, '  FAIL'
        END IF
      END DO
C
      WRITE(*,'(A)') ' -----------------------------------------------'
      WRITE(*,*)
      WRITE(*,'(A,I2,A,I2)') ' Eigenvalues passed: ', NPASS,
     1                        ' / ', NDIM
      WRITE(*,'(A,ES12.3)') ' Maximum absolute error: ', MAXERR
      WRITE(*,'(A,ES12.3)') ' Tolerance: ', TOL
      WRITE(*,*)
C
C     Print LOC array (original ordering)
C
      WRITE(*,'(A)') ' LOC array (original positions):'
      WRITE(*,'(5I6)') (LOC(I), I=1,NDIM)
      WRITE(*,*)
C
C     Print convergence info
C
      IF (NNEVER .EQ. 0) THEN
        WRITE(*,'(A)') ' Convergence: All eigenvalues converged.'
      ELSE
        WRITE(*,'(A,I4,A)') ' Convergence: ', NNEVER, 
     1    ' eigenvalue(s) had convergence issues.'
      END IF
      WRITE(*,*)
C
C     ================================================================
C     FINAL VERDICT
C     ================================================================
C
      IF (NFAIL .EQ. 0) THEN
        WRITE(*,'(A)') ' *** B10 QRITER TEST: ALL PASSED ***'
      ELSE
        WRITE(*,'(A,I2,A)') ' *** B10 QRITER TEST: ', NFAIL,
     1    ' EIGENVALUE(S) FAILED ***'
      END IF
      WRITE(*,'(A)') '=============================================='
C
      STOP
      END
