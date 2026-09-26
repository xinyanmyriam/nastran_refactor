C ============================================================
C     TEST_TRD1C_NATIVE.F - Main Driver Program
C
C     Runs NASTRAN's actual TRD1C subroutine with a mini-GINO
C     runtime to perform transient analysis on a 3-DOF
C     spring-mass system using central difference method.
C
C     Test Problem:
C       3 masses (m=1.0 each) connected by springs (k=1.0)
C       Fixed-fixed boundary: wall--k--m1--k--m2--k--m3--k--wall
C       K = [2 -1  0; -1  2 -1;  0 -1  2]
C       M = I (identity, lumped mass)
C       C = 0 (no damping)
C       IC: u0 = [1, 0, 0], udot0 = [0, 0, 0]
C       dt = 0.01 (well within stability limit)
C
C     Central Difference Scheme:
C       (M/dt^2) * u_{n+1} = P*_n
C       where P*_n = F_n + (2M/dt^2 - K)*u_n + (-M/dt^2)*u_{n-1}
C
C       Effective stiffness: K_eff = M/dt^2 = 10000*I
C       SCR1 matrix (for u_n): 2M/dt^2 - K
C       SCR4 matrix (for u_{n-1}): -M/dt^2
C
C     Expected: Free vibration with frequencies:
C       omega_1 = 2*sin(pi/8) = 0.7654
C       omega_2 = 2*sin(2*pi/8) = sqrt(2) = 1.4142
C       omega_3 = 2*sin(3*pi/8) = 1.8478
C
C     Author: Mini-NASTRAN Runtime Experiment
C     Date: 2025
C ============================================================
C
      PROGRAM TEST_TRD1C
      IMPLICIT NONE
C
C     --- NASTRAN COMMON BLOCKS (must match trd1c.f expectations) ---
C
C     /BLANK/ - first 5 words used, including NCOL at position 5
C               and ISTART at position 6
      REAL    BLANK_R(200)
      INTEGER BLANK_I(200)
      COMMON /BLANK/ BLANK_R
      EQUIVALENCE (BLANK_R, BLANK_I)
C
C     /SYSTEM/ - SYSBUF at position 1, NNOUT at 2, ISYSTM(79) at 3-81,
C                ICPFLG at 82
      INTEGER SYSTEM_I(100)
      COMMON /SYSTEM/ SYSTEM_I
C
C     /ZZZZZZ/ - Open core (big work array)
      REAL    Z(100000)
      INTEGER IZ(100000)
      COMMON /ZZZZZZ/ Z
      EQUIVALENCE (Z, IZ)
C
C     /PACKX/ - Packing control
      INTEGER PACKX_IT1, PACKX_IT2, PACKX_II, PACKX_JJ, PACKX_INCR
      COMMON /PACKX/ PACKX_IT1, PACKX_IT2, PACKX_II, PACKX_JJ,
     1               PACKX_INCR
C
C     /UNPAKX/ - Unpacking control
      INTEGER UNPAKX_IT3, UNPAKX_III, UNPAKX_JJJ, UNPAKX_INCR1
      COMMON /UNPAKX/ UNPAKX_IT3, UNPAKX_III, UNPAKX_JJJ,
     1                UNPAKX_INCR1
C
C     /TRDXX/ - TRD control: IK(7), IDUM(14), ISCR1-6, IOPEN, ISYM, 
C               TO, NOPD, ISPNL
C     Layout: IK(7)=words 1-7, IM(7)=words 8-14 (IDUM uses 8-21),
C             IB(7)=words 15-21, ISCR1=22, ISCR2=23, ISCR3=24,
C             ISCR4=25, ISCR5=26, ISCR6=27, IOPEN=28, ISYM=29,
C             TO=30, NOPD=31, ISPNL=32
      INTEGER TRDXX(32)
      REAL    TRDXX_R(32)
      COMMON /TRDXX/ TRDXX
      EQUIVALENCE (TRDXX, TRDXX_R)
C
C     /TRDD1/ - TRD1C internal state  
      INTEGER TRDD1(50)
      REAL    TRDD1_R(50)
      COMMON /TRDD1/ TRDD1
      EQUIVALENCE (TRDD1, TRDD1_R)
C
C     /NAMES/ - NASTRAN constant names
C     RD, RDREW, WRT, WRTREW, REW, NOREW, EOFNRW, RSP, RDP,
C     CSP, CDP, SQR, RECT, DIAG, LOWTRI, UPRTRI, SYM, ROW, IDENTY
      INTEGER NAMES(19)
      COMMON /NAMES/ NAMES
C
C     /INFBSX/ - FBS file control
      INTEGER INFBSX(14)
      COMMON /INFBSX/ INFBSX
C
C     /MACHIN/ - Machine type
      INTEGER MACHIN_MACH
      COMMON /MACHIN/ MACHIN_MACH
C
C     --- MINI-GINO COMMON ---
      INTEGER NROW_G, FILE_POS, FILE_NCOL, MCB_STORE
      REAL    MAT_STORE, VEC_STORE
      COMMON /MINIGINO/
     1  MAT_STORE(3,3,10),
     2  VEC_STORE(3,30),
     3  MCB_STORE(7,30),
     4  FILE_POS(30),
     5  FILE_NCOL(30),
     6  NROW_G
C
C     --- LOCAL VARIABLES ---
      INTEGER I, J, N, NSTEPS, NOUT, NGROUP
      REAL    DT, KMAT(3,3), MMAT(3,3), KEFF(3,3), KEFF_INV(3,3)
      REAL    A1(3,3), A2(3,3)
      REAL    U0(3), UDOT0(3)
      REAL    PI, OMEGA(3), PHI(3,3), T, UEXACT
      REAL    DETERM
      DOUBLE PRECISION DKEFF(3,3), DB(3,1), DDETERM
      INTEGER ISING, INDEX(3,3)
      INTEGER IC_FILE, PD_FILE, UDV_FILE
      INTEGER SCR1_FILE, SCR2_FILE, SCR3_FILE, SCR4_FILE, SCR5_FILE
C
C     ============================================================
C     PROBLEM PARAMETERS
C     ============================================================
      N = 3
      DT = 0.01
      NSTEPS = 100
      NOUT = 10
      NGROUP = 1
C
C     File numbers
      IC_FILE  = 101
      PD_FILE  = 102
      UDV_FILE = 201
      SCR1_FILE = 301
      SCR2_FILE = 302
      SCR3_FILE = 303
      SCR4_FILE = 304
      SCR5_FILE = 305
C
      WRITE(*,*) '=============================================='
      WRITE(*,*) ' B12: TRD1C Native - Mini-NASTRAN Runtime'
      WRITE(*,*) '=============================================='
      WRITE(*,*) ''
      WRITE(*,*) 'Test: 3-DOF spring-mass, central difference'
      WRITE(*,*) 'N =', N, '  dt =', DT, '  steps =', NSTEPS
      WRITE(*,*) ''
C
C     ============================================================
C     SETUP SYSTEM MATRICES
C     ============================================================
C
C     Stiffness: K = [2 -1 0; -1 2 -1; 0 -1 2]
      DO I = 1, N
        DO J = 1, N
          KMAT(I,J) = 0.0
          MMAT(I,J) = 0.0
        ENDDO
        KMAT(I,I) = 2.0
        MMAT(I,I) = 1.0
      ENDDO
      KMAT(1,2) = -1.0
      KMAT(2,1) = -1.0
      KMAT(2,3) = -1.0
      KMAT(3,2) = -1.0
C
C     Initial conditions
      U0(1) = 1.0
      U0(2) = 0.0
      U0(3) = 0.0
      UDOT0(1) = 0.0
      UDOT0(2) = 0.0
      UDOT0(3) = 0.0
C
C     ============================================================
C     COMPUTE EFFECTIVE MATRICES FOR CENTRAL DIFFERENCE
C     ============================================================
C
C     K_eff = M/dt^2 (effective stiffness to decompose/invert)
C     SCR1 = 2*M/dt^2 - K (multiply with u_n, result ADDED to P)
C     SCR4 = -M/dt^2 (multiply with u_{n-1}, result ADDED to P)
C
C     Note: In TRD1C, STEP does:
C       P += SCR1 * u_n   (this adds the 2M/dt^2 - K contribution)
C       P += SCR4 * u_{n-1} (this adds the -M/dt^2 contribution)
C       Then solves K_eff * u_{n+1} = P
C
      DO I = 1, N
        DO J = 1, N
          KEFF(I,J) = MMAT(I,J) / (DT*DT)
          A1(I,J) = 2.0*MMAT(I,J)/(DT*DT) - KMAT(I,J)
          A2(I,J) = -MMAT(I,J) / (DT*DT)
        ENDDO
      ENDDO
C
C     Compute K_eff inverse using INVERD
      DO I = 1, N
        DO J = 1, N
          DKEFF(I,J) = DBLE(KEFF(I,J))
        ENDDO
      ENDDO
      DB(1,1) = 0.0D0
      ISING = -1
      CALL INVERD(3, DKEFF, 3, DB, 0, DDETERM, ISING, INDEX)
      IF (ISING .EQ. 2) THEN
        WRITE(*,*) 'ERROR: K_eff is singular!'
        STOP
      ENDIF
C
C     Store K_eff_inv as single precision in matrix slot 9
      DO I = 1, N
        DO J = 1, N
          KEFF_INV(I,J) = REAL(DKEFF(I,J))
        ENDDO
      ENDDO
C
      WRITE(*,*) 'K_eff (M/dt^2):'
      DO I = 1, N
        WRITE(*,'(A,3F12.2)') '  ', (KEFF(I,J), J=1,N)
      ENDDO
      WRITE(*,*) ''
      WRITE(*,*) 'SCR1 (2M/dt^2 - K):'
      DO I = 1, N
        WRITE(*,'(A,3F12.2)') '  ', (A1(I,J), J=1,N)
      ENDDO
      WRITE(*,*) ''
      WRITE(*,*) 'K_eff_inv (diagonal 1/10000):'
      DO I = 1, N
        WRITE(*,'(A,3E14.6)') '  ', (KEFF_INV(I,J), J=1,N)
      ENDDO
      WRITE(*,*) ''
C
C     ============================================================
C     STORE MATRICES IN MINI-GINO
C     ============================================================
C
C     Slot 1 (SCR1=301): A1 = 2M/dt^2 - K
      DO I = 1, N
        DO J = 1, N
          MAT_STORE(I,J,1) = A1(I,J)
        ENDDO
      ENDDO
C
C     Slot 4 (SCR4=304): A2 = -M/dt^2
      DO I = 1, N
        DO J = 1, N
          MAT_STORE(I,J,4) = A2(I,J)
        ENDDO
      ENDDO
C
C     Slot 9: K_eff_inv (used by INTFBS/FBSINT)
      DO I = 1, N
        DO J = 1, N
          MAT_STORE(I,J,9) = KEFF_INV(I,J)
        ENDDO
      ENDDO
C
C     Slot 7 (IK=401): Stiffness K (used by FORM1/FORM2 via MATVEC)
      DO I = 1, N
        DO J = 1, N
          MAT_STORE(I,J,7) = KMAT(I,J)
        ENDDO
      ENDDO
C
C     Slot 8 (IM=402): Mass M (used by FORM2)
      DO I = 1, N
        DO J = 1, N
          MAT_STORE(I,J,8) = MMAT(I,J)
        ENDDO
      ENDDO
C
C     Slot 10 (IB=403): Damping B = 0 (used by FORM1/FORM2)
      DO I = 1, N
        DO J = 1, N
          MAT_STORE(I,J,10) = 0.0
        ENDDO
      ENDDO
C
C     ============================================================
C     STORE INITIAL CONDITIONS IN "IC FILE" (file 101)
C     ============================================================
C
C     IC file stores: record 1 = u0, record 2 = udot0
      DO I = 1, N
        VEC_STORE(I, 1) = U0(I)
        VEC_STORE(I, 2) = UDOT0(I)
      ENDDO
      FILE_NCOL(1) = 2
      FILE_POS(1)  = 1
C
C     ============================================================
C     STORE LOAD VECTORS IN "PD FILE" (file 102)
C     ============================================================
C
C     PD stores load vectors: first skip 1 record (header),
C     then P0, P1, P2, ...
C     For free vibration: all loads = 0
C     Need at least NSTEPS+2 zero vectors
C     Vector slot for PD starts at slot 3
      DO J = 1, 20
        DO I = 1, N
          VEC_STORE(I, 2+J) = 0.0
        ENDDO
      ENDDO
      FILE_NCOL(3) = 20
      FILE_POS(3) = 1
C
C     ============================================================
C     SETUP UDV OUTPUT FILE (file 201)
C     ============================================================
C     Starts empty, TRD1C will PACK results into it
      FILE_NCOL(5) = 0
      FILE_POS(5) = 1
C
C     ============================================================
C     SETUP MCBs (Matrix Control Blocks)
C     ============================================================
C
C     MCB for SCR2 (L factor) - needed by STEP/RDTRL
      MCB_STORE(1,16) = SCR2_FILE
      MCB_STORE(2,16) = N
      MCB_STORE(3,16) = N
      MCB_STORE(4,16) = 6
      MCB_STORE(5,16) = 1
      MCB_STORE(6,16) = 0
      MCB_STORE(7,16) = 0
C
C     MCB for SCR3 (U factor)
      MCB_STORE(1,17) = SCR3_FILE
      MCB_STORE(2,17) = N
      MCB_STORE(3,17) = N
      MCB_STORE(4,17) = 7
      MCB_STORE(5,17) = 1
      MCB_STORE(6,17) = 0
      MCB_STORE(7,17) = 0
C
C     ============================================================
C     SETUP NASTRAN COMMON BLOCKS
C     ============================================================
C
C     /BLANK/: position 4 = ?, position 5 = NCOL (initially 0)
C              position 6 = ISTART
      DO I = 1, 200
        BLANK_I(I) = 0
      ENDDO
      BLANK_I(5) = 0
      BLANK_I(6) = -1
C
C     /SYSTEM/: SYSBUF=100, NNOUT=6, ICPFLG=0
      DO I = 1, 100
        SYSTEM_I(I) = 0
      ENDDO
      SYSTEM_I(1) = 100
      SYSTEM_I(2) = 6
      SYSTEM_I(82) = 0
C
C     /NAMES/: Set standard NASTRAN constants
C     RD=0, RDREW=1, WRT=2, WRTREW=3, REW=1, NOREW=0, EOFNRW=7
C     RSP=1, RDP=2, CSP=3, CDP=4, SQR=1, RECT=2, DIAG=6
C     LOWTRI=5, UPRTRI=7, SYM=8, ROW=9, IDENTY=3
      NAMES(1)  = 0
      NAMES(2)  = 1
      NAMES(3)  = 2
      NAMES(4)  = 3
      NAMES(5)  = 1
      NAMES(6)  = 0
      NAMES(7)  = 7
      NAMES(8)  = 1
      NAMES(9)  = 2
      NAMES(10) = 3
      NAMES(11) = 4
      NAMES(12) = 1
      NAMES(13) = 2
      NAMES(14) = 6
      NAMES(15) = 5
      NAMES(16) = 7
      NAMES(17) = 8
      NAMES(18) = 9
      NAMES(19) = 3
C
C     /TRDXX/: Set up file IDs
C     IK(7) = MCB for effective stiffness K
C     Word layout: IK(1-7), then IFILM(1-7)=words 8-14,
C                  IFILB(1-7)=words 15-21
C     Then: ISCR1=22, ISCR2=23, ISCR3=24, ISCR4=25,
C           ISCR5=26, ISCR6=27, IOPEN=28, ISYM=29
C           TO=30(real), NOPD=31(logical), ISPNL=32
C
C     IK MCB (stiffness file 401)
      TRDXX(1) = 401
      TRDXX(2) = N
      TRDXX(3) = N
      TRDXX(4) = 1
      TRDXX(5) = 1
      TRDXX(6) = 0
      TRDXX(7) = 0
C     IM MCB (mass file 402) -- words 8-14
      TRDXX(8)  = 402
      TRDXX(9)  = N
      TRDXX(10) = N
      TRDXX(11) = 1
      TRDXX(12) = 1
      TRDXX(13) = 0
      TRDXX(14) = 0
C     IB MCB (damping file 403) -- words 15-21
      TRDXX(15) = 403
      TRDXX(16) = N
      TRDXX(17) = N
      TRDXX(18) = 1
      TRDXX(19) = 1
      TRDXX(20) = 0
      TRDXX(21) = 0
C     Scratch files
      TRDXX(22) = SCR1_FILE
      TRDXX(23) = SCR2_FILE
      TRDXX(24) = SCR3_FILE
      TRDXX(25) = SCR4_FILE
      TRDXX(26) = SCR5_FILE
      TRDXX(27) = 306
C     IOPEN = 0 (files not open initially)
      TRDXX(28) = 0
C     ISYM = 1 (use INTFBS, not FBSINT)
      TRDXX(29) = 1
C     TO = 0.0 (start time)
      TRDXX_R(30) = 0.0
C     NOPD = .FALSE. (we have loads, even if zero)
      TRDXX(31) = 0
C     ISPNL = 0
      TRDXX(32) = 0
C
C     /MACHIN/ - Machine type (1=generic)
      MACHIN_MACH = 1
C
C     ============================================================
C     SETUP OPEN CORE (Z array) FOR TRD1C
C     ============================================================
C
C     TRD1C uses Z as working storage. It accesses:
C     - Z(IGROUP) through Z(IGROUP+3*NGROUP-1) for time step groups
C     - Various buffer areas
C     - Vector storage areas allocated from bottom
C
C     The group data is stored at the END of available core.
C     NZ = KORSZ(Z) = 100000
C     IGROUP = NZ - 3*NGROUP + 1 = 99998
C     Group 1: NSTEP, DELTAT, NOUT
C
C     Store group parameters at the end of Z
C     Z(99998) = NSTEP (integer)
C     Z(99999) = DELTAT (real)
C     Z(100000)= NOUT (integer)
      IZ(99998) = NSTEPS
      Z(99999)  = DT
      IZ(100000)= NOUT
C
C     ============================================================
C     CALL TRD1C
C     ============================================================
C
      WRITE(*,*) 'Calling TRD1C...'
      WRITE(*,*) '  IC=', IC_FILE, ' PD=', PD_FILE
      WRITE(*,*) '  NGROUP=', NGROUP, ' NLFTP=0'
      WRITE(*,*) '  UDV=', UDV_FILE, ' ILOOP=1'
      WRITE(*,*) ''
C
      CALL TRD1C(IC_FILE, PD_FILE, NGROUP, 0, UDV_FILE, 1,
     1           SCR5_FILE, 0, 0, 0, 0, 0, 0)
C
      WRITE(*,*) ''
      WRITE(*,*) 'TRD1C returned successfully!'
      WRITE(*,*) ''
C
C     ============================================================
C     EXTRACT AND PRINT RESULTS FROM UDV FILE
C     ============================================================
C
C     UDV contains triplets: (u, udot, uddot) at each output time
C     Output every NOUT steps, so we have NSTEPS/NOUT output times
C     Plus the initial output
C
      WRITE(*,*) '=============================================='
      WRITE(*,*) ' DISPLACEMENT HISTORY (from UDV file)'
      WRITE(*,*) '=============================================='
      WRITE(*,*) ''
      WRITE(*,'(A)') '  Step    Time       u1          u2'//
     1               '          u3'
      WRITE(*,'(A)') '  ----  ------  ----------  ----------'//
     1               '  ----------'
C
C     Read back from UDV vector storage (slot 5)
C     Each output has 3 columns: u, udot, uddot
      J = 0
      DO I = 1, FILE_NCOL(5), 3
        T = J * NOUT * DT
        WRITE(*,'(I6,F8.4,3F12.6)') J*NOUT, T,
     1    VEC_STORE(1,4+I), VEC_STORE(2,4+I), VEC_STORE(3,4+I)
        J = J + 1
      ENDDO
C
C     ============================================================
C     ANALYTICAL SOLUTION COMPARISON
C     ============================================================
C
      WRITE(*,*) ''
      WRITE(*,*) '=============================================='
      WRITE(*,*) ' ANALYTICAL COMPARISON (modal superposition)'
      WRITE(*,*) '=============================================='
      WRITE(*,*) ''
C
C     Natural frequencies for fixed-fixed 3-DOF chain:
C     omega_k = 2*sin(k*pi/(2*(N+1))), k=1,2,3
      PI = 4.0*ATAN(1.0)
      DO I = 1, N
        OMEGA(I) = 2.0*SIN(REAL(I)*PI/REAL(2*(N+1)))
      ENDDO
      WRITE(*,*) 'Natural frequencies:'
      WRITE(*,'(A,F10.6)') '  omega_1 = ', OMEGA(1)
      WRITE(*,'(A,F10.6)') '  omega_2 = ', OMEGA(2)
      WRITE(*,'(A,F10.6)') '  omega_3 = ', OMEGA(3)
      WRITE(*,*) ''
C
C     Mode shapes: phi_k(i) = sin(i*k*pi/(N+1))
      DO I = 1, N
        DO J = 1, N
          PHI(I,J) = SIN(REAL(I)*REAL(J)*PI/REAL(N+1))
        ENDDO
      ENDDO
C
C     Analytical u1(t) at final output time
      T = NSTEPS * DT
      WRITE(*,'(A,F8.4,A)') ' At t = ', T, ':'
C
C     u(t) = sum_k [ (phi_k . u0)/(phi_k . phi_k) * phi_k * cos(omega_k*t) ]
C     For N=3, phi_k.phi_k = 2.0 (normalized), phi_k.u0 = phi_k(1)
      UEXACT = 0.0
      DO I = 1, N
        UEXACT = UEXACT + (PHI(1,I)/2.0) * PHI(1,I) * COS(OMEGA(I)*T)
      ENDDO
      WRITE(*,'(A,F12.6)') '  Analytical u1 = ', UEXACT
C
C     Compare with numerical (last displacement output)
      IF (FILE_NCOL(5) .GE. 3) THEN
        I = FILE_NCOL(5) - 2
        WRITE(*,'(A,F12.6)') '  Numerical  u1 = ', VEC_STORE(1,4+I)
        WRITE(*,'(A,E12.4)') '  Error         = ',
     1    ABS(VEC_STORE(1,4+I) - UEXACT)
      ENDIF
C
      WRITE(*,*) ''
      WRITE(*,*) '=============================================='
      WRITE(*,*) ' TEST COMPLETE'
      WRITE(*,*) '=============================================='
C
      STOP
      END
C
C ============================================================
C     KORSZ - Function version (TRD1C calls it as a function)
C     Returns available core size
C ============================================================
      INTEGER FUNCTION KORSZ (Z)
      REAL Z(1)
      KORSZ = 100000
      RETURN
      END
