C ============================================================
C     GINO_SIM.F - Mini-GINO Simulation Layer
C
C     Provides in-memory replacements for NASTRAN's GINO I/O
C     subsystem, sufficient to run TRD1C for a small dense system.
C
C     Key idea: Instead of tape/disk files, store matrices and
C     vectors in COMMON /MINIGINO/ arrays. File operations become
C     array index manipulations.
C
C     Author: Mini-NASTRAN Runtime Experiment
C     Date: 2025
C ============================================================
C
C     ---- COMMON BLOCK FOR MINI-GINO STORAGE ----
C
      BLOCK DATA MINIGINO_INIT
      IMPLICIT NONE
      INTEGER NROW_G, FILE_POS, FILE_NCOL, MCB_STORE
      REAL    MAT_STORE, VEC_STORE
      COMMON /MINIGINO/
     1  MAT_STORE(3,3,10),
     2  VEC_STORE(3,30),
     3  MCB_STORE(7,30),
     4  FILE_POS(30),
     5  FILE_NCOL(30),
     6  NROW_G
      DATA NROW_G /3/
      DATA FILE_POS /30*1/
      DATA FILE_NCOL /30*0/
      END
C
C ============================================================
C     MATVEC - Matrix-Vector Multiply (REPLACES NASTRAN's MATVEC)
C
C     X = X + A*Y where A is identified by FILEA(1)
C
C     In mini-GINO, we look up which matrix slot corresponds
C     to the file number in FILEA(1), then do dense multiply.
C ============================================================
C
      SUBROUTINE MATVEC (Y, X, FILEA, BUF)
      IMPLICIT NONE
      INTEGER FILEA(7), BUF(1)
      REAL    Y(1), X(1)
C
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
      INTEGER IFILE, ISLOT, I, J, N
      REAL    SUM
C
C     Get file number
      IFILE = FILEA(1)
      IF (IFILE .EQ. 0) RETURN
C
C     Map file number to matrix storage slot (1-10)
      CALL MG_FILE2SLOT(IFILE, ISLOT)
      IF (ISLOT .LE. 0 .OR. ISLOT .GT. 10) THEN
        WRITE(*,*) 'MATVEC: Unknown file ', IFILE
        RETURN
      ENDIF
C
      N = NROW_G
C
C     X = X + MAT_STORE(:,:,ISLOT) * Y
      DO I = 1, N
        SUM = 0.0
        DO J = 1, N
          SUM = SUM + MAT_STORE(I,J,ISLOT) * Y(J)
        ENDDO
        X(I) = X(I) + SUM
      ENDDO
C
      RETURN
      END
C
C ============================================================
C     INTFBS - Forward/Backward Substitution (REPLACES NASTRAN's)
C
C     Solves K_eff * DY = DX using pre-computed inverse.
C     In our mini-GINO, we store K_eff_inv in slot for ISCR2.
C     Result: DY = K_eff_inv * DX
C ============================================================
C
      SUBROUTINE INTFBS (DX, DY, IOBUF)
      IMPLICIT NONE
      REAL    DX(1), DY(1)
      INTEGER IOBUF(1)
C
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
C     INFBSX common used by STEP to pass file IDs
      INTEGER FILEL(7), FILEU(7)
      COMMON /INFBSX/ FILEL, FILEU
C
C     We use INVERD to solve: store the effective stiffness
C     in a special slot (slot 9 = "K_eff_inv")
C     Actually simpler: DY = K_eff_inv * DX
C
      INTEGER I, J, N
      REAL    SUM
C
      N = NROW_G
C
C     Slot 9 stores K_eff_inverse
      DO I = 1, N
        SUM = 0.0
        DO J = 1, N
          SUM = SUM + MAT_STORE(I,J,9) * DX(J)
        ENDDO
        DY(I) = SUM
      ENDDO
C
      RETURN
      END
C
C ============================================================
C     FBSINT - Symmetric forward/backward (for ISYM=0 path)
C     Same as INTFBS for our purposes
C ============================================================
C
      SUBROUTINE FBSINT (DX, DY)
      IMPLICIT NONE
      REAL    DX(1), DY(1)
C
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
      INTEGER I, J, N
      REAL    SUM
C
      N = NROW_G
C
C     Slot 9 stores K_eff_inverse
      DO I = 1, N
        SUM = 0.0
        DO J = 1, N
          SUM = SUM + MAT_STORE(I,J,9) * DX(J)
        ENDDO
        DY(I) = SUM
      ENDDO
C
      RETURN
      END
C
C ============================================================
C     GINO FILE I/O STUBS
C     These provide the minimal interface that TRD1C expects
C ============================================================
C
C     OPEN - Open a file for reading/writing
      SUBROUTINE OPEN (*, IFILE, BUF, MODE)
      INTEGER IFILE, BUF(1), MODE
C     In mini-GINO, opening always succeeds
C     Reset file position to start
      INTEGER FILE_POS, FILE_NCOL, MCB_STORE, NROW_G
      REAL    MAT_STORE, VEC_STORE
      COMMON /MINIGINO/
     1  MAT_STORE(3,3,10),
     2  VEC_STORE(3,30),
     3  MCB_STORE(7,30),
     4  FILE_POS(30),
     5  FILE_NCOL(30),
     6  NROW_G
      INTEGER IVSLOT
      CALL MG_FILE2VSLOT(IFILE, IVSLOT)
      IF (IVSLOT .GT. 0 .AND. IVSLOT .LE. 30) THEN
        FILE_POS(IVSLOT) = 1
      ENDIF
      RETURN
      END
C
C     GOPEN - General open (same as OPEN for us)
      SUBROUTINE GOPEN (IFILE, BUF, MODE)
      INTEGER IFILE, BUF(1), MODE
      INTEGER FILE_POS, FILE_NCOL, MCB_STORE, NROW_G
      REAL    MAT_STORE, VEC_STORE
      COMMON /MINIGINO/
     1  MAT_STORE(3,3,10),
     2  VEC_STORE(3,30),
     3  MCB_STORE(7,30),
     4  FILE_POS(30),
     5  FILE_NCOL(30),
     6  NROW_G
      INTEGER IVSLOT
      CALL MG_FILE2VSLOT(IFILE, IVSLOT)
      IF (IVSLOT .GT. 0 .AND. IVSLOT .LE. 30) THEN
        IF (MODE .EQ. 0) FILE_POS(IVSLOT) = 1
C       mode=1: write from start
C       mode=2: append (keep position)
C       mode=3: append
      ENDIF
      RETURN
      END
C
C     CLOSE - Close a file
      SUBROUTINE CLOSE (IFILE, MODE)
      INTEGER IFILE, MODE
C     No-op in mini-GINO
      RETURN
      END
C
C     REWIND - Rewind a file  
      SUBROUTINE REWIND (IFILE)
      INTEGER IFILE
      INTEGER FILE_POS, FILE_NCOL, MCB_STORE, NROW_G
      REAL    MAT_STORE, VEC_STORE
      COMMON /MINIGINO/
     1  MAT_STORE(3,3,10),
     2  VEC_STORE(3,30),
     3  MCB_STORE(7,30),
     4  FILE_POS(30),
     5  FILE_NCOL(30),
     6  NROW_G
      INTEGER IVSLOT
      CALL MG_FILE2VSLOT(IFILE, IVSLOT)
      IF (IVSLOT .GT. 0 .AND. IVSLOT .LE. 30) THEN
        FILE_POS(IVSLOT) = 1
      ENDIF
      RETURN
      END
C
C     PACK - Store a vector as next column of a file
      SUBROUTINE PACK (VEC, IFILE, MCB)
      REAL    VEC(1)
      INTEGER IFILE, MCB(7)
      INTEGER FILE_POS, FILE_NCOL, MCB_STORE, NROW_G
      REAL    MAT_STORE, VEC_STORE
      COMMON /MINIGINO/
     1  MAT_STORE(3,3,10),
     2  VEC_STORE(3,30),
     3  MCB_STORE(7,30),
     4  FILE_POS(30),
     5  FILE_NCOL(30),
     6  NROW_G
      INTEGER IVSLOT, ICOL, I, N
      N = NROW_G
      CALL MG_FILE2VSLOT(IFILE, IVSLOT)
      IF (IVSLOT .LE. 0 .OR. IVSLOT .GT. 30) RETURN
C     Get current column count, increment
      ICOL = FILE_NCOL(IVSLOT) + 1
      FILE_NCOL(IVSLOT) = ICOL
C     Store vector (use offset: (IVSLOT-1)*0 + column index)
C     We use a flat scheme: vector slot = IVSLOT base + column
C     But we have limited storage. Use direct mapping:
C     For UDV (file slot), store columns sequentially
      DO I = 1, N
        VEC_STORE(I, IVSLOT + ICOL - 1) = VEC(I)
      ENDDO
C     Update MCB column count
      MCB(2) = ICOL
      RETURN
      END
C
C     UNPACK - Read next column from a file into vector
      SUBROUTINE UNPACK (*, IFILE, VEC)
      REAL    VEC(1)
      INTEGER IFILE
      INTEGER FILE_POS, FILE_NCOL, MCB_STORE, NROW_G
      REAL    MAT_STORE, VEC_STORE
      COMMON /MINIGINO/
     1  MAT_STORE(3,3,10),
     2  VEC_STORE(3,30),
     3  MCB_STORE(7,30),
     4  FILE_POS(30),
     5  FILE_NCOL(30),
     6  NROW_G
      INTEGER IVSLOT, IPOS, I, N
      N = NROW_G
      CALL MG_FILE2VSLOT(IFILE, IVSLOT)
      IF (IVSLOT .LE. 0 .OR. IVSLOT .GT. 30) THEN
C       Return alternate (EOF)
        RETURN 1
      ENDIF
      IPOS = FILE_POS(IVSLOT)
      IF (IPOS .GT. FILE_NCOL(IVSLOT)) THEN
C       No more columns - return alternate (EOF)
        RETURN 1
      ENDIF
C     Read vector
      DO I = 1, N
        VEC(I) = VEC_STORE(I, IVSLOT + IPOS - 1)
      ENDDO
      FILE_POS(IVSLOT) = IPOS + 1
      RETURN
      END
C
C     FREAD - Read N words from file
      SUBROUTINE FREAD (IFILE, VEC, NWORDS, IEOR)
      INTEGER IFILE, NWORDS, IEOR
      REAL    VEC(1)
      INTEGER FILE_POS, FILE_NCOL, MCB_STORE, NROW_G
      REAL    MAT_STORE, VEC_STORE
      COMMON /MINIGINO/
     1  MAT_STORE(3,3,10),
     2  VEC_STORE(3,30),
     3  MCB_STORE(7,30),
     4  FILE_POS(30),
     5  FILE_NCOL(30),
     6  NROW_G
      INTEGER IVSLOT, IPOS, I, N
      N = NROW_G
      CALL MG_FILE2VSLOT(IFILE, IVSLOT)
      IF (IVSLOT .LE. 0 .OR. IVSLOT .GT. 30) RETURN
      IPOS = FILE_POS(IVSLOT)
C     Read as next record (= next vector column)
      DO I = 1, MIN(NWORDS, N)
        VEC(I) = VEC_STORE(I, IVSLOT + IPOS - 1)
      ENDDO
      FILE_POS(IVSLOT) = IPOS + 1
      RETURN
      END
C
C     READ - Read with EOF/EOR returns (used in restart code)
      SUBROUTINE READ (*, *, IFILE, VEC, NWORDS, IFLG, IFLAG)
      INTEGER IFILE, NWORDS, IFLG, IFLAG
      REAL    VEC(1)
C     For our simple case, just return EOF
      IFLAG = 0
      RETURN 1
      END
C
C     WRITE - Write N words to file
      SUBROUTINE WRITE (IFILE, VEC, NWORDS, IEOR)
      INTEGER IFILE, NWORDS, IEOR
      REAL    VEC(1)
      INTEGER FILE_POS, FILE_NCOL, MCB_STORE, NROW_G
      REAL    MAT_STORE, VEC_STORE
      COMMON /MINIGINO/
     1  MAT_STORE(3,3,10),
     2  VEC_STORE(3,30),
     3  MCB_STORE(7,30),
     4  FILE_POS(30),
     5  FILE_NCOL(30),
     6  NROW_G
      INTEGER IVSLOT, IPOS, I, N
      N = NROW_G
      CALL MG_FILE2VSLOT(IFILE, IVSLOT)
      IF (IVSLOT .LE. 0 .OR. IVSLOT .GT. 30) RETURN
      IPOS = FILE_POS(IVSLOT)
      DO I = 1, MIN(NWORDS, N)
        VEC_STORE(I, IVSLOT + IPOS - 1) = VEC(I)
      ENDDO
      IF (IEOR .NE. 0) THEN
        FILE_NCOL(IVSLOT) = MAX(FILE_NCOL(IVSLOT), IPOS)
        FILE_POS(IVSLOT) = IPOS + 1
      ENDIF
      RETURN
      END
C
C     SKPREC - Skip N records
      SUBROUTINE SKPREC (IFILE, N)
      INTEGER IFILE, N
      INTEGER FILE_POS, FILE_NCOL, MCB_STORE, NROW_G
      REAL    MAT_STORE, VEC_STORE
      COMMON /MINIGINO/
     1  MAT_STORE(3,3,10),
     2  VEC_STORE(3,30),
     3  MCB_STORE(7,30),
     4  FILE_POS(30),
     5  FILE_NCOL(30),
     6  NROW_G
      INTEGER IVSLOT
      CALL MG_FILE2VSLOT(IFILE, IVSLOT)
      IF (IVSLOT .GT. 0 .AND. IVSLOT .LE. 30) THEN
        FILE_POS(IVSLOT) = FILE_POS(IVSLOT) + N
        IF (FILE_POS(IVSLOT) .LT. 1) FILE_POS(IVSLOT) = 1
      ENDIF
      RETURN
      END
C
C     BCKREC - Back up one record (same as SKPREC -1)
      SUBROUTINE BCKREC (IFILE)
      INTEGER IFILE
      CALL SKPREC(IFILE, -1)
      RETURN
      END
C
C     FWDREC - Forward one record (alternate return on EOF)
      SUBROUTINE FWDREC (*, IFILE)
      INTEGER IFILE
C     No-op for mini-GINO (MATVEC is replaced)
      RETURN
      END
C
C     RDTRL - Read trailer (MCB) for a file
      SUBROUTINE RDTRL (MCB)
      INTEGER MCB(7)
      INTEGER FILE_POS, FILE_NCOL, MCB_STORE, NROW_G
      REAL    MAT_STORE, VEC_STORE
      COMMON /MINIGINO/
     1  MAT_STORE(3,3,10),
     2  VEC_STORE(3,30),
     3  MCB_STORE(7,30),
     4  FILE_POS(30),
     5  FILE_NCOL(30),
     6  NROW_G
      INTEGER IFILE, IVSLOT, I
      IFILE = MCB(1)
      CALL MG_FILE2VSLOT(IFILE, IVSLOT)
      IF (IVSLOT .GT. 0 .AND. IVSLOT .LE. 30) THEN
        DO I = 1, 7
          MCB(I) = MCB_STORE(I, IVSLOT)
        ENDDO
      ENDIF
      RETURN
      END
C
C     WRTTRL - Write trailer (MCB) for a file
      SUBROUTINE WRTTRL (MCB)
      INTEGER MCB(7)
      INTEGER FILE_POS, FILE_NCOL, MCB_STORE, NROW_G
      REAL    MAT_STORE, VEC_STORE
      COMMON /MINIGINO/
     1  MAT_STORE(3,3,10),
     2  VEC_STORE(3,30),
     3  MCB_STORE(7,30),
     4  FILE_POS(30),
     5  FILE_NCOL(30),
     6  NROW_G
      INTEGER IFILE, IVSLOT, I
      IFILE = MCB(1)
      CALL MG_FILE2VSLOT(IFILE, IVSLOT)
      IF (IVSLOT .GT. 0 .AND. IVSLOT .LE. 30) THEN
        DO I = 1, 7
          MCB_STORE(I, IVSLOT) = MCB(I)
        ENDDO
      ENDIF
      RETURN
      END
C
C     MAKMCB - Create a Matrix Control Block
      SUBROUTINE MAKMCB (MCB, IFILE, NROW, IFORM, ITYPE)
      INTEGER MCB(7), IFILE, NROW, IFORM, ITYPE
      INTEGER FILE_POS2, FILE_NCOL2, MCB_STORE2, NROW_G2
      REAL    MAT_STORE2, VEC_STORE2
      COMMON /MINIGINO/
     1  MAT_STORE2(3,3,10),
     2  VEC_STORE2(3,30),
     3  MCB_STORE2(7,30),
     4  FILE_POS2(30),
     5  FILE_NCOL2(30),
     6  NROW_G2
      INTEGER IVSLOT
      MCB(1) = IFILE
      MCB(2) = 0
      MCB(3) = NROW
      MCB(4) = IFORM
      MCB(5) = ITYPE
      MCB(6) = 0
      MCB(7) = 0
      CALL MG_FILE2VSLOT(IFILE, IVSLOT)
      IF (IVSLOT .GT. 0 .AND. IVSLOT .LE. 30) THEN
        MCB_STORE2(1, IVSLOT) = IFILE
        MCB_STORE2(2, IVSLOT) = 0
        MCB_STORE2(3, IVSLOT) = NROW
        MCB_STORE2(4, IVSLOT) = IFORM
        MCB_STORE2(5, IVSLOT) = ITYPE
        MCB_STORE2(6, IVSLOT) = 0
        MCB_STORE2(7, IVSLOT) = 0
      ENDIF
      RETURN
      END
C
C     KORSZ - defined as INTEGER FUNCTION in driver program
C     (removed subroutine version to avoid conflict)
C
C     TMTOGO - Return time remaining
      SUBROUTINE TMTOGO (ITLEFT)
      INTEGER ITLEFT
C     Always have time
      ITLEFT = 999
      RETURN
      END
C
C     SSWTCH - Sense switch (always off)
      SUBROUTINE SSWTCH (N, IVAL)
      INTEGER N, IVAL
      IVAL = 0
      RETURN
      END
C
C     MESAGE - Error message handler
      SUBROUTINE MESAGE (IPARM, IFILE, SUBNAM)
      INTEGER IPARM, IFILE
      INTEGER SUBNAM(2)
      WRITE(*,100) IPARM, IFILE
  100 FORMAT(' *** NASTRAN MESAGE: CODE=',I4,' FILE=',I6)
      IF (IPARM .LT. 0) THEN
        WRITE(*,*) 'FATAL ERROR - STOPPING'
        STOP
      ENDIF
      RETURN
      END
C
C     TRD1D - Nonlinear load computation (stub - not used)
      SUBROUTINE TRD1D
C     No-op for linear analysis
      RETURN
      END
C
C     INTPK - Initialize unpacking (stub - not called in our path)
      SUBROUTINE INTPK (*, IFILE, I1, I2, I3)
      INTEGER IFILE, I1, I2, I3
C     Return alternate (EOF) since our MATVEC doesn't use this
      RETURN 1
      END
C
C     ZNTPKI - Unpack next element (stub)
      SUBROUTINE ZNTPKI
C     No-op - not used in our replacement MATVEC
      RETURN
      END
C
C ============================================================
C     FILE NUMBER TO SLOT MAPPING
C
C     Convention:
C       IC   = 101 -> vector slot 1 (cols 1-2: u0, udot0)
C       PD   = 102 -> vector slot 3 (cols: load vectors)
C       UDV  = 201 -> vector slot 5 (cols: output u,v,a)
C       SCR1 = 301 -> matrix slot 1, vector slot 15
C       SCR2 = 302 -> matrix slot 2 (L factor)
C       SCR3 = 303 -> matrix slot 3 (U factor)
C       SCR4 = 304 -> matrix slot 4
C       SCR5 = 305 -> vector slot 17
C       SCR6 = 306 -> vector slot 18
C       ISCR9= 309 -> vector slot 19
C       K_eff_inv  -> matrix slot 9 (used by INTFBS)
C ============================================================
C
      SUBROUTINE MG_FILE2SLOT (IFILE, ISLOT)
C     Map file number to matrix storage slot (1-10)
      INTEGER IFILE, ISLOT
      IF (IFILE .EQ. 301) THEN
        ISLOT = 1
      ELSE IF (IFILE .EQ. 302) THEN
        ISLOT = 2
      ELSE IF (IFILE .EQ. 303) THEN
        ISLOT = 3
      ELSE IF (IFILE .EQ. 304) THEN
        ISLOT = 4
      ELSE IF (IFILE .EQ. 305) THEN
        ISLOT = 5
      ELSE IF (IFILE .EQ. 306) THEN
        ISLOT = 6
      ELSE
C       Unknown file - check if it's K, M, or B matrix
C       IK file = 401, IM = 402, IB = 403
        IF (IFILE .EQ. 401) THEN
          ISLOT = 7
        ELSE IF (IFILE .EQ. 402) THEN
          ISLOT = 8
        ELSE IF (IFILE .EQ. 403) THEN
          ISLOT = 10
        ELSE
          ISLOT = 0
        ENDIF
      ENDIF
      RETURN
      END
C
      SUBROUTINE MG_FILE2VSLOT (IFILE, IVSLOT)
C     Map file number to vector storage slot (1-30)
      INTEGER IFILE, IVSLOT
      IF (IFILE .EQ. 101) THEN
        IVSLOT = 1
      ELSE IF (IFILE .EQ. 102) THEN
        IVSLOT = 3
      ELSE IF (IFILE .EQ. 201) THEN
        IVSLOT = 5
      ELSE IF (IFILE .EQ. 301) THEN
        IVSLOT = 15
      ELSE IF (IFILE .EQ. 302) THEN
        IVSLOT = 16
      ELSE IF (IFILE .EQ. 303) THEN
        IVSLOT = 17
      ELSE IF (IFILE .EQ. 305) THEN
        IVSLOT = 18
      ELSE IF (IFILE .EQ. 306) THEN
        IVSLOT = 19
      ELSE IF (IFILE .EQ. 309) THEN
        IVSLOT = 20
      ELSE IF (IFILE .EQ. 203) THEN
        IVSLOT = 21
      ELSE
        IVSLOT = 0
      ENDIF
      RETURN
      END
