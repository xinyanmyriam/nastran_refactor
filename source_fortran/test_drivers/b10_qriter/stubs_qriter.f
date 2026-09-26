C     ================================================================
C     STUBS FOR QRITER TEST DRIVER (B10)
C     These stub routines replace NASTRAN runtime dependencies
C     so that QRITER can be tested standalone.
C     ================================================================
C
      SUBROUTINE GOPEN (IFILE, BUF, MODE)
C     Stub: no-op (QRITER uses this to open output file)
      INTEGER IFILE, MODE
      DOUBLE PRECISION BUF(1)
      RETURN
      END
C
      SUBROUTINE CLOSE (IFILE, MODE)
C     Stub: no-op
      INTEGER IFILE, MODE
      RETURN
      END
C
      SUBROUTINE WRITE (IFILE, DATA, NWDS, EOR)
C     Stub: no-op (QRITER writes eigenvalues to file)
      INTEGER IFILE, NWDS, EOR
      RETURN
      END
C
      SUBROUTINE WRTTRL (MSG)
C     Stub: no-op (writes trailer record)
      INTEGER MSG(1)
      RETURN
      END
C
      SUBROUTINE MESAGE (MSGID, VAL, NITER)
C     Stub: print warning about convergence failure
      INTEGER MSGID, NITER
      DOUBLE PRECISION VAL
      WRITE(*,*) 'WARNING: QRITER convergence issue at value=',VAL,
     1           ' after', NITER, ' iterations'
      RETURN
      END
C
      INTEGER FUNCTION KORSZ (A)
C     Stub: return a large buffer size
C     QRITER uses KORSZ to determine available memory for buffer
      DOUBLE PRECISION A(1)
      KORSZ = 10000
      RETURN
      END
