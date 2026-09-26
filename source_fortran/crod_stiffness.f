C     ==================================================================
C     NASTRAN-95 STYLE CROD ELEMENT STIFFNESS MATRIX COMPUTATION
C     ==================================================================
C     This is a representative reconstruction of the CROD element
C     stiffness computation as it would appear in NASTRAN-95.
C     The actual NASTRAN-95 code is available at:
C     https://github.com/nasa/NASTRAN-95 (mis/ folder)
C     ==================================================================
C
      SUBROUTINE RODSTF(EK,ECPT,MATBUF)
C
C     COMPUTE ELEMENT STIFFNESS MATRIX FOR CROD ELEMENT
C
C     EK     - OUTPUT ELEMENT STIFFNESS MATRIX (2X2 IN LOCAL)
C              OR (6X6 IN GLOBAL FOR 3D)
C     ECPT   - ELEMENT CONNECTION AND PROPERTY TABLE
C     MATBUF - MATERIAL PROPERTY BUFFER
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION EK(6,6), ECPT(20), MATBUF(10)
C
C     ECPT CONTENTS FOR CROD:
C       ECPT(1)  = ELEMENT ID
C       ECPT(2)  = GRID POINT A (NODE 1)
C       ECPT(3)  = GRID POINT B (NODE 2)
C       ECPT(4)  = PROPERTY ID
C       ECPT(5)  = X1A (X COORD OF NODE A)
C       ECPT(6)  = Y1A (Y COORD OF NODE A)
C       ECPT(7)  = Z1A (Z COORD OF NODE A)
C       ECPT(8)  = X1B (X COORD OF NODE B)
C       ECPT(9)  = Y1B (Y COORD OF NODE B)
C       ECPT(10) = Z1B (Z COORD OF NODE B)
C
C     MATBUF CONTENTS:
C       MATBUF(1) = YOUNG'S MODULUS (E)
C       MATBUF(2) = CROSS-SECTIONAL AREA (A)
C       MATBUF(3) = TORSIONAL CONSTANT (J) - NOT USED FOR CROD
C       MATBUF(4) = NON-STRUCTURAL MASS (NSM)
C
C     --- EXTRACT COORDINATES ---
      X1 = ECPT(5)
      Y1 = ECPT(6)
      Z1 = ECPT(7)
      X2 = ECPT(8)
      Y2 = ECPT(9)
      Z2 = ECPT(10)
C
C     --- COMPUTE ELEMENT LENGTH ---
      DX = X2 - X1
      DY = Y2 - Y1
      DZ = Z2 - Z1
      ELEN = DSQRT(DX*DX + DY*DY + DZ*DZ)
C
C     --- CHECK FOR ZERO LENGTH ---
      IF (ELEN .LT. 1.0D-10) THEN
          WRITE(6,*) 'FATAL ERROR: ZERO LENGTH ROD ELEMENT'
          RETURN
      ENDIF
C
C     --- EXTRACT MATERIAL PROPERTIES ---
      E = MATBUF(1)
      A = MATBUF(2)
C
C     --- COMPUTE AXIAL STIFFNESS ---
      AK = E * A / ELEN
C
C     --- DIRECTION COSINES ---
      CX = DX / ELEN
      CY = DY / ELEN
      CZ = DZ / ELEN
C
C     --- INITIALIZE STIFFNESS MATRIX TO ZERO ---
      DO 10 I = 1, 6
          DO 10 J = 1, 6
              EK(I,J) = 0.0D0
   10 CONTINUE
C
C     --- FORM 6X6 STIFFNESS IN GLOBAL COORDINATES ---
C     K = AK * [T]^T * [k_local] * [T]
C     WHERE k_local = [[1,-1],[-1,1]] (in axial direction)
C     AND T = [cx, cy, cz] (direction cosine vector)
C
C     UPPER LEFT 3X3 (NODE A - NODE A)
      EK(1,1) = AK * CX * CX
      EK(1,2) = AK * CX * CY
      EK(1,3) = AK * CX * CZ
      EK(2,1) = EK(1,2)
      EK(2,2) = AK * CY * CY
      EK(2,3) = AK * CY * CZ
      EK(3,1) = EK(1,3)
      EK(3,2) = EK(2,3)
      EK(3,3) = AK * CZ * CZ
C
C     UPPER RIGHT 3X3 (NODE A - NODE B) = -[UPPER LEFT]
      DO 20 I = 1, 3
          DO 20 J = 1, 3
              EK(I,J+3) = -EK(I,J)
   20 CONTINUE
C
C     LOWER LEFT 3X3 (NODE B - NODE A) = -[UPPER LEFT]
      DO 30 I = 1, 3
          DO 30 J = 1, 3
              EK(I+3,J) = -EK(I,J)
   30 CONTINUE
C
C     LOWER RIGHT 3X3 (NODE B - NODE B) = [UPPER LEFT]
      DO 40 I = 1, 3
          DO 40 J = 1, 3
              EK(I+3,J+3) = EK(I,J)
   40 CONTINUE
C
      RETURN
      END
C
C     ==================================================================
      SUBROUTINE RODSTR(STRESS,DISP,ECPT,MATBUF)
C     ==================================================================
C
C     COMPUTE STRESS IN CROD ELEMENT
C
C     STRESS - OUTPUT AXIAL STRESS
C     DISP   - NODAL DISPLACEMENTS (6 DOF: UXA,UYA,UZA,UXB,UYB,UZB)
C     ECPT   - ELEMENT CONNECTION AND PROPERTY TABLE
C     MATBUF - MATERIAL PROPERTY BUFFER
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION DISP(6), ECPT(20), MATBUF(10)
C
C     --- EXTRACT COORDINATES ---
      X1 = ECPT(5)
      Y1 = ECPT(6)
      Z1 = ECPT(7)
      X2 = ECPT(8)
      Y2 = ECPT(9)
      Z2 = ECPT(10)
C
C     --- COMPUTE LENGTH AND DIRECTION COSINES ---
      DX = X2 - X1
      DY = Y2 - Y1
      DZ = Z2 - Z1
      ELEN = DSQRT(DX*DX + DY*DY + DZ*DZ)
C
      CX = DX / ELEN
      CY = DY / ELEN
      CZ = DZ / ELEN
C
C     --- EXTRACT MATERIAL PROPERTIES ---
      E = MATBUF(1)
C
C     --- COMPUTE AXIAL DISPLACEMENT DIFFERENCE ---
C     PROJECT GLOBAL DISPLACEMENTS ONTO ELEMENT AXIS
      UA = CX*DISP(1) + CY*DISP(2) + CZ*DISP(3)
      UB = CX*DISP(4) + CY*DISP(5) + CZ*DISP(6)
C
C     --- COMPUTE STRAIN ---
      STRAIN = (UB - UA) / ELEN
C
C     --- COMPUTE STRESS (HOOKE'S LAW) ---
      STRESS = E * STRAIN
C
      RETURN
      END
