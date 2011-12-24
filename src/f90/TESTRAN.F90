!   This is the test program 'testran.f' for TRANSMP.
!
!   David H. Bailey     June 2, 1994
!
!MP+ PRECISION LEVEL 100
!MP+ MIXED MODE FAST
!MP+ OUTPUT PRECISION 56
!MP+ EPSILON 1E-90
!
      PROGRAM TESTRAN
!MP+ IMPLICIT MULTIP REAL (A-H, O-Z)
!MP+ MULTIP INTEGER IA, IB, IC
!MP+ MULTIP REAL A, B
!MP+ MULTIP COMPLEX C, D, E
      PARAMETER (N = 25)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      DOUBLE PRECISION A(N), B(N)
      DOUBLE COMPLEX C, D, E, DPCMPL
!
!   MP parameter definitions.
!
      PARAMETER (DPEPS = 1D-15, DPPIC = 3.141592653589793D0)
!
      EE = EXP (1.D0+0)
      WRITE (6, *) DPPIC, EE
      S = 0.D0
!
!   Loop with subscripted MP variables.
!
      DO 100 I = 1, N
        A(I) = 2 * I + 1
        B(I) = 2.D0 * A(I) * (A(I) + 1.D0)
        S = S + B(I) ** 2
 100  CONTINUE
!
      WRITE (6, *) S
!
!   An expression with mixed MPI and MPR entities.
!
      IA = S
      IB = 262144
      S = (S + 327.25D0) * MOD (IA, 4 * IB)
      WRITE (6, *) S
!
!   A complex square root reference.
!
      E = SQRT (DPCMPL (2.D0 * S, S))
      WRITE (6, *) E
!
!   External and intrinsic MP function references in expressions.
!
      S = DOT (N, A, B)
      T = 2.D0 * SQRT (S) ** 2
      WRITE (6, *) S, T
      S = S / 1048576.D0
      T = S + 2.D0 * LOG (S)
      X = 3 + NINT (T) * 5
      WRITE (6, *) S, T, X
!
!   Deeply nested expressions and function references.
!
      X = (S + (2 * (S - 5) + 3 * (T - 5))) * EXP (COS (LOG (S)))
      WRITE (6, *) X
!
!   A "special" subroutine call (computes both cos and sin of S).
!
      CALL DPCSSN (S, X, Y)
      T = 1.D0 - (X ** 2 + Y ** 2)
!
!   IF-THEN-ELSE construct involving MP variables.
!
      IF (ABS (T) .LT. DPEPS) THEN
        WRITE (6, *) T
      ELSE
        WRITE (6, *) DPEPS
      ENDIF
!
      STOP
      END
!
!   MP function subprogram.
!
!MP+ MULTIP REAL A, B, DOT, S
      FUNCTION DOT (N, A, B)
      DOUBLE PRECISION A(N), B(N), DOT, S
!
      S = 0.D0
!
      DO 100 I = 1, N
        S = S + A(I) * B(I)
 100  CONTINUE
!
      DOT = S
      RETURN
      END
!
!   DP equivalent of special subroutine DPCSSN.
!
      SUBROUTINE DPCSSN (A, X, Y)
      DOUBLE PRECISION A, X, Y
      X = COS (A)
      Y = SIN (A)
      RETURN
      END
!
!   DP equivalent is special function DPCMPL.
!
      FUNCTION DPCMPL (A, B)
      DOUBLE COMPLEX DPCMPL
      DOUBLE PRECISION A, B
      DPCMPL = DCMPLX (A, B)
      RETURN
      END
