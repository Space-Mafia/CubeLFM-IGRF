C Pass through the unchanged donor solver, recording its actual residual.
      SUBROUTINE LFM_MI_Z2ICG(FI,A1,A2,A3,A4,A5,A6,A7,A8,A9,SI,
     $ NS5,NX,NY,NDX,NDY,LBOUN,ITMAX,ITER,ERROR,SER,AB1,S1,AB2,S2)
      USE, INTRINSIC :: IEEE_ARITHMETIC
      INTEGER NS5,NX,NY,NDX,NDY,LBOUN(4),ITMAX,ITER
      REAL FI(NDX,NDY),A1(NDX,NDY),A2(NDX,NDY),A3(NDX,NDY)
      REAL A4(NDX,NDY),A5(NDX,NDY),A6(NDX,NDY),A7(NDX,NDY)
      REAL A8(NDX,NDY),A9(NDX,NDY),SI(NDX,NDY),ERROR,SER
      REAL AB1(NDY),AB2(NDY),S1,S2
      INTEGER CALLS,FAILED,ZERO_RHS
      COMMON /CUBE_YY_EXACT_ZERO_RHS/ ZERO_RHS
      REAL WORST,LASTSER
      INTEGER LASTITER,LASTMAX
      COMMON /CUBE_YY_SOLVE_LAST_R/ LASTSER
      COMMON /CUBE_YY_SOLVE_LAST_I/ LASTITER,LASTMAX
      COMMON /LFM_MI_SOLVE_CHECK_I/ CALLS,FAILED
      COMMON /LFM_MI_SOLVE_CHECK_R/ WORST
      IF(ZERO_RHS.EQ.1) THEN
C       Exact homogeneous solution avoids the original Krylov RS=0 -> 0/0.
C       Original PSISLV has still formed/validated the physical coefficients
C       and consumed its FIRST initialization. Never threshold a small RHS.
        IF(NX.LT.3.OR.NY.LT.3.OR.NX.GT.NDX.OR.NY.GT.NDY)
     $    ERROR STOP 'invalid homogeneous native MI dimensions'
        IF(NS5.NE.0.OR.ANY(LBOUN.NE.(/-1,-1,0,0/)))
     $    ERROR STOP 'unsupported homogeneous native MI topology'
        IF(.NOT.ALL(IEEE_IS_FINITE(FI)).OR.
     $     .NOT.ALL(IEEE_IS_FINITE(SI)).OR.
     $     .NOT.ALL(IEEE_IS_FINITE(A1)).OR.
     $     .NOT.ALL(IEEE_IS_FINITE(A2)).OR.
     $     .NOT.ALL(IEEE_IS_FINITE(A3)).OR.
     $     .NOT.ALL(IEEE_IS_FINITE(A4)).OR.
     $     .NOT.ALL(IEEE_IS_FINITE(A5)).OR.
     $     .NOT.ALL(IEEE_IS_FINITE(A6)).OR.
     $     .NOT.ALL(IEEE_IS_FINITE(A7)).OR.
     $     .NOT.ALL(IEEE_IS_FINITE(A8)).OR.
     $     .NOT.ALL(IEEE_IS_FINITE(A9)).OR.
     $     .NOT.ALL(IEEE_IS_FINITE(AB1(2:NY-1))).OR.
     $     .NOT.ALL(IEEE_IS_FINITE(AB2(2:NY-1))).OR.
     $     .NOT.IEEE_IS_FINITE(S1).OR..NOT.IEEE_IS_FINITE(S2))
     $    ERROR STOP 'nonfinite homogeneous native MI input'
        IF(ANY(SI.NE.0.).OR.S1.NE.0..OR.S2.NE.0.)
     $    ERROR STOP 'homogeneous native MI source mismatch'
        FI=0.
        ITER=0
        SER=0.
      ELSE IF(ZERO_RHS.EQ.0) THEN
      CALL Z2ICG(FI,A1,A2,A3,A4,A5,A6,A7,A8,A9,SI,
     $ NS5,NX,NY,NDX,NDY,LBOUN,ITMAX,ITER,ERROR,SER,AB1,S1,AB2,S2)
      ELSE
        ERROR STOP 'invalid exact zero RHS lifecycle'
      ENDIF
      LASTSER=SER
      LASTITER=ITER
      LASTMAX=ITMAX
      CALLS=CALLS+1
      IF(.NOT.IEEE_IS_FINITE(SER)) THEN
        FAILED=1
      ELSE
        WORST=MAX(WORST,SER)
        IF(SER.GT.ERROR) THEN
          WRITE(9,'(A,ES24.16,A,ES24.16)')
     $     'MI_SOLVE_NOT_CONVERGED residual=',SER,
     $     ' requested=',ERROR
        ENDIF
      END IF
      WRITE(9,'(A,I0,A,I0,A,ES16.8,A,ES16.8)')
     $ 'MI_SOLVE call=',CALLS,' iterations=',ITER,
     $ ' residual=',SER,' tolerance=',ERROR
      END

      SUBROUTINE LFM_MI_CHECK_RESET
      INTEGER ZERO_RHS
      COMMON /CUBE_YY_EXACT_ZERO_RHS/ ZERO_RHS
      INTEGER CALLS,FAILED
      REAL WORST
      COMMON /LFM_MI_SOLVE_CHECK_I/ CALLS,FAILED
      COMMON /LFM_MI_SOLVE_CHECK_R/ WORST
      ZERO_RHS=0
      CALLS=0
      FAILED=0
      WORST=0.0
      END

      SUBROUTINE LFM_MI_CHECK_RESULT(STATUS)
      INTEGER STATUS,CALLS,FAILED
      REAL WORST
      COMMON /LFM_MI_SOLVE_CHECK_I/ CALLS,FAILED
      COMMON /LFM_MI_SOLVE_CHECK_R/ WORST
      STATUS=0
      IF(CALLS.NE.2.OR.FAILED.NE.0) THEN
        STATUS=2
        WRITE(*,*) 'MI rejected: calls/failure/residual=',
     $    CALLS,FAILED,WORST
      END IF
      END

      SUBROUTINE LFM_MI_INTERPOLATION_RESET
#include "ion-param.inc"
#include "curint.inc"
      CURCOEF=0.0
      XINTERP=0.0
      YINTERP=0.0
      LMIN=0
      LMAX=0
      LCINT=0
      END
