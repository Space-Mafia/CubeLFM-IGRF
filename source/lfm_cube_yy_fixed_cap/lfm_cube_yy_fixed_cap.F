C Driver adapter only: retain original COEF1/PSISLV/Z2ICG discretization.
C RN/RS are already mapped to the donor fixed-cap nodes and carry the
C SAME sign/units as CURNORTH/CURSOUTH, not unsigned FAC magnitude.
C This routine deliberately adds NO corotation and does NOT install edges.
      SUBROUTINE LFM_CUBE_YY_SOLVE_CAP(RN,RS,PN,PS,STATUS)
      USE, INTRINSIC :: IEEE_ARITHMETIC
#include "ion-param.inc"
#include "ion_bridge_common.inc"
#include "volts.inc"
      REAL RN(NJP1,NK2P3),RS(NJP1,NK2P3)
      REAL PN(NJP1,NK2P3),PS(NJP1,NK2P3),ERR
      INTEGER STATUS,H,NSOLVE,CALLS,FAILED,ZERO_RHS
      COMMON /CUBE_YY_EXACT_ZERO_RHS/ ZERO_RHS
      REAL WORST
      COMMON /LFM_MI_SOLVE_CHECK_I/ CALLS,FAILED
      COMMON /LFM_MI_SOLVE_CHECK_R/ WORST
      STATUS=1
C     Only the requested constant-conductance branch is implemented.
      IF(ION_MODEL_BR.OR.ION_TABLE_BR.OR.ION_PED0_BR.LE.0.) RETURN
      IF(ION_HALLMIN_BR.NE.0..OR.ION_SIGRATIO_BR.NE.0.) RETURN
      IF(ION_PHARDY_BR.NE.0..OR.ION_CHARDY_BR.NE.0.) RETURN
      IF(.NOT.ALL(IEEE_IS_FINITE(RN)).OR.
     $   .NOT.ALL(IEEE_IS_FINITE(RS))) RETURN
      IF(.NOT.ALL(IEEE_IS_FINITE(PN)).OR.
     $   .NOT.ALL(IEEE_IS_FINITE(PS))) RETURN
C     Boundary/guard RHS entries must already follow the donor topology.
C     Reject instead of silently deleting incorrectly mapped boundary FAC.
      IF(ANY(RN(1,:).NE.0.).OR.ANY(RN(NJP1,:).NE.0.).OR.
     $   ANY(RS(1,:).NE.0.).OR.ANY(RS(NJP1,:).NE.0.)) RETURN
      IF(ANY(RN(:,1:2).NE.0.).OR.ANY(RN(:,NK2+2:).NE.0.).OR.
     $   ANY(RS(:,1:2).NE.0.).OR.ANY(RS(:,NK2+2:).NE.0.)) RETURN
      CALL LFM_MI_CHECK_RESET
      NSOLVE=0
      DO H=1,2
        IF(H.EQ.1) THEN
          CURNT=RN
          PSI=PN
        ELSE
          CURNT=RS
          PSI=PS
        END IF
C       Always enter original PSISLV so its saved FIRST/warm-start lifecycle
C       is identical on cold initialization and own checkpoint restoration.
C       The wrapper handles only exactly zero homogeneous RHS, not small FAC.
        ZERO_RHS=0
        IF(ALL(CURNT.EQ.0.)) ZERO_RHS=1
        CALL PSISLV(ERR,H.EQ.1)
        ZERO_RHS=0
        CALL LFM_CUBE_YY_RECORD_SOLVE(H)
        NSOLVE=NSOLVE+1
        IF(.NOT.ALL(IEEE_IS_FINITE(PSI))) RETURN
        IF(H.EQ.1) THEN
          PN=PSI
        ELSE
          PS=PSI
        END IF
      END DO
      STATUS=2
      IF(CALLS.EQ.NSOLVE.AND.FAILED.EQ.0) STATUS=0
      END
