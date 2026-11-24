C Retain donor INTRP's inner-sphere nodal potential and zero-node closure.
C Export a COMPLETE sphere in the fixed magnetic frame. Corotation is NOT
C added here: its existing time-averaged Cartesian edge EMF remains separate.
      SUBROUTINE LFM_MI_INNER_TABLE(PN,PS,TABLE,STATUS)
      USE, INTRINSIC :: IEEE_ARITHMETIC
#include "ion-param.inc"
#include "volts.inc"
      REAL PN(NJP1,NK2P3),PS(NJP1,NK2P3),TABLE(NJP1,NKP1)
      INTEGER STATUS,J,K
      STATUS=1
      PSI=PN
      CALL INTRP
      TABLE(:,1:NK2P1)=PSII
      PSI=PS
      CALL INTRP
      DO K=1,NK2P1
      DO J=1,NJP1
        TABLE(J,NKP2-K)=PSII(J,K)
      END DO
      END DO
      IF(.NOT.ALL(IEEE_IS_FINITE(TABLE))) RETURN
      IF(ANY(TABLE(:,1).NE.0.).OR.ANY(TABLE(:,NK2P1).NE.0.))
     $ RETURN
      IF(ANY(TABLE(:,NKP1).NE.TABLE(:,1))) RETURN
      STATUS=0
      END
