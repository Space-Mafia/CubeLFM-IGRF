      SUBROUTINE LFM_SW_MI_SNAPSHOT_PACK(NJ_OUT,NK_OUT,A,STATUS)
      USE, INTRINSIC :: IEEE_ARITHMETIC
      USE LFM_ROTATING_MI_STATE, ONLY: ROTATING_MI
#include "ion-param.inc"
#include "iongeo.inc"
#include "ion_bridge_common.inc"
#include "lfm_sw_state.inc"
      INTEGER NJ_OUT,NK_OUT,STATUS,I,J,H
      REAL A(NJ_OUT,NK_OUT,10,2),PI,THETA,PHI,LFM_PLANET_PHI
      STATUS=1
      IF (.NOT.LFM_IONOSPHERE_ENABLED) RETURN
      IF (NJ_OUT.NE.ION_NJP1.OR.NK_OUT.NE.ION_NK2P3) RETURN
      PI=4.0*ATAN(1.0)
      DO H=1,2
        DO J=1,NK_OUT
          DO I=1,NJ_OUT
            THETA=SQRT(X2ION(I,J)**2+Y2ION(I,J)**2)
            PHI=ATAN2(Y2ION(I,J),X2ION(I,J))
            A(I,J,1,H)=X2ION(I,J)
            A(I,J,2,H)=Y2ION(I,J)
            A(I,J,3,H)=(3-2*H)*(90.0-180.0*THETA/PI)
            A(I,J,4,H)=MODULO(12.0+12.0*PHI/PI,24.0)
            A(I,J,10,H)=LFM_PLANET_PHI(THETA,ION_TIME_BR)
C           In the dynamic path rigid rotation is an independent inductive
C           edge EMF, not this aligned-dipole scalar corotation potential.
            IF(ROTATING_MI) A(I,J,10,H)=0.
          END DO
        END DO
      END DO
      A(:,:,5,1)=ION_PSI_N
      A(:,:,5,2)=ION_PSI_S
      A(:,:,9,:)=A(:,:,5,:)
      A(:,:,5,:)=A(:,:,9,:)+A(:,:,10,:)
      A(:,:,6,1)=ION_CUR_N
      A(:,:,6,2)=ION_CUR_S
      A(:,:,7,1)=ION_SIGP_N
      A(:,:,7,2)=ION_SIGP_S
      A(:,:,8,1)=ION_SIGH_N
      A(:,:,8,2)=ION_SIGH_S
      IF (ALL(IEEE_IS_FINITE(A))) STATUS=0
      END
