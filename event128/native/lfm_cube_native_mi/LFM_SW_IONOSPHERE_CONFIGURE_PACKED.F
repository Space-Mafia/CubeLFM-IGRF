      SUBROUTINE LFM_SW_IONOSPHERE_CONFIGURE_PACKED(PEDERSEN_MHO,
     $ GEOMETRY,COUNT,STATUS)
      USE, INTRINSIC :: IEEE_ARITHMETIC
#include "lfm_sw_state.inc"
#include "ion_bridge_common.inc"
      INTEGER COUNT,STATUS,I,J,K,N,FIELD
      REAL PEDERSEN_MHO,GEOMETRY(COUNT),LFM_PLANET_RION

      STATUS=1
      LFM_IONOSPHERE_ENABLED=.FALSE.
      IF (.NOT.LFM_DIPOLE_BACKGROUND) THEN
        WRITE(*,'(A)') 'ionosphere requires the dipole magnetosphere mode'
        RETURN
      END IF
      IF (COUNT .NE. 3*ION_NIP1*ION_NJP1*ION_NKP1) THEN
        WRITE(*,'(A,2(1X,I0))') 'packed ionosphere geometry mismatch:',
     $    COUNT,3*ION_NIP1*ION_NJP1*ION_NKP1
        RETURN
      END IF
      IF (.NOT.IEEE_IS_FINITE(PEDERSEN_MHO) .OR.
     $    PEDERSEN_MHO .LE. 0.0) THEN
        WRITE(*,'(A)') 'constant Pedersen conductance must be positive'
        RETURN
      END IF

      N=0
      DO FIELD=1,3
        DO K=1,ION_NKP1
          DO J=1,ION_NJP1
            DO I=1,ION_NIP1
              N=N+1
              IF (FIELD .EQ. 1) ION_X(I,J,K)=GEOMETRY(N)
              IF (FIELD .EQ. 2) ION_Y(I,J,K)=GEOMETRY(N)
              IF (FIELD .EQ. 3) ION_Z(I,J,K)=GEOMETRY(N)
            END DO
          END DO
        END DO
      END DO

      ION_PE_BR=.TRUE.
      ION_NI_BR=ION_NJP1
      ION_NJ_BR=ION_NK2P3
      ION_MI_BR=ION_NI_BR
      ION_MJ_BR=ION_NJ_BR
      ION_TABLE_BR=.FALSE.
      ION_SM_BR=.FALSE.
      ION_MODEL_BR=.FALSE.
      ION_ITABLE_BR=1
      ION_ITMAX_BR=8000
      ION_LOOPTRIES_BR=1
      ION_ITER_BR=0
      ION_NLOOP_BR=1
      ION_TILT_BR=0.0
      ION_TZERO_BR=0.0
      ION_DELTAT_BR=1.0
      ION_ZANGLE_BR=0.0
      ION_S107_BR=1.0
      ION_PED0_BR=PEDERSEN_MHO
      ION_PEDMIN_BR=PEDERSEN_MHO
      ION_HALLMIN_BR=0.0
      ION_SIGRATIO_BR=0.0
      ION_PHARDY_BR=0.0
      ION_CHARDY_BR=0.0
      ION_HHARDY_BR=0.0
      ION_ERRMAX_BR=1.0E-3
      ION_TIME_BR=0.0
      ION_DT_BR=0.0
      ION_LSTEP_BR=0
      ION_RIONM_BR=LFM_PLANET_RION()
      ION_BX=0.0
      ION_BY=0.0
      ION_BZ=0.0
      ION_RHO=0.0
      ION_CS=0.0
      ION_EJ=0.0
      ION_EK=0.0
      ION_VX=0.0
      ION_VY=0.0
      ION_VZ=0.0
      ION_PSI_N=1.0E-4
      ION_PSI_S=1.0E-4
      ION_CUR_N=0.0
      ION_CUR_S=0.0
      ION_SIGP_N=PEDERSEN_MHO
      ION_SIGH_N=0.0
      ION_SIGP_S=PEDERSEN_MHO
      ION_SIGH_S=0.0
      ION_CURRENT=0.0

      OPEN(UNIT=9,FILE='lfm_ionosphere_legacy.log',STATUS='REPLACE',
     $ ACTION='WRITE')
      CALL LFM_MI_INTERPOLATION_RESET
      CALL IONINIT
      FLUSH(9)
      LFM_IONOSPHERE_ENABLED=.TRUE.
      STATUS=0
      WRITE(*,'(A,F8.3,A)') 'original LFM ionosphere initialized: SigmaP=',
     $ PEDERSEN_MHO,' mho, SigmaH=0 mho'
      RETURN
      END
