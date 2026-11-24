C Thin ABI adapter only. Native IONINIT/IONOSF/PSISLV/Z2ICG are linked
C from the separately qualified original Earth build without solver edits.
      SUBROUTINE LFM_CUBE_MI_SETUP(GEOMETRY,COUNT,PEDERSEN,STATUS)
      USE LFM_ROTATING_MI_STATE, ONLY: ROTATING_MI
      USE, INTRINSIC :: IEEE_ARITHMETIC, ONLY: IEEE_IS_FINITE
#include "lfm_sw_state.inc"
#include "dipole.inc"
      INTEGER COUNT,STATUS
      REAL GEOMETRY(COUNT),PEDERSEN
C VCOEF consumes the existing native GEOQMU COMMON. The caller must have
C initialized the same Earth background as the MHD driver; do not reset it.
      STATUS=11
      IF(.NOT.IEEE_IS_FINITE(GEOQMU)) THEN
        WRITE(*,*) 'Cube MI requires finite initialized native GEOQMU'
        RETURN
      ENDIF
      IF(ABS(GEOQMU-.30*(6.38E8)**3).GT.
     $ 1.E-12*ABS(.30*(6.38E8)**3)) THEN
        WRITE(*,*) 'Cube MI requires original Earth GEOQMU setup first'
        RETURN
      ENDIF
      CALL LFM_PLANET_CONFIGURE(6.38E8,.30,6.5E6,0.,0.,
     $ 5.5E-5,STATUS)
      IF(STATUS.NE.0) RETURN
      ROTATING_MI=.FALSE.
      LFM_DIPOLE_BACKGROUND=.TRUE.
      CALL LFM_SW_IONOSPHERE_CONFIGURE_PACKED(PEDERSEN,
     $ GEOMETRY,COUNT,STATUS)
      END

      SUBROUTINE LFM_CUBE_MI_GET(PN,PS,TAB,WORST,STATUS)
#include "ion-param.inc"
#include "ion_bridge_common.inc"
      REAL PN(NJP1,NK2P3),PS(NJP1,NK2P3),TAB(NJP1,NKP1)
      REAL WORST,CHECK_WORST
      INTEGER STATUS,CALLS,FAILED
      COMMON /LFM_MI_SOLVE_CHECK_I/ CALLS,FAILED
      COMMON /LFM_MI_SOLVE_CHECK_R/ CHECK_WORST
      PN=ION_PSI_N
      PS=ION_PSI_S
      WORST=CHECK_WORST
      CALL LFM_MI_INNER_TABLE(PN,PS,TAB,STATUS)
      END

C Import only after the same cold-initial MI solve used by the native
C restart lifecycle. That primes PSISLV's private SAVE FIRST branch.
      SUBROUTINE LFM_CUBE_MI_PUT(PN,PS,STATUS)
      USE, INTRINSIC :: IEEE_ARITHMETIC
#include "ion-param.inc"
#include "ion_bridge_common.inc"
#include "volts.inc"
      REAL PN(NJP1,NK2P3),PS(NJP1,NK2P3)
      INTEGER STATUS
      STATUS=1
      IF(.NOT.ALL(IEEE_IS_FINITE(PN)).OR.
     $ .NOT.ALL(IEEE_IS_FINITE(PS))) RETURN
      ION_PSI_N=PN
      ION_PSI_S=PS
      PSI=PS
      STATUS=0
      END

C Observational configuration/result scalars for bounded direct probes.
      SUBROUTINE LFM_CUBE_MI_METRICS(CONFIG,COUNTS)
#include "ion_bridge_common.inc"
      REAL CONFIG(5),WORST
      INTEGER COUNTS(3),CALLS,FAILED
      COMMON /LFM_MI_SOLVE_CHECK_I/ CALLS,FAILED
      COMMON /LFM_MI_SOLVE_CHECK_R/ WORST
      CONFIG=[ION_PED0_BR,ION_HALLMIN_BR,ION_ERRMAX_BR,
     $ ION_RIONM_BR,REAL(ION_ITMAX_BR)]
      COUNTS=[CALLS,FAILED,ION_ITER_BR]
      END

C IO reconstruction only: reuse existing INTRP and native SNAPSHOT_PACK.
C No IONOSF/PSISLV call and no dependency on science output files.
      SUBROUTINE LFM_CUBE_MI_RESTORE_HELD(PN,PS,CN,CS,STEP,TIME,
     $ TAB,SPATIAL,STATUS)
#include "ion-param.inc"
#include "ion_bridge_common.inc"
      REAL PN(NJP1,NK2P3),PS(NJP1,NK2P3)
      REAL CN(NJP1,NK2P3),CS(NJP1,NK2P3),TIME
      REAL TAB(NJP1,NKP1),SPATIAL(NJP1,NK2P3,10,2)
      INTEGER STEP,STATUS
      CALL LFM_CUBE_MI_PUT(PN,PS,STATUS)
      IF(STATUS.NE.0) RETURN
      ION_CUR_N=CN
      ION_CUR_S=CS
      ION_TIME_BR=TIME
      ION_LSTEP_BR=STEP
C Constant conductances are set by the unchanged configure/init path.
      CALL LFM_MI_INNER_TABLE(PN,PS,TAB,STATUS)
      IF(STATUS.NE.0) RETURN
      CALL LFM_SW_MI_SNAPSHOT_PACK(NJP1,NK2P3,SPATIAL,STATUS)
      END
