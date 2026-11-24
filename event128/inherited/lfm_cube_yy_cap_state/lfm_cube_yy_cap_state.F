C State/snapshot transport only. No solve or dipole INTRP mapping here.
      SUBROUTINE LFM_CUBE_YY_INSTALL_CAPS(PN,PS,CN,CS,STEP,TIME,
     $ SPATIAL,STATUS)
      USE, INTRINSIC :: IEEE_ARITHMETIC
#include "ion-param.inc"
#include "ion_bridge_common.inc"
#include "volts.inc"
      REAL PN(NJP1,NK2P3),PS(NJP1,NK2P3)
      REAL CN(NJP1,NK2P3),CS(NJP1,NK2P3),TIME
      REAL SPATIAL(NJP1,NK2P3,10,2)
      INTEGER STEP,STATUS
      STATUS=1
      IF(.NOT.ALL(IEEE_IS_FINITE(PN)).OR.
     $ .NOT.ALL(IEEE_IS_FINITE(PS)).OR.
     $ .NOT.ALL(IEEE_IS_FINITE(CN)).OR.
     $ .NOT.ALL(IEEE_IS_FINITE(CS))) RETURN
      ION_PSI_N=PN
      ION_PSI_S=PS
      ION_CUR_N=CN
      ION_CUR_S=CS
      PSI=PS
      ION_TIME_BR=TIME
      ION_LSTEP_BR=STEP
      CALL LFM_SW_MI_SNAPSHOT_PACK(NJP1,NK2P3,SPATIAL,STATUS)
      END
