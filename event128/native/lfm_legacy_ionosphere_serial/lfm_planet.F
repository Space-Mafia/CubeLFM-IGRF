      SUBROUTINE LFM_PLANET_CONFIGURE(RCM,BEQ,RION,PERIOD,RAMP,
     $ FAC_T,STATUS)
      USE, INTRINSIC :: IEEE_ARITHMETIC
#include "planet.inc"
      REAL RCM,BEQ,RION,PERIOD,RAMP,FAC_T
      INTEGER STATUS
      STATUS=1
      IF (.NOT.ALL(IEEE_IS_FINITE([RCM,BEQ,RION,PERIOD,
     $ RAMP,FAC_T]))) RETURN
      IF (RCM.LE.0.OR.RION.LE.0.OR.PERIOD.LT.0.OR.RAMP.LT.0
     $ .OR.FAC_T.LE.0) RETURN
      PL_RION=RION
C dipole.F uses minus the usual +Z magnetic moment.
      PL_MSI=-BEQ*RCM**3*1.E-10
      PL_OMEGA=0.
      IF (PERIOD.GT.0) PL_OMEGA=8.*ATAN(1.)/PERIOD
      PL_RAMP=RAMP
C Convert J dot B to the donor radial FAC convention with reversed B.
      PL_FAC=4376.76*(FAC_T/5.5E-5)*SIGN(1.,BEQ)
      STATUS=0
      END

      REAL FUNCTION LFM_PLANET_RION()
#include "planet.inc"
      LFM_PLANET_RION=PL_RION
      END

      REAL FUNCTION LFM_PLANET_FAC()
#include "planet.inc"
      LFM_PLANET_FAC=PL_FAC
      END

      REAL FUNCTION LFM_PLANET_PHI(THETA,TIME)
#include "ion-param.inc"
#include "iongeo.inc"
#include "planet.inc"
      REAL THETA,TIME,W
      W=1.
      IF (PL_RAMP.GT.0.) THEN
        W=MAX(0.,MIN(1.,TIME/PL_RAMP))
        W=W**3*(10.-15.*W+6.*W**2)
      ENDIF
C Gauge zero at spherical MHD inner-equator rim: INTRP enforces zero there.
      LFM_PLANET_PHI=W*PL_OMEGA*PL_MSI/PL_RION*
     $ (SIN(THETA)**2-SIN(TETA1(1,1))**2)
      END

      SUBROUTINE LFM_PLANET_MAP_COROTATION
#include "ion-param.inc"
#include "iongeo.inc"
#include "volts.inc"
#include "ion_bridge_common.inc"
#include "planet.inc"
      INTEGER J,K
      REAL LFM_PLANET_PHI,R(3),B(3),E(3),V(3),RSQ,RR,S2,W
C FAC PSI and its ordinary INTRP results are not replaced. The analytic
C contribution is constant on dipole L and evaluated at target footprints.
C Original IONOSF below this hook forms EJ/EK from total PSII differences.
      IF (PL_OMEGA.EQ.0.) RETURN
      DO K=1,NK2P1
        DO J=1,NJP1
          PSII(J,K)=PSII(J,K)+
     $      LFM_PLANET_PHI(TETA1(J,K),ION_TIME_BR)
        END DO
      END DO
C Keep diagnostic drift consistent: differentiate that SAME scalar
C analytically and calculate E cross B at original VCOEF sample locations.
C No velocity is imposed on MHD; MPI installs EJ/EK only.
      W=1.
      IF (PL_RAMP.GT.0.) THEN
        W=MAX(0.,MIN(1.,ION_TIME_BR/PL_RAMP))
        W=W**3*(10.-15.*W+6.*W**2)
      ENDIF
      DO K=1,NK2
        DO J=1,NJ
          R(1)=SUM(ION_X(1:2,J:J+1,K:K+1))*1.E-2/8.
          R(2)=SUM(ION_Y(1:2,J:J+1,K:K+1))*1.E-2/8.
          R(3)=SUM(ION_Z(1:2,J:J+1,K:K+1))*1.E-2/8.
          RSQ=SUM(R**2)
          RR=SQRT(RSQ)
          S2=(R(1)**2+R(2)**2)/RSQ
          B=PL_MSI/RR**5*
     $      [3.*R(1)*R(3),3.*R(2)*R(3),3.*R(3)**2-RSQ]
          E=W*PL_OMEGA*PL_MSI/RR**3*
     $      [(3.*S2-2.)*R(1),(3.*S2-2.)*R(2),3.*S2*R(3)]
          V=100.*[E(2)*B(3)-E(3)*B(2),
     $       E(3)*B(1)-E(1)*B(3),E(1)*B(2)-E(2)*B(1)]/SUM(B**2)
          VXINT(J,K)=VXINT(J,K)+V(1)
          VYINT(J,K)=VYINT(J,K)+V(2)
          VZINT(J,K)=VZINT(J,K)+V(3)
        END DO
      END DO
      END

      BLOCK DATA LFM_PLANET_DEFAULTS
#include "planet.inc"
#include "prefill.inc"
      DATA PL_RION/6.5E6/,PL_MSI/0./,PL_OMEGA/0./,
     $ PL_RAMP/0./,PL_FAC/4376.76/
      DATA PF_ENABLED/.FALSE./,PF_X/0./,PF_WIDTH/1./,
     $ PF_BY/0./,PF_BZ/0./
      END

      SUBROUTINE LFM_PREFILL_CONFIGURE(ENABLED,RCM,XRJ,WRJ,
     $ BX,BY,BZ,STATUS)
      USE, INTRINSIC :: IEEE_ARITHMETIC
#include "prefill.inc"
      LOGICAL ENABLED
      REAL RCM,XRJ,WRJ,BX,BY,BZ
      INTEGER STATUS
      STATUS=1
      IF(.NOT.ALL(IEEE_IS_FINITE([RCM,XRJ,WRJ,BX,BY,BZ]))) RETURN
      IF(ENABLED.AND.(WRJ.LE.0..OR.BX.NE.0.)) RETURN
      PF_ENABLED=ENABLED
      PF_X=XRJ*RCM
      PF_WIDTH=MAX(WRJ*RCM,1.)
      PF_BY=BY
      PF_BZ=BZ
      STATUS=0
      END

      REAL FUNCTION LFM_PREFILL_WEIGHT(X)
#include "prefill.inc"
      REAL X,T
      LFM_PREFILL_WEIGHT=0.
      IF(.NOT.PF_ENABLED) RETURN
      T=MAX(0.,MIN(1.,(X-PF_X)/PF_WIDTH))
      LFM_PREFILL_WEIGHT=T**3*(10.-15.*T+6.*T*T)
      END

      REAL FUNCTION LFM_PREFILL_PRIMITIVE(X)
#include "prefill.inc"
      REAL X,T
      LFM_PREFILL_PRIMITIVE=0.
      IF(.NOT.PF_ENABLED.OR.X.LE.PF_X) RETURN
      T=(X-PF_X)/PF_WIDTH
      IF(T.GE.1.) THEN
        LFM_PREFILL_PRIMITIVE=X-PF_X-0.5*PF_WIDTH
      ELSE
        LFM_PREFILL_PRIMITIVE=PF_WIDTH*(2.5*T**4-3.*T**5+T**6)
      ENDIF
      END
