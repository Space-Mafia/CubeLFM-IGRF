C Initialize the existing centered dipole ABI for its exact native functions.
C Infinite-background split: no radius taper and no initial-envelope cutoff.
C The evolved total field is not held fixed by this decomposition.
      SUBROUTINE LFM_CUBE_BACKGROUND_SETUP(MOMENT)
#include "dipole.inc"
      REAL MOMENT
      GEOQMU=MOMENT
      RINNER=1.E30
      ROUTER=2.E30
      A1=1.
      B1=0.
      AQ=0.
      BQ=0.
      XQ_OFF=0.
      BXQ=0.
      SINNER=0.96
      SOUTER=1.
      COS_TILT=1.
      SIN_TILT=0.
      END

C Same analytic Earth reference and initial-potential parameters as the
C existing LFM Earth namelist, without replacing any native field function.
      SUBROUTINE LFM_CUBE_EARTH_BACKGROUND(MOMENT,RP)
#include "dipole.inc"
      REAL MOMENT,RP,AHALF
      CALL LFM_CUBE_BACKGROUND_SETUP(MOMENT)
      RINNER=4.5*RP
      ROUTER=8.*RP
      A1=ROUTER/(ROUTER-RINNER)
      B1=-1./(ROUTER-RINNER)
      AHALF=62.5*RP
      AQ=1./AHALF
      BQ=1./(80.*RP)
      XQ_OFF=25.*RP-AHALF
      END
