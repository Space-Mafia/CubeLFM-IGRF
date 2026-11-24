C Recovered dipole parameter contract used by the analytic background field
C and the original LFM initial-load vector potential.  The historical COMMON
C name/order was lost with dipole.inc, so this remains a named rescue ABI.
      REAL GEOQMU,RINNER,ROUTER,A1,B1
      REAL AQ,BQ,XQ_OFF,BXQ,SINNER,SOUTER,COS_TILT,SIN_TILT
      COMMON /LFM_RESCUE_DIPOLE_R/ GEOQMU,RINNER,ROUTER,A1,B1,
     $ AQ,BQ,XQ_OFF,BXQ,SINNER,SOUTER,COS_TILT,SIN_TILT
      SAVE /LFM_RESCUE_DIPOLE_R/
