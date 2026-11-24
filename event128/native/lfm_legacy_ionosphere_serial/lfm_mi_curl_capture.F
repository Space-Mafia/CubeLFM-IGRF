C Observational adapter for the UNMODIFIED JPARA Cartesian curl.
C Arrays contain curl(B) in G/cm, NOT amperes/m^2. Convert only on export.
C No static dipole projection, hemisphere split, smoothing, or B ratio here.
      SUBROUTINE LFM_MI_CAPTURE_CURL(JX,JY,JZ)
#include "ion-param.inc"
      REAL JX(NJP1,NKP1),JY(NJP1,NKP1),JZ(NJP1,NKP1)
      REAL RAW(3,NJP1,NKP1)
      INTEGER VALID,J,K
      COMMON /LFM_MI_CURL_RAW/ RAW
      COMMON /LFM_MI_CURL_VALID/ VALID
      SAVE /LFM_MI_CURL_RAW/,/LFM_MI_CURL_VALID/
      DO K=1,NKP1
      DO J=1,NJP1
        RAW(1,J,K)=JX(J,K)
        RAW(2,J,K)=JY(J,K)
        RAW(3,J,K)=JZ(J,K)
      END DO
      END DO
      VALID=1
      END

      SUBROUTINE LFM_MI_CURL_RESET
      INTEGER VALID
      COMMON /LFM_MI_CURL_VALID/ VALID
      SAVE /LFM_MI_CURL_VALID/
      VALID=0
      END

      SUBROUTINE LFM_MI_CURL_EXPORT_SI(OUT,COUNT,STATUS)
      USE, INTRINSIC :: IEEE_ARITHMETIC
#include "ion-param.inc"
      INTEGER COUNT,STATUS,VALID,J,K,N,F
      REAL OUT(COUNT),RAW(3,NJP1,NKP1),FACTOR
      COMMON /LFM_MI_CURL_RAW/ RAW
      COMMON /LFM_MI_CURL_VALID/ VALID
      SAVE /LFM_MI_CURL_RAW/,/LFM_MI_CURL_VALID/
      STATUS=1
      IF(COUNT.NE.3*NJP1*NKP1.OR.VALID.NE.1) RETURN
      IF(.NOT.ALL(IEEE_IS_FINITE(RAW))) RETURN
C     curl(B)[G/cm] * (1e-4 T/G)/(1e-2 m/cm) / mu0.
      FACTOR=0.01/(16.*ATAN(1.)*1.E-7)
      N=0
      DO K=1,NKP1
      DO J=1,NJP1
      DO F=1,3
        N=N+1
        OUT(N)=FACTOR*RAW(F,J,K)
      END DO
      END DO
      END DO
      STATUS=0
      END

      BLOCK DATA LFM_MI_CURL_INIT
      INTEGER VALID
      COMMON /LFM_MI_CURL_VALID/ VALID
      DATA VALID /0/
      END
