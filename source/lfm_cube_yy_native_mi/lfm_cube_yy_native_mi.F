      SUBROUTINE LFM_CUBE_YY_LOAD_PACKED(STEP,TIME_S,DT_S,
     $ STATE,COUNT,STATUS,CPCP_N_KV,CPCP_S_KV,MAX_FAC_UA_M2)
      USE, INTRINSIC :: IEEE_ARITHMETIC
      USE LFM_ROTATING_MI_STATE, ONLY: ROTATING_MI
#include "lfm_sw_state.inc"
#include "ion_bridge_common.inc"
      INTEGER STEP,COUNT,STATUS,I,J,K,N,FIELD,NFIELD,KSRC
      REAL TIME_S,DT_S,STATE(COUNT)
      REAL CPCP_N_KV,CPCP_S_KV,MAX_FAC_UA_M2

      STATUS=1
      CPCP_N_KV=0.0
      CPCP_S_KV=0.0
      MAX_FAC_UA_M2=0.0
      IF (.NOT.LFM_IONOSPHERE_ENABLED) RETURN
      IF (COUNT .NE. (2*ION_NI+3*ION_NIP1)*ION_NJ*ION_NK) THEN
        WRITE(*,'(A,2(1X,I0))') 'packed ionosphere state mismatch:',
     $    COUNT,(2*ION_NI+3*ION_NIP1)*ION_NJ*ION_NK
        RETURN
      END IF

      N=0
      DO FIELD=1,5
        IF (FIELD .LE. 2) THEN
          NFIELD=ION_NI
        ELSE
          NFIELD=ION_NIP1
        END IF
        DO K=1,ION_NK
          DO J=1,ION_NJ
            DO I=1,NFIELD
              N=N+1
              IF (FIELD .EQ. 1) ION_RHO(I,J,K)=STATE(N)
              IF (FIELD .EQ. 2) ION_CS(I,J,K)=STATE(N)
              IF (FIELD .EQ. 3) ION_BX(I,J,K)=STATE(N)
              IF (FIELD .EQ. 4) ION_BY(I,J,K)=STATE(N)
              IF (FIELD .EQ. 5) ION_BZ(I,J,K)=STATE(N)
            END DO
          END DO
        END DO
      END DO

      DO K=1,ION_NK
        KSRC=1+MODULO(K-1+ION_NK/2,ION_NK)
        DO I=1,ION_NIP1
          ION_BX(I,ION_NJP1,K)=ION_BX(I,ION_NJ,KSRC)
          ION_BY(I,ION_NJP1,K)=ION_BY(I,ION_NJ,KSRC)
          ION_BZ(I,ION_NJP1,K)=ION_BZ(I,ION_NJ,KSRC)
        END DO
      END DO
      DO J=1,ION_NJP1
        DO I=1,ION_NIP1
          ION_BX(I,J,ION_NKP1)=ION_BX(I,J,1)
          ION_BY(I,J,ION_NKP1)=ION_BY(I,J,1)
          ION_BZ(I,J,ION_NKP1)=ION_BZ(I,J,1)
        END DO
      END DO

      ION_TIME_BR=TIME_S
      ION_DT_BR=DT_S
      ION_LSTEP_BR=STEP
C     Current is supplied by the qualified released YY operator on Cube.
C     Keep original common-state loads/time; do not also run native JPARA.
      STATUS=0
      END

C Cap geometry is the existing native LFM mesh, not a replacement Poisson grid.
      SUBROUTINE LFM_CUBE_YY_CAP_GEOMETRY(X,Y,TH,RAD,PH,STATUS)
#include "ion-param.inc"
#include "ion_bridge_common.inc"
#include "iongeo.inc"
      REAL X(NJP1,NK2P3),Y(NJP1,NK2P3),TH(NJP1),RAD(NJP1)
      REAL PH(NKP1)
      INTEGER STATUS,K
      X=X2ION
      Y=Y2ION
C Released mi_current_shell=1 is zero-based: SECOND active cell centre.
C The former ION_*(2) vertex shell mixed first-cell wall curl into FAC.
C Midpoint of vertices 1/2 is the existing mapped radial cell-2 surface.
      TH=ATAN2(SQRT((ION_Y(2,:,1)+ION_Y(3,:,1))**2+
     $ (ION_Z(2,:,1)+ION_Z(3,:,1))**2),
     $ ION_X(2,:,1)+ION_X(3,:,1))
      RAD=.005*SQRT((ION_X(2,:,1)+ION_X(3,:,1))**2+
     $ (ION_Y(2,:,1)+ION_Y(3,:,1))**2+
     $ (ION_Z(2,:,1)+ION_Z(3,:,1))**2)
      DO K=1,NKP1
        PH(K)=ATAN2(ION_Z(1,2,K),ION_Y(1,2,K))
        IF(K.GT.1) THEN
          DO WHILE(PH(K).LE.PH(K-1))
            PH(K)=PH(K)+8.*ATAN(1.)
          END DO
        ENDIF
      END DO
      STATUS=0
      END

C Explicit per-run user contract after native geometry init, before any solve.
      SUBROUTINE LFM_CUBE_YY_MI100_POLICY(STATUS)
#include "ion_bridge_common.inc"
      INTEGER STATUS
      STATUS=1
      IF(ION_PED0_BR.NE.5..OR.ION_HALLMIN_BR.NE.0.) RETURN
      ION_ITMAX_BR=100
      ION_ERRMAX_BR=1.E-3
      STATUS=0
      END

      SUBROUTINE LFM_CUBE_YY_RECORD_SOLVE(H)
      USE, INTRINSIC :: IEEE_ARITHMETIC
      INTEGER H,LASTITER,LASTMAX,ITER(2),MAXIMUM(2),CONVERGED(2)
      REAL LASTSER,SER(2)
      COMMON /CUBE_YY_SOLVE_LAST_R/ LASTSER
      COMMON /CUBE_YY_SOLVE_LAST_I/ LASTITER,LASTMAX
      COMMON /CUBE_YY_SOLVE_METRICS_R/ SER
      COMMON /CUBE_YY_SOLVE_METRICS_I/ ITER,MAXIMUM,CONVERGED
      SER(H)=LASTSER
      ITER(H)=LASTITER
      MAXIMUM(H)=LASTMAX
      CONVERGED(H)=0
      IF(IEEE_IS_FINITE(LASTSER).AND.LASTSER.LE.1.E-3)
     $ CONVERGED(H)=1
      END

      SUBROUTINE LFM_CUBE_YY_SOLVE_METRICS(SER_OUT,ITER_OUT,
     $ MAX_OUT,CONV_OUT)
      INTEGER ITER_OUT(2),MAX_OUT(2),CONV_OUT(2)
      INTEGER ITER(2),MAXIMUM(2),CONVERGED(2)
      REAL SER_OUT(2),SER(2)
      COMMON /CUBE_YY_SOLVE_METRICS_R/ SER
      COMMON /CUBE_YY_SOLVE_METRICS_I/ ITER,MAXIMUM,CONVERGED
      SER_OUT=SER
      ITER_OUT=ITER
      MAX_OUT=MAXIMUM
      CONV_OUT=CONVERGED
      END
