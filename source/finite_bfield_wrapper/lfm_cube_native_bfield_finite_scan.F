C Storage/topology adapter; BFIELD,HIGH2,PDM2,XBHI,YBHI remain byte exact.
C ISWEEP=0 excludes the donor polar/periodic topology branches. Cube MPI
C owns shared edges; physical radial edge conditions belong to the driver.
      SUBROUTINE LFM_CUBE_BFIELD(NI_IN,NJ_IN,NG,BASE,PRED,XYZ,
     $ EDGE,FRAME,AREA,NORMAL,BFACE,BPRED,BG,EMF,STATUS)
      USE, INTRINSIC :: IEEE_ARITHMETIC, ONLY: IEEE_IS_FINITE
#include "param.inc"
#include "mdims.inc"
#include "global_dims.inc"
      INTEGER NI_IN,NJ_IN,NG,STATUS
      INTEGER LFM_CUBE_ALL_FINITE_FAST
      EXTERNAL LFM_CUBE_ALL_FINITE_FAST
      REAL BASE(1-NG:NI_IN+NG,1-NG:NJ_IN+NG,8)
      REAL PRED(1-NG:NI_IN+NG,1-NG:NJ_IN+NG,8)
      REAL XYZ(1-NG:NI_IN+NG,1-NG:NJ_IN+NG,3)
      REAL AREA(1-NG:NI_IN+NG,1-NG:NJ_IN+NG,2)
      REAL NORMAL(1-NG:NI_IN+NG,1-NG:NJ_IN+NG,3,2)
      REAL BFACE(1-NG:NI_IN+NG,1-NG:NJ_IN+NG,2)
      REAL BPRED(1-NG:NI_IN+NG,1-NG:NJ_IN+NG,2)
      REAL BG(NI_IN+1,NJ_IN+1,3),EDGE(NI_IN+1,NJ_IN+1)
      REAL FRAME(NI_IN+1,NJ_IN+1,9),EMF(NI_IN+1,NJ_IN+1)
      REAL QB(LLOW:LIHIGH,LLOW:LJHIGH,8)
      REAL QP(LLOW:LIHIGH,LLOW:LJHIGH,8)
      REAL XX(LLOW:LIHIGH,LLOW:LJHIGH,3)
      REAL AA(LLOW:LIHIGH,LLOW:LJHIGH,2)
      REAL NN(LLOW:LIHIGH,LLOW:LJHIGH,3,2)
      REAL BB(LLOW:LIHIGH,LLOW:LJHIGH,2)
      REAL PP(LLOW:LIHIGH,LLOW:LJHIGH,2)
      REAL GG(LIP1,LJP1,3),EE(LIP1,LJP1),TT(LIP1,LJP1,9)
      REAL LENGTH(LIP1,LJP1),ACC1(LIP1,LJ),ACC2(LI,LJP1)
      REAL BOUND(LIP1)
C Reuse process-local workspace, as the existing native scratch COMMONs.
C Active inputs/accumulators are still explicitly defined below; original
C BFIELD and sNaN/FPE/dynamic checks remain unchanged. Not thread-reentrant,
C just like the already shared native scratch/MDIMS/GDIMS path.
      SAVE QB,QP,XX,AA,NN,BB,PP,GG,EE,TT,LENGTH,ACC1,ACC2,BOUND
      STATUS=1
      EMF=0.
      IF(NG.NE.NO2.OR.NI_IN.LT.1.OR.NJ_IN.LT.1.OR.
     $ NI_IN.GT.LI.OR.NJ_IN.GT.LJ) RETURN
C Equivalent integer-AVX2 finite screening over ALL passed array words.
C Retain exact original IEEE invalid path, output-zeroing and traps.
      IF(LFM_CUBE_ALL_FINITE_FAST(BASE,SIZE(BASE,KIND=8)).EQ.0
     $ .OR.LFM_CUBE_ALL_FINITE_FAST(PRED,SIZE(PRED,KIND=8)).EQ.0
     $ .OR.LFM_CUBE_ALL_FINITE_FAST(BFACE,SIZE(BFACE,KIND=8)).EQ.0
     $ .OR.LFM_CUBE_ALL_FINITE_FAST(BPRED,SIZE(BPRED,KIND=8)).EQ.0)
     $ THEN
      IF(.NOT.ALL(IEEE_IS_FINITE(BASE)).OR.
     $ .NOT.ALL(IEEE_IS_FINITE(PRED)).OR.
     $ .NOT.ALL(IEEE_IS_FINITE(BFACE)).OR.
     $ .NOT.ALL(IEEE_IS_FINITE(BPRED))) RETURN
      END IF
      CALL MIJSET(NI_IN,NJ_IN)
      NL_I=NI_IN
      NL_J=NJ_IN
      ISWEEP=0
      LEVEL3=0
      MYPE=0
      ILOW_BOUND=.FALSE.
      IHIGH_BOUND=.FALSE.
      JLOW_BOUND=.FALSE.
      JHIGH_BOUND=.FALSE.
      KLOW_BOUND=.FALSE.
      KHIGH_BOUND=.FALSE.
C Storage-only: all locations read by the active MI/MJ stencil are
C overwritten below. BFIELD consumes only rho,Vx,Vy,Vz primitive channels.
C ISWEEP=0 disables branches that consume input EE/BOUND boundary values.
C Accumulators still require defined active values even for the DT=0 call.
      ACC1(1:NI_IN+1,1:NJ_IN)=0.
      ACC2(1:NI_IN,1:NJ_IN+1)=0.
      BOUND=0.
      QB(1-NG:NI_IN+NG,1-NG:NJ_IN+NG,1:4)=BASE(:,:,1:4)
      QP(1-NG:NI_IN+NG,1-NG:NJ_IN+NG,1:4)=PRED(:,:,1:4)
      XX(1-NG:NI_IN+NG,1-NG:NJ_IN+NG,:)=XYZ
      AA(1-NG:NI_IN+NG,1-NG:NJ_IN+NG,:)=AREA
      NN(1-NG:NI_IN+NG,1-NG:NJ_IN+NG,:,:)=NORMAL
      BB(1-NG:NI_IN+NG,1-NG:NJ_IN+NG,:)=BFACE
      PP(1-NG:NI_IN+NG,1-NG:NJ_IN+NG,:)=BPRED
      GG(1:NI_IN+1,1:NJ_IN+1,:)=BG
      TT(1:NI_IN+1,1:NJ_IN+1,:)=FRAME
      LENGTH(1:NI_IN+1,1:NJ_IN+1)=EDGE
      CALL BFIELD(0.,LENGTH,TT,AA(:,:,1),AA(:,:,2),
     $ XX(:,:,1),XX(:,:,2),XX(:,:,3),
     $ QB(:,:,1),QB(:,:,2),QB(:,:,3),QB(:,:,4),BB(:,:,1),BB(:,:,2),
     $ QP(:,:,1),QP(:,:,2),QP(:,:,3),QP(:,:,4),PP(:,:,1),PP(:,:,2),
     $ ACC1,ACC2,EE,BOUND,NN(:,:,1,2),NN(:,:,2,2),NN(:,:,3,2),
     $ NN(:,:,1,1),NN(:,:,2,1),NN(:,:,3,1),GG(:,:,1),GG(:,:,2),GG(:,:,3))
      EMF=EE(1:NI_IN+1,1:NJ_IN+1)
      IF (.NOT.ALL(IEEE_IS_FINITE(EMF))) RETURN
      STATUS=0
      END
