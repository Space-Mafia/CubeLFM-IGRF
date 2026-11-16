C Cube-only storage/flux adapter around unchanged original HYDRO/HIGH/PDM.
C Units and all physical coefficients are explicit caller inputs.
C Channels: mass, hydro momentum XYZ, energy, magnetic momentum XYZ.
C The last three are the original separate Boris-force accumulators, not
C magnetic-induction fluxes or an independently closed MHD update.
      SUBROUTINE LFM_CUBE_HYDRO_FLUX(NI_IN,NJ_IN,NG,DT_IN,
     $ GAMMA_IN,CA_IN,PDM_IN,BASE,PRED,VOL,AREA,FRAME,BN,BG,
     $ FLUX,ACC,STATUS)
      USE, INTRINSIC :: IEEE_ARITHMETIC, ONLY: IEEE_IS_FINITE
#include "param.inc"
#include "mdims.inc"
#include "help.inc"
      INTEGER NI_IN,NJ_IN,NG,STATUS
      INTEGER LFM_CUBE_ALL_FINITE_FAST
      EXTERNAL LFM_CUBE_ALL_FINITE_FAST
      INTEGER LFM_CUBE_FRAME_VALIDATION_HIT
      EXTERNAL LFM_CUBE_FRAME_VALIDATION_HIT
      EXTERNAL LFM_CUBE_FRAME_VALIDATION_RECORD
      REAL DT_IN,GAMMA_IN,CA_IN,PDM_IN
      REAL BASE(1-NG:NI_IN+NG,1-NG:NJ_IN+NG,8)
      REAL PRED(1-NG:NI_IN+NG,1-NG:NJ_IN+NG,8)
      REAL VOL(1-NG:NI_IN+NG,1-NG:NJ_IN+NG)
      REAL AREA(NI_IN+1,NJ_IN),FRAME(NI_IN+1,NJ_IN,9)
      REAL BN(NI_IN+1,NJ_IN),BG(NI_IN+1,NJ_IN,8)
      REAL FLUX(NI_IN+1,NJ_IN,8),ACC(NI_IN,NJ_IN,8)
      REAL QB(LLOW:LIHIGH,LLOW:LJHIGH,8)
      REAL QP(LLOW:LIHIGH,LLOW:LJHIGH,8)
      REAL VV(LLOW:LIHIGH,LLOW:LJHIGH)
      REAL AA(LIP1,LJ),TT(LIP1,LJ,9),BB(LIP1,LJ)
      REAL QQ(LIP1,LJ,8),UU(LI,LJ,8),CROSS(3),DOT
C Same persistent-workspace lifetime as the original scratch COMMONs.
C Every consumed active slab below is overwritten before use; keep the
C compiler sNaN/FPE checks, without refilling unused capacity every call.
      SAVE QB,QP,VV,AA,TT,BB,QQ,UU
      LOGICAL COEFS_READY
      SAVE COEFS_READY
      DATA COEFS_READY/.FALSE./
C Verbatim HYDRO scratch declaration follows: inspect source receipt.
      COMMON /SCRACH/
     $     RHO(LIP1,LJ),RHO0(LIP4,LJ),VX(LIP1,LJ),VX0(LIP4,LJ),
     $     VY0(LIP4,LJ),VZ0(LIP4,LJ),C(LIP1,LJ),C0(LIP4,LJ),
     $     VY(LIP1,LJ),P0(LIP4,LJ),BX0(LIP4,LJ),BY0(LIP4,LJ),
     $     BZ0(LIP4,LJ),VOLq(LIP1,LJ),VOLQ0(LIP4,LJ),
     $     VPLPL(LIP1,LJ),VPLMN(LIP1,LJ),VMNMN(LIP1,LJ),VMNPL(LIP1,LJ),
     $     FPLPL(LIP1,LJ),FPLMN(LIP1,LJ),FMNMN(LIP1,LJ),FMNPL(LIP1,LJ),
     $     RHOPL(LIP1,LJ),RHOMN(LIP1,LJ),VXPL(LIP1,LJ),VXMN(LIP1,LJ),
     $     VYPL(LIP1,LJ),VYMN(LIP1,LJ),CPL(LIP1,LJ),CMN(LIP1,LJ),
     $     BZPL(LIP1,LJ),BZMN(LIP1,LJ),APL(LIP1,LJ),AMN(LIP1,LJ),
     $     BX(LIP1,LJ),BY(LIP1,LJ),
     $     BXD(LLOW:LIHIGH,LJ),BYD(LLOW:LIHIGH,LJ),
     $     BSQPL(LIP1,LJ),BXPL(LIP1,LJ),BYPL(LIP1,LJ),BZ(LIP1,LJ),
     $     BSQMN(LIP1,LJ),BXMN(LIP1,LJ),BYMN(LIP1,LJ),
     $     BZD(LLOW:LIHIGH,LJ),VZ(LIP1,LJ),
     $     FRHOMN(LIP1,LJ),FRHOPL(LIP1,LJ),FRHO1(LIP1,LJ),FPX1(LIP1,LJ),
     $     FPY1(LIP1,LJ),FPZ1(LIP1,LJ),FE1(LIP1,LJ),FBPX1(LIP1,LJ),
     $     FBPY1(LIP1,LJ),FBPZ1(LIP1,LJ),VZMN(LIP1,LJ),VZPL(LIP1,LJ),
     $     FPZ2(LIP1,LJ)

      STATUS=1
      FLUX=0.
      ACC=0.
      IF (NG .NE. NO2 .OR. NI_IN.LT.1 .OR. NI_IN.GT.LI .OR.
     $ NJ_IN.LT.1 .OR. NJ_IN.GT.LJ) RETURN
C Integer-vector screening checks every passed word of the same arrays.
C Any nonfinite input follows the EXACT original IEEE rejection expression
C below, retaining its status/output and floating-trap behavior.
      IF (LFM_CUBE_ALL_FINITE_FAST(BASE,SIZE(BASE,KIND=8)).EQ.0
     $ .OR.LFM_CUBE_ALL_FINITE_FAST(PRED,SIZE(PRED,KIND=8)).EQ.0
     $ .OR.LFM_CUBE_ALL_FINITE_FAST(VOL,SIZE(VOL,KIND=8)).EQ.0
     $ .OR.LFM_CUBE_ALL_FINITE_FAST(AREA,SIZE(AREA,KIND=8)).EQ.0
     $ .OR.LFM_CUBE_ALL_FINITE_FAST(FRAME,SIZE(FRAME,KIND=8)).EQ.0
     $ .OR.LFM_CUBE_ALL_FINITE_FAST(BN,SIZE(BN,KIND=8)).EQ.0
     $ .OR.LFM_CUBE_ALL_FINITE_FAST(BG,SIZE(BG,KIND=8)).EQ.0) THEN
      IF (.NOT.ALL(IEEE_IS_FINITE(BASE)) .OR.
     $ .NOT.ALL(IEEE_IS_FINITE(PRED)) .OR.
     $ .NOT.ALL(IEEE_IS_FINITE(VOL)) .OR.
     $ .NOT.ALL(IEEE_IS_FINITE(AREA)) .OR.
     $ .NOT.ALL(IEEE_IS_FINITE(FRAME)) .OR.
     $ .NOT.ALL(IEEE_IS_FINITE(BN)) .OR.
     $ .NOT.ALL(IEEE_IS_FINITE(BG))) RETURN
      END IF
      IF (.NOT.ALL(IEEE_IS_FINITE([DT_IN,GAMMA_IN,CA_IN,PDM_IN])))
     $ RETURN
      IF (GAMMA_IN.LE.1. .OR. CA_IN.LE.0. .OR. PDM_IN.LT.0.
     $ .OR. DT_IN.LT.0.) RETURN
      IF (ANY(BASE(:,:,1).LE.0.) .OR. ANY(PRED(:,:,1).LE.0.) .OR.
     $ ANY(BASE(:,:,5).LE.0.) .OR.
     $ ANY(VOL.LE.0.) .OR. ANY(AREA.LE.0.)) RETURN
C Native MHD_HACK's extrapolated C1 is signed. Exact HYDRO squares it;
C its physical interface sound speed is clipped against positive BASE C.
      IF(LFM_CUBE_FRAME_VALIDATION_HIT(FRAME,
     $ SIZE(FRAME,KIND=8)).EQ.0) THEN
      DO J=1,NJ_IN
        DO I=1,NI_IN+1
          DO K=1,3
            DO L=1,3
              DOT=SUM(FRAME(I,J,3*K-2:3*K)*
     $                    FRAME(I,J,3*L-2:3*L))
              IF(K.EQ.L) DOT=DOT-1.
              IF(ABS(DOT).GT.2.E-12) RETURN
            END DO
          END DO
          CROSS(1)=FRAME(I,J,2)*FRAME(I,J,6)-
     $             FRAME(I,J,3)*FRAME(I,J,5)
          CROSS(2)=FRAME(I,J,3)*FRAME(I,J,4)-
     $             FRAME(I,J,1)*FRAME(I,J,6)
          CROSS(3)=FRAME(I,J,1)*FRAME(I,J,5)-
     $             FRAME(I,J,2)*FRAME(I,J,4)
          IF(MAXVAL(ABS(CROSS-FRAME(I,J,7:9))).GT.2.E-12) RETURN
        END DO
      END DO
      CALL LFM_CUBE_FRAME_VALIDATION_RECORD(FRAME,
     $ SIZE(FRAME,KIND=8))
      END IF
      CALL MIJSET(NI_IN,NJ_IN)
      ORDER=NO
      ORD2=ORDER/2
      GAMMA=GAMMA_IN
      GAMM1=GAMMA-1.
      GMM1IN=1./GAMM1
      GMQGM1=GAMMA/GAMM1
      ALPHA=SQRT(3./GAMMA)
      PI=ACOS(0.)*2.
      PI4=4.*PI
      PI8=8.*PI
      PI4INV=1./PI4
      PI8INV=1./PI8
      CA=CA_IN
      CAINV=1./CA
      ALCON=PI4INV*CAINV**2
      TWELVI=1./12.
      GMINV=1./GAMMA
      GM1=GAMMA-1.
      GP1=GAMMA+1.
      GM1INV=1./GM1
      WIDTH=SQRT(3./GAMMA)
      GGM1=GAMMA*GM1
      GGM1IN=1./GGM1
      BETA=(3.-GAMMA)*0.5*GGM1IN
      PDMB=PDM_IN
      IF (.NOT.COEFS_READY) THEN
        CALL DERCO(FCOEF,ORDER,1)
        COEFS_READY=.TRUE.
      END IF
C Storage-only: HYDRO/HIGH sweep J=1:MJ and never read J ghosts.
C Every input location that the unchanged kernels read is copied below.
C Do not zero unused64x64 capacity padding on each small active plane.
      UU(1:NI_IN,1:NJ_IN,:)=0.
      QB(1-NG:NI_IN+NG,1:NJ_IN,:)=BASE(:,1:NJ_IN,:)
      QP(1-NG:NI_IN+NG,1:NJ_IN,:)=PRED(:,1:NJ_IN,:)
      VV(1-NG:NI_IN+NG,1:NJ_IN)=VOL(:,1:NJ_IN)
      AA(1:NI_IN+1,1:NJ_IN)=AREA
      TT(1:NI_IN+1,1:NJ_IN,:)=FRAME
      BB(1:NI_IN+1,1:NJ_IN)=BN
      QQ(1:NI_IN+1,1:NJ_IN,:)=BG
      CALL HYDRO(DT_IN,BB,VV,AA,TT,
     $ QQ(:,:,1),QQ(:,:,2),QQ(:,:,3),QQ(:,:,4),QQ(:,:,5),
     $ QQ(:,:,6),QQ(:,:,7),QQ(:,:,8),
     $ QB(:,:,1),QB(:,:,2),QB(:,:,3),QB(:,:,4),QB(:,:,5),
     $ QB(:,:,6),QB(:,:,7),QB(:,:,8),
     $ QP(:,:,1),QP(:,:,2),QP(:,:,3),QP(:,:,4),QP(:,:,5),
     $ QP(:,:,6),QP(:,:,7),QP(:,:,8),
     $ UU(:,:,1),UU(:,:,2),UU(:,:,3),UU(:,:,4),UU(:,:,5),
     $ UU(:,:,6),UU(:,:,7),UU(:,:,8))
      DO J=1,NJ_IN
        DO I=1,NI_IN+1
          FLUX(I,J,1)=AREA(I,J)*FRHO1(I,J)
          FLUX(I,J,2)=AREA(I,J)*FPX1(I,J)
          FLUX(I,J,3)=AREA(I,J)*FPY1(I,J)
          FLUX(I,J,4)=AREA(I,J)*FPZ1(I,J)
          FLUX(I,J,5)=AREA(I,J)*FE1(I,J)
          FLUX(I,J,6)=AREA(I,J)*FBPX1(I,J)
          FLUX(I,J,7)=AREA(I,J)*FBPY1(I,J)
          FLUX(I,J,8)=AREA(I,J)*FBPZ1(I,J)
        END DO
      END DO
      ACC=UU(1:NI_IN,1:NJ_IN,:)
      IF(.NOT.ALL(IEEE_IS_FINITE(FLUX)) .OR.
     $ .NOT.ALL(IEEE_IS_FINITE(ACC))) RETURN
      STATUS=0
      END
