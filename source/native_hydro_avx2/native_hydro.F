      SUBROUTINE HYDRO(DT,BNORM,VOL1,FACEI,
     $       TRANS,
     $       BNqI,bsqqi,BXqI,BYqI,BZqI,
     $       bxqni,byqni,bzqni,
     $       RHO1,VX1,VY1,VZ1,C1,BX1,BY1,BZ1,
     $       RHO2,VX2,VY2,VZ2,C2,BX2,BY2,BZ2,
     $       RHO3,VX3,VY3,VZ3,C3,VX3B,VY3B,VZ3B)
#include "param.inc"
#include "mdims.inc"
#include "help.inc"
#include "global_dims.inc"
C
      DIMENSION BNORM(LIP1,LJ),VOL1(llow:LIhigh,llow:ljhigh),
     $  FACEI(LIP1,LJ),TRANS(LIP1,LJ,9),bsqqi(lip1,lj),
     $  BNqI(LIP1,LJ),BXqI(LIP1,LJ),BYqI(LIP1,LJ),BZqI(LIP1,LJ),
     $  BXqni(LIP1,LJ),BYqni(LIP1,LJ),BZqni(LIP1,LJ),
     $  RHO1(LLOW:LIHIGH,LLOW:LJHIGH),RHO2(LLOW:LIHIGH,LLOW:LJHIGH),
     $     RHO3(LI,LJ),
     $    VX1(LLOW:LIHIGH,LLOW:LJHIGH),VX2(LLOW:LIHIGH,LLOW:LJHIGH),
     $     VX3(LI,LJ),
     $   VY1(LLOW:LIHIGH,LLOW:LJHIGH),VY2(LLOW:LIHIGH,LLOW:LJHIGH),
     $     VY3(LI,LJ),
     $  VZ1(LLOW:LIHIGH,LLOW:LJHIGH),VZ2(LLOW:LIHIGH,LLOW:LJHIGH),
     $     VZ3(LI,LJ),
     $  C1(LLOW:LIHIGH,LLOW:LJHIGH),C2(LLOW:LIHIGH,LLOW:LJHIGH),
     $     C3(LI,LJ),
     $  BX1(LLOW:LIHIGH,LLOW:LJHIGH),BX2(LLOW:LIHIGH,LLOW:LJHIGH),
     $     BY1(LLOW:LIHIGH,LLOW:LJHIGH),BY2(LLOW:LIHIGH,LLOW:LJHIGH),
     $  BZ1(LLOW:LIHIGH,LLOW:LJHIGH),BZ2(LLOW:LIHIGH,LLOW:LJHIGH),
     $     VX3B(LI,LJ),VY3B(LI,LJ),VZ3B(LI,LJ)
C
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
C
C
      DIMENSION RH2(LLOW:LIHIGH,LJ),PX2(LLOW:LIHIGH,LJ),
     $             RHOINV(LIP1,LJ),VL2(LLOW:LIHIGH,LJ),
     $          E2(LLOW:LIHIGH,LJ),E(LIP1,LJ),P(LIP1,LJ)
     $            ,PY2(LLOW:LIHIGH,LJ),FPX2(LIP1,LJ),FPY2(LIP1,LJ),
     $          pz2(llow:lihigh,lj)
      DIMENSION VPLDUM(LIP1,LJ),VMNDUM(LIP1,LJ),SGNPL(LIP1,LJ),
     $          SGNMN(LIP1,LJ)
*      EQUIVALENCE (FMNMN,VPLDUM),(FMNPL,VMNDUM),(FPLMN,SGNPL),
*     $             (FPLPL,SGNMN)
C
*      EQUIVALENCE (RHOPL,RH2),(RHOMN,PX2),!   (RHOINV,VXPL),
*     $          (E2,VXMN),(E,CPL),(CMN,P)
*      EQUIVALENCE (VYPL,PY2)
      DATA PXSUM,PYSUM,PZSUM,BPXSUM,BPYSUM,BPZSUM/0.,0.,0.,0.,0.,0./
#ifdef DEBUG_MODE_ON
      write(*,*) "DEBUG:  in hydro.F::HYDRO(...)"
#endif
C
C     LOAD BASE TIME LEVEL ARRAYS
C
      DO 100 J = 1,MJ
      DO 100 I = -1,MI+2
         VOLQ0(I+2,J) = VOL1(I,J)
         RHO0(I+2,J) = RHO1(I,J)
         VX0(I+2,J) = VX1(I,J)
         VY0(I+2,J) = VY1(I,J)
         VZ0(I+2,J) = VZ1(I,J)
 100     C0(I+2,J) = C1(I,J)
C
*        CALL BR2VOL
*      CALL BR2RH
*      CALL BR2V(TRANS)
*      CALL BR2C
C
C
C     FIND INTERPOLATED VALUES AT INTERFACES
C
C
        DO 102 J=1,MJ
        DO 102 I=-ORD2+1,MI+ORD2
  102   VL2(I,J) = VOL1(I,J)
        CALL HIGH(VL2,VOLq,1)
        CALL CLIP(VOLQ,VOLQ0)
C
      DO 105 J = 1,MJ
      DO 105 I = -ORD2+1,MI+ORD2
  105 RH2(I,J) = RHO2(I,J)*vol1(i,j)
      CALL HIGH(RH2,RHO,2)
        do 107 j=1,mj
        do 107 i=1,mip1
  107   rho(i,j) = rho(i,j)/volq(i,j)
      CALL CLIP(RHO,RHO0)
      DO 110 J = 1,MJ
      DO 110 I = 1,MIP1
  110 RHOINV(I,J) = 1./(RHO(I,J)*volq(i,j))
C
C
      DO 120 J = 1,MJ
      DO 120 I = -ORD2+1,MI+ORD2
  120 PX2(I,J) = VX2(I,J)*rh2(i,j)
      CALL HIGH(PX2,VX,3)
        do 121 j=1,mj
        do 121 i=1,mip1
  121   vx(i,j) = vx(i,j)*rhoinv(I,j)
      CALL CLIP(VX,VX0)
C
      DO 133 J = 1,MJ
      DO 133 I = -ORD2+1,MI+ORD2
  133 PY2(I,J) = VY2(I,J)*rh2(I,j)
      CALL HIGH(PY2,VY,4)
        do 134 j=1,mj
        do 134 i=1,mip1
  134   vy(i,j) = vy(i,j)*rhoinv(I,j)
      CALL CLIP(VY,VY0)
C
      DO 138 J = 1,MJ
      DO 138 I = -ORD2+1,MI+ORD2
  138 PZ2(I,J) = VZ2(I,J)*rh2(i,j)
      CALL HIGH(PZ2,VZ,5)
        do 139 j=1,mj
        do 139 i=1,mip1
  139   vz(i,j) = vz(i,j)*rhoinv(I,j)
      CALL CLIP(VZ,VZ0)

      DO 140 J = 1,MJ
      DO 140 I = -ORD2+1,MI+ORD2
  140 E2(I,J) = RHo2(I,J)*(0.5*(VX2(I,J)**2+VY2(I,J)**2+
     $  VZ2(I,J)**2) + GGM1IN*C2(I,J)**2)
      CALL HIGH(E2,E,6)
        DO 141 J=1,MJ
        DO 141 I=1,MIP1
  141   E(I,J) = E(I,J)
      DO 150 J = 1,MJ
      DO 150 I = 1,MIP1
  150 P(I,J) = GGM1*(E(I,J)-0.5*RHO(I,J)*(VX(I,J)**2+VY(I,J)**2
     $          +vz(i,j)**2))
        DO 151 J=1,MJ
        DO 151 I=1,MIP4
  151   P0(I,J) = RHO0(I,J)*C0(I,J)**2
      CALL CLIP(P,P0)
C
      DO 170 J = 1,MJ
      DO 170 I = 1,MIP1
  170 C(I,J) = SQRT(P(I,J)*RHOINV(I,J))
      CALL CLIP(C,C0)
C
C
C    FIND VARIABLES ON BOTH SIDES OF THE INTERFACE:
C         IF FLOW DISCONTINUOUS, VARIABLES WILL DIFFER ON EITHER
C           SIDE.
C

      CALL PDM(RHO,RHO0,RHOPL,RHOMN)
      CALL PDM(VX,VX0,CPL,CMN)
      CALL PDM(VY,VY0,BXPL,BXMN)
      CALL PDM(VZ,VZ0,BYPL,BYMN)
*       if ( isweep .eq. 1 ) then
*         do 3200 i=1,mi
*         do 3200 j=1,mj
*           if ( (j .le. 3 .or. j .ge. mj-2 ) .and.
*     $           abs( vx(i,j) ) .lt. 7.e6 ) then
*             cpl(i,j) = vx(i,j)
*             cmn(i,j) = vx(i,j)
*             bxpl(i,j) = vy(i,j)
*             bxmn(i,j) = vy(i,j)
*             bypl(i,j) = vz(i,j)
*             bymn(i,j) = vz(i,j)
*           endif
* 3200    continue
*       endif
C
C   TRANSFORM TO INTERFACE NORMAL
C
        DO 195 J=1,MJ
        DO 195 I=1,MIP1
        VXPL(I,J) = TRANS(I,J,1)*CPL(I,J) + TRANS(I,J,2)*BXPL(I,J) +
     $                  TRANS(I,J,3)*BYPL(I,J)
        VYPL(I,J) = TRANS(I,J,4)*CPL(I,J) + TRANS(I,J,5)*BXPL(I,J) +
     $                  TRANS(I,J,6)*BYPL(I,J)
        VZPL(I,J) = TRANS(I,J,7)*CPL(I,J) + TRANS(I,J,8)*BXPL(I,J) +
     $                  TRANS(I,J,9)*BYPL(I,J)
  195   CONTINUE
        DO 196 J=1,MJ
        DO 196 I=1,MIP1
        VXMN(I,J) = TRANS(I,J,1)*CMN(I,J) + TRANS(I,J,2)*BXMN(I,J) +
     $                  TRANS(I,J,3)*BYMN(I,J)
        VYMN(I,J) = TRANS(I,J,4)*CMN(I,J) + TRANS(I,J,5)*BXMN(I,J) +
     $                  TRANS(I,J,6)*BYMN(I,J)
        VZMN(I,J) = TRANS(I,J,7)*CMN(I,J) + TRANS(I,J,8)*BXMN(I,J) +
     $                  TRANS(I,J,9)*BYMN(I,J)
  196 CONTINUE
C
C   make sure we don't have flow across the surface of the sphere
C
*        if (isweep .eq. 1) then
*          do 197 j=1,mj
*              vxpl(1,j) = 0.
*              vxmn(1,j) = 0.
*  197     continue
*        endif
C
      CALL PDM(C,C0,CPL,CMN)
*       if ( isweep .eq. 1 ) then
*         do 3201 i=1,mi
*         do 3201 j=1,mj
*           if ( (j .le. 3 .or. j .ge. mj-2 ) .and.
*     $           abs( vx(i,j) ) .lt. 7.e6 ) then
*              if ( vxpl(i,j) .gt. 0. ) then
*                 rhopl(i,j) = rhomn(i,j)
*                 cpl(i,j) = c(i,j)*sqrt(rho(i,j)/rhomn(i,j))
*                 cmn(i,j) = cpl(I,j)
*              else
*                 rhomn(i,j) = rhopl(i,j)
*                 cmn(i,j) = c(i,j)*sqrt(rho(i,j)/rhopl(i,j))
*                 cpl(i,j) = cmn(i,j)
*              endif
*           endif
* 3201    continue
*       endif
C
C
C   DEFINE VELOCITY BOUNDS OF DISTRIBUTIONS
C
C
      DO 200 J = 1,MJ
      DO 200 I = 1,MIP1
      AMN(I,J) = ALPHA*CMN(I,J)
      APL(I,J) = ALPHA*CPL(I,J)
      VPLDUM(I,J) = VXPL(I,J) - APL(I,J)
      VMNDUM(I,J) = VXMN(I,J) + AMN(I,J)
      VPLPL(I,J)= AMIN1(0.,VXPL(I,J)+APL(I,J))
      VPLMN(I,J) = AMIN1(0.,VPLDUM(I,J))
      VMNPL(I,J) = AMAX1(0.,VMNDUM(I,J))
  200 VMNMN(I,J)= AMAX1(0.,VXMN(I,J)-AMN(I,J))
C
C
C
      DO 250 J = 1,MJ
      DO 250 I = 1,MIP1
      APL(I,J) = 1./APL(I,J)
      AMN(I,J) = 1./AMN(I,J)
 250  continue
C
C
C     DETERMINE FLUXES
C
C
      DO 290 J=1,MJ
      DO 290 I=1,MIP1
      FMNPL(I,J)=AMN(I,J)*RHOMN(I,J)*VMNPL(I,J)
      FMNMN(I,J)=AMN(I,J)*RHOMN(I,J)*VMNMN(I,J)
      FPLPL(I,J)=APL(I,J)*RHOPL(I,J)*VPLPL(I,J)
      FPLMN(I,J)=APL(I,J)*RHOPL(I,J)*VPLMN(I,J)
  290 CONTINUE
C
C
      DO 300 J = 1,MJ
      DO 300 I = 1,MIP1
      FMNPL(I,J) = FMNPL(I,J)*VMNPL(I,J)
      FMNMN(I,J) = FMNMN(I,J)*VMNMN(I,J)
      FPLPL(I,J) = FPLPL(I,J)*VPLPL(I,J)
      FPLMN(I,J) = FPLMN(I,J)*VPLMN(I,J)
      FRHOMN(I,J) = 0.25*(FMNPL(I,J)-FMNMN(I,J))
      FRHOPL(I,J) = -0.25*(FPLMN(I,J)-FPLPL(I,J))
      FRHO1(I,J) = FRHOPL(I,J) + FRHOMN(I,J)
      FPY1(I,J) = FRHOPL(I,J)*VYPL(I,J)+FRHOMN(I,J)*VYMN(I,J)
      FPZ1(I,J) = FRHOPL(I,J)*VZPL(I,J)+FRHOMN(I,J)*VZMN(I,J)
  300 FE1(I,J)=FRHOPL(I,J)*(BETA*CPL(I,J)**2+
     $          0.5*(VYPL(I,J)**2+VZPL(I,J)**2))
     $   +FRHOMN(I,J)*(BETA*CMN(I,J)**2+0.5*(VYMN(I,J)**2+VZMN(I,J)**2))
C
C
C      MOMENTUM FLUX
C
C
      DO 400 J = 1,MJ
      DO 400 I = 1,MIP1
      FMNPL(I,J) = FMNPL(I,J)*VMNPL(I,J)
      FMNMN(I,J) = FMNMN(I,J)*VMNMN(I,J)
      FPLPL(I,J) = FPLPL(I,J)*VPLPL(I,J)
      FPLMN(I,J) = FPLMN(I,J)*VPLMN(I,J)
      FPX1(I,J) = 0.1666666*((FMNPL(I,J)-FMNMN(I,J))-(FPLMN(I,J)-
     $          FPLPL(I,J)))
  400 CONTINUE
C
C
C       ENERGY FLUX
C
C
      DO 500 J = 1,MJ
      DO 500 I = 1,MIP1
  500 FE1(I,J) = FE1(I,J) + 0.0625*(
     $       ( FMNPL(I,J)*VMNPL(I,J)-FMNMN(I,J)*VMNMN(I,J))
     $    -  ( FPLMN(I,J)*VPLMN(I,J)- FPLPL(I,J)*VPLPL(I,J) ) )
C
C
C   NOW TRANSFORM BACK TO THE BASE COORDINATE SYSTEM
C
        DO 600 J=1,MJ
        DO 600 I=1,MIP1
        FPZ2(I,J) = FPZ1(I,J)
        FPX2(I,J) = FPX1(I,J)
  600   FPY2(I,J) = FPY1(I,J)
C
        DO 650 J=1,MJ
        DO 650 I=1,MIP1
        FPX1(I,J) = TRANS(I,J,1)*FPX2(I,J) + TRANS(I,J,4)*FPY2(I,J)
     $          + TRANS(I,J,7)*FPZ2(I,J)
        FPY1(I,J) = TRANS(I,J,2)*FPX2(I,J) + TRANS(I,J,5)*FPY2(I,J)
     $          + TRANS(I,J,8)*FPZ2(I,J)
        FPZ1(I,J) = TRANS(I,J,3)*FPX2(I,J) + TRANS(I,J,6)*FPY2(I,J)
     $          + TRANS(I,J,9)*FPZ2(I,J)
  650   CONTINUE
C
C
C
C       NOW LET'S FIND THE MAGNETIC FORCES -- WITH THE ALFVEN CORRECTION
C               INCLUDED
C
#ifdef HYDROBUG
        if ( isweep .eq. 1 .and. level3 .eq. 5 ) then
           iunit = 15 + mype
           write (iunit,*) '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>',lstep
           lstep = lstep + 1 
           write (iunit,*) '<<<<<<<   rho  >>>>>>>',mip1
           write (iunit,953) (rho2(i,5),i=mip1-5,mip1+1)
           write (iunit,953) (rho1(i,5),i=mip1-5,mip1+1)
           write (iunit,953) (rho0(i+2,5),i=mip1-5,mip1+1)
           write (iunit,954) (rho(i,5),i=mip1-5,mip1)
           write (iunit,954) (rhopl(i,5),i=mip1-5,mip1)
           write (iunit,954) (rhomn(i,5),i=mip1-5,mip1)
           write (iunit,954) (frhopl(i,5),i=mip1-5,mip1)
           write (iunit,954) (frhomn(i,5),i=mip1-5,mip1)
           write (iunit,954) (frho1(i,5),i=mip1-5,mip1)
 953       format(8x,7(2x,1pe14.6))
 954       format(7(2x,1pe14.6))
           write (iunit,*) '<<<<<<<   c  >>>>>>>',mip1
           write (iunit,953) (c2(i,5),i=mip1-5,mip1+1)
           write (iunit,953) (c1(i,5),i=mip1-5,mip1+1)
           write (iunit,953) (c0(i+2,5),i=mip1-5,mip1+1)
           write (iunit,954) (c(i,5),i=mip1-5,mip1)
           write (iunit,954) (cpl(i,5),i=mip1-5,mip1)
           write (iunit,954) (cmn(i,5),i=mip1-5,mip1)
           write (iunit,954) (fe1(i,5),i=mip1-5,mip1)
c
           write (iunit,*) '<<<<<<<   vx  >>>>>>>',mip1
           write (iunit,953) (vx2(i,5),i=mip1-5,mip1+1)
           write (iunit,953) (vx1(i,5),i=mip1-5,mip1+1)
           write (iunit,953) (vx0(i+2,5),i=mip1-5,mip1+1)
           write (iunit,954) (vx(i,5),i=mip1-5,mip1)
           write (iunit,954) (vxpl(i,5),i=mip1-5,mip1)
           write (iunit,954) (vxmn(i,5),i=mip1-5,mip1)
           write (iunit,954) (fpx1(i,5),i=mip1-5,mip1)
        endif
#endif
        DO 700 J=1,MJ
        DO 700 I=-ord2+1,MI+ord2
        BZD(I,J) = BZ2(I,J)*VOL1(I,J)
        BXD(I,J) = BX2(I,J)*VOL1(I,J)
  700   BYD(I,J) = BY2(I,J)*VOL1(I,J)
C
        CALL HIGH(BXD,BX,7)
        CALL HIGH(BYD,BY,8)
        CALL HIGH(BZD,BZ,9)
C
        do 705 j=1,mj
        do 705 i=1,mip1
        rhoinv(i,j) = 1./volq(i,j)
        bx(i,j) =rhoinv(i,j)*bx(i,j)
        by(i,j) =rhoinv(i,j)*by(i,j)
        bz(i,j) =rhoinv(i,j)*bz(i,j)
  705   continue
C
C
        DO 710 J=1,MJ
        DO 710 I=-1,MI+2
        BZ0(I+2,J) = BZ1(I,J)
        BX0(I+2,J) = BX1(I,J)
  710   BY0(I+2,J) = BY1(I,J)
*        CALL BR2BX
*        CALL BR2BY
*        CALL BR2BZ
C
        CALL CLIP(BX,BX0)
        CALL CLIP(BY,BY0)
        CALL CLIP(BZ,BZ0)
        CALL PDM(BX,BX0,BXPL,BXMN)
        CALL PDM(BY,BY0,BYPL,BYMN)
        CALL PDM(BZ,BZ0,BZpl,BZmn)
*       if ( isweep .eq. 1 ) then
*         do 4200 i=1,mi
*         do 4200 j=1,mj
*           if ( (j .le. 3 .or. j .ge. mj-2 ) .and.
*     $           abs( vx(i,j) ) .lt. 7.e6 ) then
*             bzpl(i,j) = bz(i,j)
*             bzmn(i,j) = bz(i,j)
*             bxpl(i,j) = bx(i,j)
*             bxmn(i,j) = bx(i,j)
*             bypl(i,j) = by(i,j)
*             bymn(i,j) = by(i,j)
*           endif
* 4200    continue
*       endif
C
        DO 750 J=1,MJ
        DO 750 I=1,MIP1
        BSQPL(I,J) = PI4INV*(BXPL(I,J)**2+BYPL(I,J)**2+BZPL(I,J)**2+2.*
     $   (BXqI(I,J)*BXPL(I,J)+BYqI(I,J)*BYPL(I,J)+BZqI(I,J)*BZPL(I,J)))
     $           + 1.E-18
        BSQMN(I,J) = PI4INV*(BXMN(I,J)**2+BYMN(I,J)**2+BZMN(I,J)**2+2.*
     $   (BXqI(I,J)*BXMN(I,J)+BYqI(I,J)*BYMN(I,J)+BZqI(I,J)*BZMN(I,J)))
     $           + 1.E-18
        APL(I,J) = SQRT(1.3333333*((BSQPL(I,J)+BSQqI(I,J))/(
     $     RHOPL(I,J)+(bsqpl(i,j)+bsqqi(i,j))*cainv**2 )
     $          +CPL(I,J)**2))
        AMN(I,J) = SQRT(1.3333333*((BSQMN(I,J)+BSQqI(I,J))/(
     $     RHOMN(I,J)+(bsqmn(i,j)+bsqqi(i,j))*cainv**2 )
     $          +CMN(I,J)**2))
      VPLDUM(I,J) = VXPL(I,J) - APL(I,J)
      VMNDUM(I,J) = VXMN(I,J) + AMN(I,J)
      VPLPL(I,J)= AMIN1(0.,VXPL(I,J)+APL(I,J))-VXPL(I,J)
      VPLMN(I,J) = AMIN1(0.,VPLDUM(I,J))-VXPL(I,J)
      VMNPL(I,J) = AMAX1(0.,VMNDUM(I,J))-VXMN(I,J)
  750 VMNMN(I,J)= AMAX1(0.,VXMN(I,J)-AMN(I,J))-VXMN(I,J)
C
        DO 800 J=1,MJ
        DO 800 I=1,MIP1
        APL(I,J) = 0.5/APL(I,J)**3
        AMN(I,J) = 0.5/AMN(I,J)**3
        FRHOPL(I,J) = APL(I,J)*
     $  (VPLPL(I,J)**3-VPLMN(I,J)**3)
  800   FRHOMN(I,J) = AMN(I,J)*
     $  (VMNPL(I,J)**3-VMNMN(I,J)**3)
C
C
        DO 900 J=1,MJ
        DO 900 I=1,MIP1
        FBPX1(I,J) = FRHOPL(I,J)*(TRANS(I,J,1)*0.5*BSQPL(I,J)-
     $          PI4INV*(bnorm(I,j)*(BXPL(I,J)+BXqi(i,j))
     $            +bnqi(i,j)*BXPL(I,J) )) +
     $          FRHOMN(I,J)*(TRANS(I,J,1)*0.5*BSQMN(I,J)
     $          - PI4INV*(bnorm(I,j)*(BXMN(I,J)+BXqi(i,j))
     $            +bnqi(i,j)*bxmn(i,j) ) )
        FBPY1(I,J) = FRHOPL(I,J)*(TRANS(I,J,2)*0.5*BSQPL(I,J)-
     $          PI4INV*(bnorm(I,j)*(BYPL(I,J)+BYqi(i,j))
     $            +bnqi(i,j)*BYPL(I,J) )) +
     $           FRHOMN(I,J)*(TRANS(I,J,2)*0.5*BSQMN(I,J)
     $                  - PI4INV*(bnorm(I,j)*(BYMN(I,J)+BYqi(i,j))
     $            +bnqi(i,j)*BYMN(I,J) ))
  900   FBPZ1(I,J) = FRHOPL(I,J)*(TRANS(I,J,3)*0.5*BSQPL(I,J)-
     $                  PI4INV*(bnorm(I,j)*(BZPL(I,J)+BZqi(i,j))
     $            +bnqi(i,j)*BZPL(I,J) )) +
     $           FRHOMN(I,J)*(TRANS(I,J,3)*0.5*BSQMN(I,J)
     $                  - PI4INV*(bnorm(I,j)*(BZMN(I,J)+BZqi(i,j))
     $            +bnqi(i,j)*BZMN(I,J) ) )
C
        DO 950 J=1,MJ
        DO 950 I=1,MIP1
        FBPX1(I,J)=FBPX1(I,J)+PI4INV*(0.5*TRANS(I,J,1)*BSQqI(I,J)
     $          -BNqI(I,J)*BXqnI(I,J))
        FBPY1(I,J)=FBPY1(I,J)+PI4INV*(0.5*TRANS(I,J,2)*BSQqI(I,J)
     $          -BNqI(I,J)*BYqnI(I,J))
        FBPZ1(I,J)=FBPZ1(I,J)+PI4INV*(0.5*TRANS(I,J,3)*BSQqI(I,J)
     $          -BNqI(I,J)*BZqnI(I,J))
  950   CONTINUE
C
        DO 960 J=1,MJ
        DO 960 I=1,MIP1
        VX(I,J) = TRANS(I,J,1)*(VXmn(I,J)-VXpl(I,J)) +
     $            TRANS(I,J,4)*(VYmn(I,J)-VYpl(I,J))
     $          + TRANS(I,J,7)*(VZmn(I,J)-VZpl(I,J))
        VY(I,J) = TRANS(I,J,2)*(VXmn(I,J)-VXpl(I,J)) +
     $            TRANS(I,J,5)*(VYmn(I,J)-VYpl(I,J))
     $          + TRANS(I,J,8)*(VZmn(I,J)-VZpl(I,J))
        VZ(I,J) = TRANS(I,J,3)*(VXmn(I,J)-VXpl(I,J)) +
     $            TRANS(I,J,6)*(VYmn(I,J)-VYpl(I,J))
     $          + TRANS(I,J,9)*(VZmn(I,J)-VZpl(I,J))
        bsqpl(i,j)=0.5*(alcon*bsqqi(i,j)+0.5*(bsqpl(i,j)+bsqmn(i,j))*
     $          cainv**2)*ca
        fbpx1(i,j) = fbpx1(i,j) + bsqpl(i,j)*vx(i,j)
        fbpy1(i,j) = fbpy1(i,j) + bsqpl(i,j)*vy(i,j)
        fbpz1(i,j) = fbpz1(i,j) + bsqpl(i,j)*vz(i,j)
  960   CONTINUE

C
      DO 1000 J=1,MJ
      DO 1000 I=1,MI
        RHO3(I,J) = RHO3(I,J) + DT*(
     $          FACEI(I,J)*FRHO1(I,J)-FACEI(I+1,J)*FRHO1(I+1,J))
        VX3(I,J) = VX3(I,J) + DT*(
     $          FACEI(I,J)*FPX1(I,J)-FACEI(I+1,J)*FPX1(I+1,J))
        VY3(I,J) = VY3(I,J) + DT*(
     $          FACEI(I,J)*FPY1(I,J)-FACEI(I+1,J)*FPY1(I+1,J))
        VZ3(I,J) = VZ3(I,J) + DT*(
     $          FACEI(I,J)*FPZ1(I,J)-FACEI(I+1,J)*FPZ1(I+1,J))
        C3(I,J) = C3(I,J) + DT*(
     $          FACEI(I,J)*FE1(I,J)-FACEI(I+1,J)*FE1(I+1,J))
 1000   CONTINUE
C
*       DO 1017 J=1,MJ
*       PXSUM = PXSUM + FACEI(1,J)*FPX1(1,J)
*       PYSUM = PYSUM + FACEI(1,J)*FPY1(1,J)
*       PZSUM = PZSUM + FACEI(1,J)*FPZ1(1,J)
* 1017  CONTINUE
*       WRITE (9,*) 'PXSUM,PYSUM,PZSUM',PXSUM,PYSUM,PZSUM
C
        DO 1100 J=1,MJ
        DO 1100 I=1,MI
        VX3B(I,J) = VX3B(I,J) + DT*(
     $          FACEI(I,J)*FBPX1(I,J)-FACEI(I+1,J)*FBPX1(I+1,J))
        VY3B(I,J) = VY3B(I,J) + DT*(
     $          FACEI(I,J)*FBPY1(I,J)-FACEI(I+1,J)*FBPY1(I+1,J))
        VZ3B(I,J) = VZ3B(I,J) + DT*(
     $          FACEI(I,J)*FBPZ1(I,J)-FACEI(I+1,J)*FBPZ1(I+1,J))
 1100   CONTINUE
C
*       DO 1117 J=1,MJ
*       BPXSUM = BPXSUM + FACEI(1,J)*FBPX1(1,J)
*       BPYSUM = BPYSUM + FACEI(1,J)*FBPY1(1,J)
*       BPZSUM = BPZSUM + FACEI(1,J)*FBPZ1(1,J)
* 1117  CONTINUE
*       WRITE (9,*) 'BPXSUM,BPYSUM,BPZSUM',BPXSUM,BPYSUM,BPZSUM
C

      RETURN
      END
