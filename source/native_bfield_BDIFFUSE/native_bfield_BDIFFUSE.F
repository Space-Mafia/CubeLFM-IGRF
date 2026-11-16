#define BDIFFUSE 1
      SUBROUTINE bfieLD(DT,EDGEK,
     $          TRANS2,facei,facej,
     $          X,Y,Z,
     $          RHO,VX,VY,VZ,BI,BJ,
     $          RHO1,VX1,VY1,VZ1,BI1,BJ1,
     $                  BI2,BJ2,e2,bound_mix,
     $          xnormj,ynormj,znormj,
     $          xnormi,ynormi,znormi,
     $          bxc0,byc0,bzc0)
!
! ... Local variables ..................................................
!
#include "param.inc"
#include "global_dims.inc"
#include "mdims.inc"
#include "help.inc"
#include "meter.inc"
        DIMENSION 
     $       EDGEK(lip1,Ljp1),
     $       TRANS2(lip1,ljp1,9),
     $          RHO(LLOW:LIHIGH,LLOW:LJHIGH),
     $       VX(LLOW:LIHIGH,LLOW:LJHIGH),
     $       VY(LLOW:LIHIGH,LLOW:LJHIGH),VZ(LLOW:LIHIGH,LLOW:LJHIGH),
     $       BI(llow:lihigh,llow:LJhigh),BJ(LLOW:LIHIGH,LLOW:LJHIGH),
     $          RHO1(LLOW:LIHIGH,LLOW:LJHIGH),
     $       VX1(LLOW:LIHIGH,LLOW:LJHIGH),VY1(LLOW:LIHIGH,LLOW:LJHIGH),
     $       VZ1(LLOW:LIHIGH,LLOW:LJHIGH),
     $       BI1(llow:lihigh,Llow:ljhigh),BJ1(LLOW:LIHIGH,LLOW:LJHIGH),
     $          BI2(lip1,LJ),BJ2(li,ljp1),e2(lip1,ljp1),
     $       bound_mix(lip1),
     $          X(llow:lihigh,LLOW:LJHIGH),Y(llow:lihigh,LLOW:LJHIGH),
     $       Z(llow:lihigh,LLOW:LJHIGH),
     $          bxc0(lip1,ljp1),byc0(lip1,ljp1),bzc0(lip1,ljp1),
     $       facei(llow:lihigh,llow:ljhigh),
     $       facej(llow:lihigh,llow:ljhigh),
     $       xnormi(llow:lihigh,llow:ljhigh),
     $       ynormi(llow:lihigh,llow:ljhigh),
     $       znormi(llow:lihigh,llow:ljhigh),
     $       xnormj(llow:lihigh,llow:ljhigh),
     $       ynormj(llow:lihigh,llow:ljhigh),
     $       znormj(llow:lihigh,llow:ljhigh)
C
      COMMON /SCRACH/
     $     VX0(LIP4,LJP4),VY0(LIP4,LJP4),VZ0(LIP4,LJP4),
     $     VXBAR(LIP1,LJP1),VYBAR(LIP1,LJP1),VZBAR(LIP1,LJP1),
     $     VXQ00(LIP1,LJP1),VXQ01(LIP1,LJP1),VXQ10(LIP1,LJP1),
     $     VXQ11(LIP1,LJP1),VYQ00(LIP1,LJP1),VYQ01(LIP1,LJP1),
     $     VYQ10(LIP1,LJP1),VYQ11(LIP1,LJP1),VZQ00(LIP1,LJP1),
     $     VZQ01(LIP1,LJP1),VZQ10(LIP1,LJP1),VZQ11(LIP1,LJP1),
     $     vxCORN(LIp1,LJp1),VYCORN(LIP1,LJP1),VZCORN(LIP1,LJP1),
     $     DET(LIP1,LJP1),VBARQX(LIP1,LJP1),VBARQY(LIP1,LJP1),
     $     MASS(LI,LJ),PX(LI,LJ),PY(LI,LJ),PZ(LI,LJ),
     $     MCORN(LIP1,LJP1),
     $     DXQI(LIP2,LJP2),DXQJ(LIP2,LJP2),DYQI(LIP2,LJP2),
     $     DYQJ(LIP2,LJP2),SDENOM(LIP1,LJP1),
     $     DXQA(LIP2,LJP2),DXQB(LIP2,LJP2),DELDX(LIP1,LJP1),
     $     DELDY(LIP1,LJP1),DYQA(LIP2,LJP2),DYQB(LIP2,LJP2),
     $     E(LIP1,LJP1),BII(LIP1,LJP1),BJJ(LIP1,LJP1),
     $     BIQ(llow:lihigh,llow:ljhigh),BJQ(llow:lihigh,llow:ljhigh),
     $     VOLQJ(LIP1,LJP1),
     $     JZ(LIP1,LJP1),BXCORN(LIP1,LJP1),BYCORN(LIP1,LJP1),
     $     BZCORN(LIP1,LJP1)
c
      COMMON /SCRACH/
     $    DZQI(LIP2,LJP2),DZQJ(LIP2,LJP2),DXDI(LIP1,LJP1),
     $    DXDJ(LIP1,LJP1),DYDI(LIP1,LJP1),DYDJ(LIP1,LJP1),
     $    SQI(LIP1,LJP1),SQJ(LIP1,LJP1),
     $     DZQA(LIP1,LJP1),DZQB(LIP1,LJP1),
     $     SQ00(LIP1,LJP1),SQ01(LIP1,LJP1),
     $    SQ10(LIP1,LJP1),SQ11(LIP1,LJP1),
     $    DXA(LIP1,LJP1),DXB(LIP1,LJP1),DYA(LIP1,LJP1),DYB(LIP1,LJP1),
     $    FII(LIP1,LJP1),FJJ(LIP1,LJP1),xnQi(LIp1,LJp1),
     $    ynQi(LIp1,LJp1),xnQj(LIp1,LJp1),ynQj(LIp1,LJp1),
     $     znqi(lip1,ljp1),znqj(lip1,ljp1),
     $    bxq0(lip1,ljp1),byq0(lip1,ljp1)
c
      DIMENSION BI0(LIP1,LJP4),BJ0(LIP4,LJP1),BILFT(LIP1,LJP1),
     $     BIRT(LIP1,LJP1),BJLFT(LIP1,LJP1),BJRT(LIP1,LJP1)
      DIMENSION DVXQ00(LIP1,LJP1),DVYQ00(LIP1,LJP1),DVZQ00(LIP1,LJP1)
      DIMENSION DVZZ(LIP1,LJP1)
*      REAL JZ
      REAL MASS,MCORN
C
*        EQUIVALENCE (BI0,VX0),(BJ0,VY0)
*        EQUIVALENCE (BILFT,VXQ00),(BIRT,VXQ10),(BJLFT,VXQ01),
*     $          (BJRT,VXQ11)
!
! ... Begin ............................................................
!
#ifdef DEBUG_MODE_ON
        write(*,*) "DEBUG:  in bfield.F::bfieLD(...)"
#endif
        DO 250 J=-1,MJ+2
        DO 250 I=-1,MI+2
        VX0(I+2,J+2) = VX(I,J)
        VZ0(I+2,J+2) = VZ(I,J)
  250   VY0(I+2,J+2) = VY(I,J)
C
*        CALL BR2qV(TRANSI)
C
C
*        CALL HIGH2(MASS,MCORN,2)
        CALL HIGH2(VX1,VXCORN,3)
        CALL HIGH2(VY1,VYCORN,4)
        CALL HIGH2(VZ1,VZCORN,5)
C
*       DO 350 J=1,MJP1
*       DO 350 I=1,MIP1
*       VXCORN(I,J) = VXCORN(I,J)/MCORN(I,J)
*       VZCORN(I,J) = VZCORN(I,J)/MCORN(I,J)
*  350  VYCORN(I,J) = VYCORN(I,J)/MCORN(I,J)
C
        CALL PDM2(VXCORN,VX0,VXQ00,VXQ01,VXQ10,VXQ11)
        CALL PDM2(VYCORN,VY0,VYQ00,VYQ01,VYQ10,VYQ11)
        CALL PDM2(VZCORN,VZ0,VZQ00,VZQ01,VZQ10,VZQ11)
C
C
        DO 400 J=1,MJP1
        DO 400 I=1,MIP1
        VXBAR(I,J) = VX0(I+1,J+1)+VX0(I+1,J+2)+VX0(I+2,J+1)+VX0(I+2,J+2)
        VYBAR(I,J) = VY0(I+1,J+1)+VY0(I+1,J+2)+VY0(I+2,J+1)+VY0(I+2,J+2)
        VZBAR(I,J) = VZ0(I+1,J+1)+VZ0(I+1,J+2)+VZ0(I+2,J+1)+VZ0(I+2,J+2)
  400   CONTINUE
C
C
        DO 500 J=0,MJ+1
        DO 500 I=0,MI+1
        DXQI(I+1,J+1) = X(I+1,J)-X(I,J)
        DXQJ(I+1,J+1) = X(I,J+1)-X(i,J)
        DYQI(I+1,J+1) = Y(I+1,J)-Y(i,J)
        DYQJ(I+1,J+1) = Y(I,J+1)-Y(I,J)
        DZQI(I+1,J+1) = Z(I+1,J)-Z(I,J)
        DZQJ(I+1,J+1) = Z(I,J+1)-Z(I,J)
  500   CONTINUE
C
C
C
*        CALL BDZZ(DXQI,DXQJ,DYQI,DYQJ,DZQI,DZQJ)
C
C  TRANSFORM TO THE NORMAL COORDINATES
C
        DO 600 J=1,MJP1
        DO 600 I=1,MIP1
        DXDI(I,J) = 0.5*(TRANS2(I,J,1)*(DXQI(I,J+1)+DXQI(I+1,J+1))
     $          + TRANS2(I,J,2)*(DYQI(I,J+1)+DYQI(I+1,J+1))
     $               +TRANS2(I,J,3)*(DZQI(I,J+1)+DZQI(I+1,J+1)))
        DYDI(I,J) = 0.5*( TRANS2(I,J,4)*(DXQI(I,J+1)+DXQI(I+1,J+1))
     $          + TRANS2(I,J,5)*(DYQI(I,J+1)+DYQI(I+1,J+1))
     $               + TRANS2(I,J,6)*(DZQI(I,J+1)+DZQI(I+1,J+1)))
        DXDJ(I,J) = 0.5*( TRANS2(I,J,1)*(DXQJ(I+1,J)+DXQJ(I+1,J+1))
     $          + TRANS2(I,J,2)*(DYQJ(I+1,J)+DYQJ(I+1,J+1))
     $               + TRANS2(I,J,3)*(DZQJ(I+1,J)+DZQJ(I+1,J+1)))
        DYDJ(I,J) = 0.5*( TRANS2(I,J,4)*(DXQJ(I+1,J)+DXQJ(I+1,J+1))
     $          + TRANS2(I,J,5)*(DYQJ(I+1,J)+DYQJ(I+1,J+1))
     $               +TRANS2(I,J,6)*(DZQJ(I+1,J)+DZQJ(I+1,J+1)))
        VBARQX(I,J) = TRANS2(I,J,1)*VXBAR(I,J) +
     $          TRANS2(I,J,2)*VYBAR(I,J) + TRANS2(I,J,3)*VZBAR(I,J)
        VBARQY(I,J) = TRANS2(I,J,4)*VXBAR(I,J) +
     $          TRANS2(I,J,5)*VYBAR(I,J) + TRANS2(I,J,6)*VZBAR(I,J)
  600   CONTINUE
C
C
        DO 620 J=1,MJP1
        DO 620 I=1,MIP1
        DET(I,J) = sign(1.,DXDI(I,J)*DYDJ(I,J)-DXDJ(I,J)*DYDI(I,J))
        SQI(I,J) = SIGN(0.5,DET(I,J)*(
     $          VBARQY(I,J)*DXDJ(I,J)-VBARQX(I,J)*DYDJ(I,J)))
        SQJ(I,J) = SIGN(0.5,DET(I,J)*(
     $          VBARQX(I,J)*DYDI(I,J)-VBARQY(I,J)*DXDI(I,J)))
  620   CONTINUE

        DO 700 J=1,MJP1
        DO 700 I=1,MIP1
        SQ00(I,J) = AMAX1(0.,-SQI(I,J)-SQJ(I,J))
        SQ01(I,J) = AMAX1(0.,-SQI(I,J)+SQJ(I,J))
        SQ10(I,J) = AMAX1(0.,SQI(I,J)-SQJ(I,J))
        SQ11(I,J) = AMAX1(0.,SQI(I,J)+SQJ(I,J))
        SDENOM(I,J) = 1./(SQ00(I,J)+SQ01(I,J)+SQ10(I,J)+SQ11(I,J))
  700   CONTINUE
C
        do 750 j=1,mjp1
           do 750 i=1,mip1 
              dvvx00 = trans2(i,j,1)*vxq00(i,j)+trans2(i,j,2)*vyq00(i,j)
     $             + trans2(i,j,3)*vzq00(i,j)
              dvvx10 = trans2(i,j,1)*vxq10(i,j)+trans2(i,j,2)*vyq10(i,j)
     $             + trans2(i,j,3)*vzq10(i,j)
              dvvx01 = trans2(i,j,1)*vxq01(i,j)+trans2(i,j,2)*vyq01(i,j)
     $             + trans2(i,j,3)*vzq01(i,j)
              dvvx11 = trans2(i,j,1)*vxq11(i,j)+trans2(i,j,2)*vyq11(i,j)
     $             + trans2(i,j,3)*vzq11(i,j)
        dvxq00(i,j)=amax1(dvvx00,dvvx01,dvvx10,dvvx11)
     $             -amin1(dvvx00,dvvx01,dvvx10,dvvx11)
              dvvy00 = trans2(i,j,4)*vxq00(i,j)+trans2(i,j,5)*vyq00(i,j)
     $             + trans2(i,j,6)*vzq00(i,j)
              dvvy10 = trans2(i,j,4)*vxq10(i,j)+trans2(i,j,5)*vyq10(i,j)
     $             + trans2(i,j,6)*vzq10(i,j)
              dvvy01 = trans2(i,j,4)*vxq01(i,j)+trans2(i,j,5)*vyq01(i,j)
     $             + trans2(i,j,6)*vzq01(i,j)
              dvvy11 = trans2(i,j,4)*vxq11(i,j)+trans2(i,j,5)*vyq11(i,j)
     $             + trans2(i,j,6)*vzq11(i,j)
        dvyq00(i,j)=amax1(dvvy00,dvvy01,dvvy10,dvvy11)
     $          -amin1(dvvy00,dvvy01,dvvy10,dvvy11)
  750   continue
C
*       DO 800 J=1,MJP1
*       DO 800 I=1,MIP1
*       VXQ00(I,J) = (VXQ00(I,J)*SQ00(I,J) + VXQ01(I,J)*SQ01(i,J)
*     $    + VXQ10(I,J)*SQ10(I,J) + VXQ11(I,J)*SQ11(I,J))*SDENOM(I,J)
*       VYQ00(I,J) = (VYQ00(I,J)*SQ00(I,J) + VYQ01(I,J)*SQ01(i,J)
*     $     + VYQ10(I,J)*SQ10(I,J) + VYQ11(I,J)*SQ11(I,J))*SDENOM(I,J)
*       VZQ00(I,J) = (VZQ00(I,J)*SQ00(I,J) + VZQ01(I,J)*SQ01(i,J)
*     $     + VZQ10(I,J)*SQ10(I,J) + VZQ11(I,J)*SQ11(I,J))*SDENOM(I,J)
*  800  CONTINUE
C
        do 800 j=1,mjp1
        do 800 i=1,mip1
        vxq00(i,j) = 0.25*(vxq00(i,j)+vxq01(i,j)+vxq10(i,j)+vxq11(i,j))
        vyq00(i,j) = 0.25*(vyq00(i,j)+vyq01(i,j)+vyq10(i,j)+vyq11(i,j))
        vzq00(i,j) = 0.25*(vzq00(i,j)+vzq01(i,j)+vzq10(i,j)+vzq11(i,j))
  800   continue
C
*       call wrote('vxq00',vxq00,1,mip1,lip1,1,mjp1,ljp1)
*       call wrote('vyq00',vyq00,1,mip1,lip1,1,mjp1,ljp1)
*       call wrote('vzq00',vzq00,1,mip1,lip1,1,mjp1,ljp1)
C
        DO 810 J=1,MJP1
        DO 810 I=1,MIP1
        dvzz(i,j) = sqrt(
     $       dvxq00(i,j)**2+dvyq00(i,j)**2)
        VXCORN(I,J) = TRANS2(I,J,1)*VXQ00(I,J) +
     $          TRANS2(I,J,2)*VYQ00(I,J) + TRANS2(I,J,3)*VZQ00(I,J)
  810   VYCORN(I,J) = TRANS2(I,J,4)*VXQ00(I,J) +
     $          TRANS2(I,J,5)*VYQ00(I,J) + TRANS2(I,J,6)*VZQ00(I,J)
C
C
        DO 840 J=-1,MJ+2
        DO 840 I=1,MIP1
  840   BI0(I,J+2) = BI(I,J)
*        CALL BR2BI
C
        DO 841 J=1,MJP1
        DO 841 I=-1,MI+2
  841   BJ0(I+2,J) = BJ(I,J)
*        CALL BR2BJ
C
C
C       COMPUTE A HIGH ORDER B AT THE INTERFACE
C
        DO 850 J=1,MJP1
        DO 850 I=llow,MI+no2
        xnormj(i,j) = xnormj(i,j)*facej(i,j)
        ynormj(i,j) = ynormj(i,j)*facej(i,j)
        znormj(i,j) = znormj(i,j)*facej(i,j)
  850   BJQ(I,J) = BJ1(I,J)*FACEJ(I,J)
        CALL XBHI(FACEJ,FJJ,1)
        CALL XBHI(BJQ,BJJ,2)
        CALL xbhi(xnorMJ,xnQj,3)
        CALL xbhi(ynorMJ,ynQj,4)
        CALL xbhi(znorMJ,znQj,4)
C
        myunit = mype+30
*        write (myunit,*) '>>>>>>>>>>>>>>>>>>>>>>>>>>'
*        write (myunit,*) 'In Bfield, isweep= ',isweep
*        write (myunit,*)     '   level3= ',level3
*        write (myunit,*) ilow,ihigh,iactive,I_global_off,nl_i
*        write (myunit,*) jlow,jhigh,jactive,J_global_off,nl_j
*        write (myunit,*) klow,khigh,kactive,K_global_off,nl_k
*        if ( mype .eq. 3 .and. level3 .ge. 24 ) then
*           lli = lihigh-llow+1
*           llj = ljhigh-llow+1
*           call print2d(xnormi,'xni',lli,llj,lli,llj)
*           call print2d(facei,'facei',lli,llj,lli,llj)
*           call print2d(xnormj,'xnj',lli,llj,lli,llj)
*           call print2d(facej,'facej',lli,llj,lli,llj)
*           call print2d(xnormk,'xnk',lli,llj,lli,llj)
*           call print2d(facek,'facek',lli,llj,lli,llj)
*        endif
#ifdef PSCTEST
      myunit = mype+50
      if (level3 .ge. 23) then
      do j=llow,mj+no2
         DO I=1,MIP1
            write(myunit,977) level3,i,j,xnormi(i,j),ynormi(i,j),
     $           znormi(i,j),facei(i,j)
            call flush(myunit)
 977        format(3i6,4(1pe13.4))
         enddo
      enddo
      endif
#endif
        DO 851 J=llow,mj+no2
        DO 851 I=1,MIP1
        xnormi(i,j) = xnormi(i,j)*facei(i,j)
        ynormi(i,j) = ynormi(i,j)*facei(i,j)
        znormi(i,j) = znormi(i,j)*facei(i,j)
  851   BIQ(I,J) = BI1(I,J)*FACEI(I,J)
        CALL YBHI(FACEI,FII,1)
        CALL YBHI(BIQ,BII,2)
        CALL ybhi(xnorMI,xnQi,3)
        CALL ybhi(ynorMI,ynQi,4)
        CALL ybhi(znorMI,znQi,4)
C
c  high order interpolation doesn't work for the very distorted cells
c  near the front boundary
#ifdef BADGRID
       if ( isweep .eq. 1 ) then
c  compute the local equivalent for the j stop
#ifdef BADGRID_CORRECTED_ALLJ
c  Isolated upper-bound experiment: avoid the donor's hemisphere cap.
c  Keep the last global boundary edge out of the active-cell range.
        nj_stop = min(nj_bin-j_global_off,jactive+1)
#else
        nj_stop = min(nj_global/2-j_global_off,nj_bin-j_global_off,
     $         jactive+1)
#endif
c  compute local bounds for the i stop
c     global i runs from ni-ni_bin to nip1
c     on local processor this range is ni-ni_bin-i_global_off to 
c      nip1-ni_bin-i_global_off. Since we're coming in from the outer boundary
c      it suffices to compute the proper ni_start
        ni_start = max(1,min(nl_i+2,ni_global-ni_bin-i_global_off))
        
        do j=1,nj_stop
        do i=ni_start, nl_i+1
#ifdef BADGRID_CORRECTED_ALLJ
c  XBHI(F,FHF) interpolates F(i-1) and F(i) at FHF(i), since
c  FDUM(i+ORD2)=F(i). The donor F(i)+F(i+1) is shifted by one cell.
c  Area, area-weighted flux and area normals share this same midpoint.
           bjj(i,j) = 0.5*(bjq(i-1,j)+bjq(i,j))
           fjj(i,j) = 0.5*(facej(i-1,j)+facej(i,j))
           xnqj(i,j) = 0.5*(xnormj(i-1,j)+xnormj(i,j))
           ynqj(i,j) = 0.5*(ynormj(i-1,j)+ynormj(i,j))
           znqj(i,j) = 0.5*(znormj(i-1,j)+znormj(i,j))
#else
           bjj(i,j) = 0.5*(bjq(i+1,j)+bjq(i,j))
           fjj(i,j) = 0.5*(facej(i+1,j)+facej(i,j))
           xnqj(i,j) = 0.5*(xnormj(i+1,j)+xnormj(i,j))
           ynqj(i,j) = 0.5*(ynormj(i+1,j)+ynormj(i,j))
           znqj(i,j) = 0.5*(znormj(i+1,j)+znormj(i,j))
#endif
        enddo
        enddo
c
        endif
c
#endif
c
        DO 852 J=1,MJP1
        DO 852 I=1,MIP1
        xnQj(i,j) =xnQj(i,j)/fjj(i,j)
        ynQj(i,j) = ynQj(i,j)/fjj(i,j)
        znQj(i,j) = znQj(i,j)/fjj(i,j)
        BJJ(I,J) = BJJ(I,J)/FJJ(I,J)
        xnQi(i,j) = xnQi(i,j)/fii(i,j)
        ynQi(i,j) = ynQi(i,j)/fii(i,j)
        znQi(i,j) = znQi(i,j)/fii(i,j)
        BII(I,J) = BII(I,J)/FII(I,J)
 852  continue
c
C
        do j=1,mjp1
           DO i  = 1,mip1
              dxti = xnqi(i,j)*trans2(i,j,1) +
     $             ynqi(i,j)*trans2(i,j,2) +
     $             znqi(i,j)*trans2(i,j,3)
              dyti = xnqi(i,j)*trans2(i,j,4) +
     $             ynqi(i,j)*trans2(i,j,5) +
     $             znqi(i,j)*trans2(i,j,6)
              dinv = 1./sqrt(dxti**2+dyti**2)
              xnqi(i,j) = dinv*dxti
              ynqi(i,j) = dinv*dyti
           enddo
        enddo
C
        do j=1,mjp1
           DO i  = 1,mip1
              dxtj = xnqj(i,j)*trans2(i,j,1) +
     $             ynqj(i,j)*trans2(i,j,2) +
     $             znqj(i,j)*trans2(i,j,3)
              dytj = xnqj(i,j)*trans2(i,j,4) +
     $             ynqj(i,j)*trans2(i,j,5) +
     $             znqj(i,j)*trans2(i,j,6)
              dinv = 1./sqrt(dxtj**2+dytj**2)
              xnqj(i,j) = dinv*dxtj
              ynqj(i,j) = dinv*dytj
           enddo
        enddo
C
        do 880 j=1,mjp1
        do 880 i=1,mip1
        BXQ0(I,J) = TRANS2(I,J,1)*BXC0(I,J) +
     $          TRANS2(I,J,2)*BYC0(I,J) + TRANS2(I,J,3)*BZC0(I,J)
        BYQ0(I,J) = TRANS2(I,J,4)*BXC0(I,J) +
     $          TRANS2(I,J,5)*BYC0(I,J) + TRANS2(I,J,6)*BZC0(I,J)
  880   continue
C
        CALL JCLIP(BII,BI0)
        CALL ICLIP(BJJ,BJ0)
C
        CALL JPDM(BII,BI0,BIRT,BILFT)
        CALL IPDM(BJJ,BJ0,BJRT,BJLFT)
*       if (isweep .eq. 2 .and.(level3 .eq. 4))then
*       call wrote('bi0',bi0,1,mip1,lip1,1,mjp4,ljp4)
*       call wrote('bii',bii,1,mip1,lip1,1,mjp1,ljp1)
*       call wrote('bilft',bilft,1,mip1,lip1,1,mjp1,ljp1)
*       call wrote('birt',birt,1,mip1,lip1,1,mjp1,ljp1)
*       call wrote('bj0',bj0,1,mip4,lip4,1,mjp1,ljp1)
*       call wrote('BJJ',BJJ,1,mip1,lip1,1,mjp1,ljp1)
*       call wrote('bjlft',bjlft,1,mip1,lip1,1,mjp1,ljp1)
*       call wrote('bjrt',bjrt,1,mip1,lip1,1,mjp1,ljp1)
*       endif
C
        DO 900 J=1,MJP1
        DO 900 I=1,MIP1
        BII(I,J) = ((SQ00(I,J)+SQ10(I,J))*BILFT(I,J)+
     $          (SQ11(I,J)+SQ01(I,J))*BIRT(I,J))*SDENOM(I,J)
  900   BJJ(I,J) = ((SQ00(I,J)+SQ01(I,J))*BJLFT(I,J)+
     $          (SQ11(I,J)+SQ10(I,J))*BJRT(I,J))*SDENOM(I,J)
C
        DO 910 J=1,MJP1
        DO 910 I=1,MIP1
        det(i,j) = 1./(xnQi(i,j)*ynQj(i,j)-xnQj(i,j)*ynQi(i,j))
        BXCORN(I,J) = det(i,j)*(ynQj(i,j)*BII(I,J) - YnQI(i,j)*BJJ(I,J))
     $          + bxq0(i,j)
  910   BYCORN(I,J) = det(i,j)*(-xnQj(i,j)*BII(I,J) +xnQi(I,J)*BJJ(I,J))
     $          + byq0(i,j)
C
C
#ifdef BDIFFUSE
        casq = ca**2
        do j=1,mjp1
           do i=1,mip1
           rhobar = 4./(rho(i-1,j-1)+rho(i,j-1)+rho(i-1,j)+rho(i,j))
           vasq = (bxcorn(i,j)**2+bycorn(i,j)**2)*pi4inv*rhobar
              dvzz(i,j) = dvzz(i,j) + sqrt( vasq*casq/(casq+vasq))
           enddo
        enddo
#endif
C            
C
        DO 1000 J=1,MJP1
        DO 1000 I=1,MIP1
 1000   E(I,J) = (VXCORN(I,J)*BYCORN(I,J) - VYCORN(I,J)*BXCORN(I,J)
     $     + 0.5*dvzz(i,j)*(birt(i,j)-bilft(i,j)-bjrt(i,j)+bjlft(i,j))
     $          )*EDGEK(I,J)
* 1000      e(i,j) = 0.
C
c
*        CALL EFIX(E,TRANS2,EDGEK)
c   insert in-line here -- kind of klugey, but it should work
        if ( isweep .eq. 1 ) then
c
         if (ilow_bound) then
           do j=1,mjp1
              e(1,j) = e2(1,j)
           enddo
         endif
         if (ihigh_bound) then
           do j=1,mjp1
              e(mip1,j) = (1.-bound_mix(j))*e2(mip1,j) +
     $             bound_mix(j)*e(mip1,j)
*     $         + 0.25*abs(vxbar(mip1,j))*edgek(mip1,j)*
*     $         (bi0(mip1,j+2)-bi0(mip1,j+1)-bj0(mi+3,j)+bj0(mi+2,j))
           enddo
         endif
c    make sure bj on axis is zero
         if ( jlow_bound ) then 
            do i=1,nl_i+1
               e(i,1) = 0.0
            enddo
         endif
         if ( jhigh_bound ) then 
            do i=1,nl_i+1
               e(i,nl_j+1) = 0.0
            enddo
         endif
c
        else if ( isweep .eq. 2 ) then
c
           if ( jlow_bound ) then
           esum1 = e(1,1)
           do 2701 j=2,mj
 2701         esum1 = esum1 + e(1,j)
              fact = 1./float(mj)
              esum1 = esum1*fact
              do 2702 j=1,mjp1
 2702            e(1,j) = esum1
           endif
           if ( jhigh_bound ) then
           esum2 = e(mip1,1)
           do 2801 j=2,mj
 2801         esum2 = esum2 + e(mip1,j)
              fact = 1./float(mj)
              esum2 = esum2*fact
              do 2802 j=1,mjp1
 2802            e(mIp1,j) = esum2
            endif
C   and let's get e constant at the seam in the k direction
              do 2705 i=1,mip1
 2705            e(i,mjp1) = 0.5*(e(i,1)+e(i,mjp1))
               do 2706 i=1,mip1
 2706             e(i,1) = e(i,mjp1)
c
         else if ( isweep .eq. 3 ) then
c
c  more upper and lower global boundaries
c
          if (ilow_bound) then
            do i=1,mip1
               e(i,1) = e2(i,1)
            enddo
          endif
          if (ihigh_bound) then
            do i=1,mip1
               e(i,mjp1) = (1.-bound_mix(i))*e2(i,mjp1) +
     $              bound_mix(i)*e(i,mjp1)
*     $           + 0.25*abs(vxbar(i,mjp1))*edgek(i,mjp1)*
*     $           (bi0(i,mj+3)-bi0(i,mj+2)-bj0(i+3,mjp1)+bj0(i+2,mjp1))
            enddo
          endif
c   same at seam
            do j=1,mjp1
               edd = 0.5*(e(mip1,j)+e(1,j))
               e(1,j) = edd
               e(mip1,j) = edd
            enddo
c            
c
         endif
C
           
C
C
        DO 1020 J=1,MJ
        DO 1020 I=1,MIP1
 1020   BI2(I,J) = BI2(I,J) + DT*(E(I,J+1)-E(I,J))
C
*        write (6,*) (bj2(mi,j),j=1,mjp1)
        DO 1040 J=1,MJP1
        DO 1040 I=1,MI
 1040   BJ2(I,J) = BJ2(I,J) - DT*(E(I+1,J)-E(I,J))
*        write (6,*) (bj2(mi,j),j=1,mjp1)
C
C
        
        if ( isweep .eq. 2 ) then
         if ( jlow_bound ) then
           d1 = bj2(1,nk2+1)+bj2(1,1)
           d2 = bj2(1,nk4+1)+bj2(1,nk2+1+nk4)
           efixit1 = axfix*0.25*(d1+d2)
           do 1141 j=1,nl_k+1
             bj2(1,j) = bj2(1,j) - efixit1
             e(1,j) = e(1,j)  -efixit1/dt
 1141   continue
         endif
         if ( jhigh_bound ) then
           d1 = bj2(nl_j,nk2+1)+bj2(nl_j,1)
           d2 = bj2(nl_j,nk4+1)+bj2(nl_j,nk2+1+nk4)
           efixit2 = axfix*0.25*(d1+d2)
           do 1142 j=1,nl_k+1
             bj2(nl_j,j) = bj2(nl_j,j) - efixit2
             e(nl_j+1,j) = e(nl_j+1,j) + efixit2/dt
 1142     continue
          endif
        endif
c
c  put electric fields in var common block for later processing
c
        do j=1,mjp1
           do i=1,mip1
              e2(i,j) = e(i,j)
           enddo
        enddo
*        call eload(e)
*       if (isweep .eq. 2 .and.(level3 .eq. 4))then
*       WRITE(9,*) '******* ISWEEP, LEVEL3 *********',ISWEEP,LEVEL3
*       call wrote('EK',E,1,MIP1,LIP1,1,MJP1,MJP1)
*       CALL WROTE('VXCORN',VXCORN,1,MIP1,LIP1,1,MJP1,MJP1)
*       CALL WROTE('VYCORN',VYCORN,1,MIP1,LIP1,1,MJP1,MJP1)
*       CALL WROTE('BXCORN',BXCORN,1,MIP1,LIP1,1,MJP1,MJP1)
*       CALL WROTE('BYCORN',BYCORN,1,MIP1,LIP1,1,MJP1,MJP1)
*       endif
C
!
! ... End SUBROUTINE bfieLD.............................................
!
        RETURN
        END
