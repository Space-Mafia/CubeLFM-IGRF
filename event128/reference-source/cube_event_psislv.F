      SUBROUTINE CUBE_EVENT_PSISLV(err,north,sp_event,sh_event)
!
! ... Local variables ..................................................
!
#include "ion-param.inc"
      REAL LFM_PLANET_RION,LFM_PLANET_FAC
#include "ion_P++_pointers.inc"
#include "ion_fortran_P++_arrays.inc"
#include "iongeo.inc"
#include "ionconst.inc"
#include "volts.inc"
#include "ion-var.inc"
#include "gcoef.inc"
#include "boundx.inc"
      common /ionscrach/
     $    curnt1(njp1,nk2p3),SUMA(NJP1,NK2P3),
     $          ERRA(NJP1,NK2P3),apot1(njp1,nk2p3,9),
     $          ab1(nk2p3),ab2(nk2p3)
      common /sigpeuv/sigpeuv(njp1,nk2p3)
      REAL sp_event(njp1,nk2p3),sh_event(njp1,nk2p3)
      LOGICAL first
      DIMENSION LBOUN(4)
      SAVE first,lboun,itmaxx,error
      DATA (LBOUN(I),I=1,4)/-1,-1,0,0/
      DATA first /.true./
      DATA ITMAXX,ERROR/400,5.e-3/
!
! ... Parameter variables ..............................................
!
      logical north
!
! ... Begin ............................................................
!
#ifdef DEBUG_MODE_ON
      write(*,*) "DEBUG:  in ionosphere.F::PSISLV(...)"
#endif
      do loopy = 1,looptries
	itmaxx = itmax
	if (first) then
	  do k=1,nk2p3
	   do j=1,njp1
	     psi(j,k) = 1.e-4
	   enddo
          enddo

	first = .false.
	endif

        pi = atan(1.0)*4.0
        pi2 = pi/2.0
        ang65 = pi/180.0*65.0
        ang100 = pi*5./9.
        pref = 2.0*250.0**(-0.666666)
        href = 1./(1.8*sqrt(250.0))
        rad2deg = 180.0/pi

        shall = 1.8*sqrt(s107)
        speder = 0.5*s107**0.666666
        pedslope = 0.24*pref*speder*rad2deg
        hallslope = 0.27*href*shall*rad2deg
        sigmap65 = speder*cos(ang65)**0.666666
        sigmah65 = shall*cos(ang65)
        sigmap100 = sigmap65 - (ang100-ang65)*pedslope
        pedslope2 = 0.13*pref*speder*rad2deg
        if ( north ) then
           tilt0 = tilt_angle
        else
           tilt0 = -tilt_angle
        endif

	if ( ionmodel ) then
c  ionospheric conductivity model

c         Put in the day-night gradient in Mhos
      do k=1,nk2p3
         do j=1,njp1
            zenith = pi2 - (x2ion(j,k)+tilt0)
            if( zenith .le. ang65) then
               sigmah(j,k) = shall*cos(zenith)
               sigmap(j,k) = speder*cos(zenith)**0.666666
            elseif (zenith .le. ang100) then
               sigmap(j,k) = sigmap65 - pedslope*(zenith-ang65)
               sigmah(j,k) = sigmah65 - hallslope*(zenith-ang65)
            else
               sigmap(j,k) = sigmap100 - pedslope2*(zenith-ang100)
               sigmah(j,k) = sigmah65 - hallslope*(zenith-ang65)
            endif

            sigpeuv(j,k) = sigmap(j,k)
         enddo
      enddo
      
c
      call pardy(curnt1,suma,dflux,denergy,p_hardy,c_hardy,h_hardy)
C note curnt1 is used as temporary storage
      do k=1,nk2p3
         do j=1,njp1
            sigmap(j,k) = sqrt(
     $           sigmap(j,k)**2 + curnt1(j,k)**2)
            sigmah(j,k) = sqrt(
     $           sigmah(j,k)**2 + suma(j,k)**2)

            sigmap(j,k) = amax1(pedmin,sigmap(j,k))
            sigmah(j,k) = -amin1(amax1(hallmin,sigmah(j,k)),
     $           sigmap(j,k)*sigma_ratio)
         enddo
      enddo
c
	else
c  constant pedersen conductivity
c
         do k=1,nk2p3
            do j=1,njp1
               sigmap(j,k) = sp_event(j,k)
               sigmah(j,k) = -sh_event(j,k)
            enddo
         enddo
c
      endif
c
#ifdef DEBUG
      write(9,*) ' PSISLV: first=',first
#endif
C
      DO 22 J=1,Nk2p3
      DO 22 I=1,Njp1
      ERRA(I,J) = 0.
   22 continue
C
C   experimental renormalization
      sumarea = 0.
      sumj = 0.
      do 23 j=2,nk2+2
      do 23 i=1,njp1
         sumarea = sumarea + aj(i,j,9)
         sumj = sumj + curnt(I,j)*aj(i,j,9)
 23   continue
      curav = sumj/sumarea
*     do 24 j=2,nk2+2
*     do 24 i=1,njp1
*     curnt(i,j) = curnt(i,j)-curav
*24   continue
C
         do j=1,nk2p3
         do i=1,njp1
            curnt1(i,j) = 0.
            do kl=1,9
               apot1(i,j,kl) = 0.
            enddo
         enddo
         enddo
      do 30 j=2,nK2P2
      do 30 i=2,nj
         CURNT1(I,J)=Aj(I,J,1)*CURNT(I-1,J-1)+Aj(I,J,2)*CURNT(I,J-1)+
     1        Aj(I,J,3)*CURNT(I+1,J-1)+Aj(I,J,4)*CURNT(I+1,J)+
     2        Aj(I,J,5)*CURNT(I+1,J+1)+Aj(I,J,6)*CURNT(I,J+1)+
     3        Aj(I,J,7)*CURNT(I-1,J+1)+Aj(I,J,8)*CURNT(I-1,J)+
     4        aj(i,j,9)*curnt(i,j)
      do 30 kl=1,9
         apot1(I,j,kl) = Apeder(I,j,kl,1)*SIGMAP(I-1,J-1)+
     $        Apeder(I,j,kl,2)*SIGMAP(I,J-1)+
     1        Apeder(I,j,kl,3)*SIGMAP(I+1,J-1)+
     $        Apeder(I,j,kl,4)*SIGMAP(I+1,J)+
     2        Apeder(I,j,kl,5)*SIGMAP(I+1,J+1)+
     $        Apeder(I,j,kl,6)*SIGMAP(I,J+1)+
     3        Apeder(I,j,kl,7)*SIGMAP(I-1,J+1)+
     $        Apeder(I,j,kl,8)*SIGMAP(I-1,J)+
     4        apeder(i,j,kl,9)*sigmap(i,j)
     $      + Ahall(I,j,kl,1)*SIGMAH(I-1,J-1)+
     $        Ahall(I,j,kl,2)*SIGMAH(I,J-1)+
     1        Ahall(I,j,kl,3)*SIGMAH(I+1,J-1)+
     $        Ahall(I,j,kl,4)*SIGMAH(I+1,J)+
     2        Ahall(I,j,kl,5)*SIGMAH(I+1,J+1)+
     $        Ahall(I,j,kl,6)*SIGMAH(I,J+1)+
     3        Ahall(I,j,kl,7)*SIGMAH(I-1,J+1)+
     $        Ahall(I,j,kl,8)*SIGMAH(I-1,J)+
     4        ahall(i,j,kl,9)*sigmah(i,j)
 30   continue
C
C     SET UP THE CONVERGENT POINTS
      DO 40 j=2,NK2P2
         CURNT1(1,J)=AJ(1,J,2)*CURNT(1,J-1)+AJ(1,J,3)*CURNT(2,J-1)
     $        +  AJ(1,J,4)*CURNT(2,J) + AJ(1,J,5)*CURNT(2,J+1) +
     $        AJ(1,J,6)*CURNT(1,J+1) + AJ(1,J,9)*CURNT(1,J)
         CURNT1(NJP1,J) = AJ(NJP1,J,1)*CURNT(NJ,J-1) +
     $        AJ(NJP1,J,2)*CURNT(NJP1,J-1)
     $        +  AJ(NJP1,J,6)*CURNT(NJP1,J+1) +
     $        AJ(NJP1,J,7)*CURNT(NJ,J+1) +
     $        AJ(NJP1,J,8)*CURNT(NJ,J)
     $        + AJ(NJP1,J,9)*CURNT(NJP1,J)
      DO 40 KL=1,9
         apot1(1,j,kl) = Apeder(1,j,kl,2)*SIGMAP(1,j-1)+
     1        Apeder(1,j,kl,3)*SIGMAP(2,J-1)+
     $        Apeder(1,j,kl,4)*SIGMAP(2,J)+
     2        Apeder(1,j,kl,5)*SIGMAP(2,J+1)+
     $        Apeder(1,j,kl,6)*SIGMAP(1,j+1)+
     4        apeder(1,j,kl,9)*sigmap(1,j)
     $        + Ahall(1,j,kl,2)*SIGMAH(1,J-1)+
     1        Ahall(1,j,kl,3)*SIGMAH(2,J-1)+
     $        Ahall(1,j,kl,4)*SIGMAH(2,J)+
     2        Ahall(1,j,kl,5)*SIGMAH(2,J+1)+
     $        Ahall(1,j,kl,6)*SIGMAH(1,J+1)+
     4        ahall(1,j,kl,9)*sigmah(1,j)
         apot1(Njp1,j,kl) = Apeder(Njp1,j,kl,1)*SIGMAP(Nj,J-1)+
     $        Apeder(Njp1,j,kl,2)*SIGMAP(Njp1,j-1)+
     $        Apeder(Njp1,j,kl,6)*SIGMAP(Njp1,j+1)+
     3        Apeder(Njp1,j,kl,7)*SIGMAP(Nj,J+1)+
     $        Apeder(Njp1,j,kl,8)*SIGMAP(Nj,J)+
     4        apeder(njp1,j,kl,9)*sigmap(njp1,j)
     $      + Ahall(njp1,j,kl,1)*SIGMAH(nj,J-1)+
     $        Ahall(njp1,j,kl,2)*SIGMAH(njp1,J-1)+
     $        Ahall(njp1,j,kl,6)*SIGMAH(njp1,J+1)+
     3        Ahall(njp1,j,kl,7)*SIGMAH(nj,J+1)+
     $        Ahall(njp1,j,kl,8)*SIGMAH(nj,J)+
     4        ahall(njp1,j,kl,9)*sigmah(njp1,j)
 40   continue
C
      do 41 j=3,nk2p2
         apot1(1,2,9) = apot1(1,2,9)+apot1(1,j,9)
         apot1(1,j,1) = apot1(1,j,4)+apot1(1,j-1,5)
         apot1(njp1,2,9) = apot1(njp1,2,9)+apot1(njp1,j,9)
         apot1(njp1,j,4) = apot1(njp1,j,8)+apot1(njp1,j-1,7)
 41   continue
      apot1(1,2,1) = apot1(1,2,4)
      apot1(njp1,2,4) = apot1(njp1,2,8)
      source00 = 0.
      source11 = 0.
      do 42 j=2,nk2p2
         source00 = source00 + curnt1(1,j)
         source11 = source11 + curnt1(njp1,j)
         ab1(j) = - apot1(1,j,1)/apot1(1,2,9)
         ab2(j) = - apot1(njp1,j,4)/apot1(njp1,2,9)
 42   continue
C
      s1 = source00/apot1(1,2,9)
      s2 = source11/apot1(njp1,2,9)
#ifdef DEBUG
      write(9,*) ' PSISLV: source00,source11=',source00,source11
      write (9,*) 'stencil at 1',apot1(1,2,9),(apot1(1,j,1),j=2,nk2p2)
      write (9,*) 'stencil at nj',
     $     apot1(njp1,2,9),(apot1(njp1,j,4),j=2,nk2p2)
C
      write (9,*) 'sources',s1,s2,s1,s2
      write(9,*) 'ab1',ab1
      write(9,*) 'ab2',ab2
#endif /*DEBUG*/
#ifdef PICTURES
      CALL PICTUR(APOT1,'A  ',1,NJP1,1,NK2P3,1,9,0,2,3,1,NJP1,
     $    NK2P3,9,NJP1*NK2P3)
      CALL PRINT2D(CURNT1,'CURN',NJP1,NK2P3,NJP1,NK2P3)
#endif /*PICTURES*/
C         
      CALL LFM_MI_Z2ICG(PSI,APOT1(1,1,8),APOT1(1,1,2),APOT1(1,1,9),
     $     APOT1(1,1,6),APOT1(1,1,4),APOT1(1,1,1),APOT1(1,1,7),
     $     APOT1(1,1,3),APOT1(1,1,5),curnt1,0,
     $     njp1,nk2p3,njp1,nk2p3,lboun,ITMAXx,ITER,ERRmax,SER,
     $     ab1,s1,ab2,s2)
#ifdef DEBUG
      write (9,*) 'psi output'
      call print2d(psi,'psi',njp1,nk2p3,njp1,nk2p3)
#endif
C
#ifdef DEBUG_MODE_ON
      write (6,*) 'ionosphere.F:  north:',
     $     north, '  loopy: ', loopy, 
     $     '   ser: ', ser
#endif
*	if ( ser .gt. errmax ) then
*	write (6,*) 'exceeded error margin in psislv',ser,errmax,
*     $	  itmax, iter
*	endif
*	if ( ser .gt. error ) then
*	write (6,*) 'exceeded stop error in psislv, stopping run',
*     $	   ser,errmax,error,itmax,iter
*	endif
c              
*      CALL Z2ICG(PSI(1,2),APOT1(1,2,8),APOT1(1,2,2),APOT1(1,2,9),
*     $     APOT1(1,2,6),APOT1(1,2,4),APOT1(1,2,1),APOT1(1,2,7),
*     $     APOT1(1,2,3),APOT1(1,2,5),curnt1(1,2),0,
*     $     njp1,nk2p1,njp1,nk2p2,lboun,ITMAX,ITER,ERROR,SER,
*     $     ab1(2),s1,ab2(2),s2)
C
#ifdef PICTURES
      call print2d(psi,'psi',njp1,nk2p3,njp1,nk2p3)
#endif /*PICTURES*/
C
      if ( ser .lt. errmax ) return
      enddo   !   end of loopy loop surrounding whole subroutine
      if ( iter .ge. itmaxx .or. ser .ge. error )
     $     WRITE(9,999) ITER,ITMAX,SER,ERROR
      RETURN
 999  FORMAT(' ******* POTENTIAL SOLUTION ******* ',/,
     1     '    AFTER   ',I5,'  ITERATIONS WITH A MAX  OF   ',
     $     I5,'  ITERATIONS  ',/,
     2     '    ERROR IN POTENTIAL =  ',E10.4,'  WITH A MAX OF ',
     $     E10.4)
!
! ... End SUBROUTINE PSISLV(err,north)..................................
!
      END
