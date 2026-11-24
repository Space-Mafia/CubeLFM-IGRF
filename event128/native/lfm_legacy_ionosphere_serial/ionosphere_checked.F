!> ionosphere.F : FORTRAN ionosphere solver
!! Excerpt from 2004
!!   -  Lyon, Fedder, Mobarry, "The Lyon_Fedder-Mobarry (LFM) global MHD magnetospheric simulation code"
!!
!!     The ionosphere is one of the largest, if not the largest energy
!!     recipient in the magnetosphere.  As far as a global MHD
!!     calculation of the magnetosphere is concerned, it is a boundary
!!     condition, but one that may have a profound effect on
!!     magnetospheric structure.  Thus, we dcided to implement a simple
!!     model fo the polar ionospheres which allowed for self-consistent
!!     electrodynamic coupling with the magnetosphere.
!!     
!! This FORTRAN code solves the ionosphere model referenced in the
!! aforementioned paper.
!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!> aveion(p,a)
!! Interpolate cell quantitites to scalar quantities above the ionosphere
!!
      subroutine aveion(p,a)
!
! ... Local variables ..................................................
!
#include "ion-param.inc"
      REAL LFM_PLANET_RION,LFM_PLANET_FAC
#include "boundx.inc"
!
! ... Parameter variables ..............................................
!
      real p(ni,nj,nk),a(njp1,nkp1)
!
! ... Begin ............................................................
!
#ifdef DEBUG_MODE_ON
      write(*,*) "DEBUG:  in ionosphere.F::aveion(p,a)"
#endif
      do k=2,nkp1
         do j=2,njp1
            a(j,k) = 0.25*(p(2,j,k)+p(2,j-1,k)
     $           +p(2,j,k-1)+p(2,j-1,k-1))
         enddo
      enddo
      do j=2,njp1
         a(j,1) = 0.25*(p(2,j,1)+p(2,j-1,1)
     $        +p(2,j,nk)+p(2,j-1,nk))
      enddo
      do k=2,nkp1
         a(1,k) = 0.25*(p(2,1,k)+p(2,nj,k)
     $        +p(2,1,k-1)+p(2,nj,k-1))
      enddo
      a(1,1) = 0.25*(p(2,1,1)+p(2,nj,1)
     $     +p(2,1,nk)+p(2,nj,nk))
!
! ... End subroutine aveion(p,a)........................................
!
      return
      end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!> COEF1 - Setup the coefficient stencil
!!  
!!    here's how the coefficient stencil works out. The
!!    center point, marked - 9 - is point (i,j) in the
!!    grid. so a(i,j,9) refers to a weight applied to
!!    central point of the stencil. a(i,j,1) applies
!!    to f(i-1,j-1) and so on. The numbers in the
!!    cells refer the sub-elements - e - that span
!!    that particular cell.
!!    
!!    7 ----- 6 ----- 5
!!    |       |       |
!!    |   2   |   1   |
!!    |       |       |
!!    8 ----- 9 ----- 4
!!    |       |       |
!!    |   4   |   3   |
!!    |       |       |
!!    1 ----- 2 ----- 3
!!    
      SUBROUTINE COEF1
!
! ... Local variables ..................................................
!
#include "ion-param.inc"
      REAL LFM_PLANET_RION,LFM_PLANET_FAC
#include "ion_P++_pointers.inc"
#include "ion_fortran_P++_arrays.inc"
#include "iongeo.inc"
#include "gcoef.inc"
#include "ionconst.inc"
      dimension gpt(ngpts),gwt(ngpts)
      data (gpt(i),i=1,ngpts)/0.069431844,0.160018957,
     $     0.839981043,0.930568155/
      data (gwt(i),i=1,ngpts)/0.173927422,0.326072577,
     $     0.326072577,0.173927422/
      dimension apot1(njp1,nk2p3,9),dom(njp1,nk2p3)
      dimension factor(njp1,nk2p3,3)
      dimension kstencl(4,4)
      SAVE GPT,GWT,KSTENCL
      data ((kstencl(i,j),i=1,4),j=1,4) /
     $     9,      8,      2,      1,
     $     4,      9,      3,      2,
     $     6,      7,      9,      8,
     $     5,      6,      4,      9 /
!
! ... Parameter variables ..............................................
!
! ... Begin ............................................................
!
#ifdef DEBUG_MODE_ON
      write(*,*) "DEBUG:  in ionosphere.F::COEF1(...)"
#endif
C For each point in a 9-point stencil
      do 200 k=1,9
C nk2p3 = (nk/2)+3
      do 200 j=1,nk2p3
C njp1 = nj + 1
      do 200 i=1,njp1
         theta = sqrt(x2ion(i,j)**2+y2ion(i,j)**2)
         sinth = sin(theta)
         if ( theta .lt. 1.e-2 ) then
            gxx(i,j) = 1.
            gyy(i,j) = 1.
            gxy(i,j) = 0.
            gh(i,j) = 1.
            gj(i,j) = rionm**2
         else
            stt = sinth/theta
            costh = cos(theta)
            gxx(i,j) = stt*(
     $           x2ion(i,j)**2*(1.+3.*costh**2)/(2.*theta*costh)**2+
     $           (y2ion(i,j)/sinth)**2 )
            gyy(i,j) = stt*(
     $           y2ion(i,j)**2*(1.+3.*costh**2)/(2.*theta*costh)**2+
     $           (x2ion(i,j)/sinth)**2 )
            gxy(i,j) = stt*x2ion(i,j)*y2ion(i,j)*(
     $           (1.+3.*costh**2)/(2.*theta*costh)**2 -
     $           1./sinth**2 )
            gh(i,j) = sqrt(1.+3.*costh**2)/(2.*costh)
            gj(i,j) = stt*rionm**2
         endif
         aj(i,j,k) = 0.
         do 200 kl=1,9
            apeder(i,j,k,kl) = 0.
            ahall(i,j,k,kl) = 0.
 200     continue
#ifdef DEBUG
         call print2d(gxx,'gxx',njp1,nk2p3,njp1,nk2p3)
         call print2d(gyy,'gyy',njp1,nk2p3,njp1,nk2p3)
         call print2d(gxy,'gxy',njp1,nk2p3,njp1,nk2p3)
         call print2d(gh,'gh',njp1,nk2p3,njp1,nk2p3)
         call print2d(gj,'gj',njp1,nk2p3,njp1,nk2p3)
#endif
C     
         do 300 m=1,ngpts
         do 300 l=1,ngpts
            e(l,m,1,1) = (1.-gpt(l))*(1.-gpt(m))
            e(l,m,2,1) =     gpt(l) *(1.-gpt(m))
            e(l,m,3,1) = (1.-gpt(l))*    gpt(m)
            e(l,m,4,1) =     gpt(l) *  gpt(m)
            eeta(l,m,1,1) = - (1.-gpt(m))
            eeta(l,m,2,1) = (1.-gpt(m))
            eeta(l,m,3,1) = -gpt(m)
            eeta(l,m,4,1) = gpt(m)
            echi(l,m,1,1) = -(1.-gpt(l))
            echi(l,m,2,1) = -gpt(l)
            echi(l,m,3,1) = (1.-gpt(l))
            echi(l,m,4,1) = gpt(l)
            wtij(l,m) = gwt(l)*gwt(m)
 300     continue
         do 301 m=1,ngpts
         do 301 l=1,ngpts
            e(l,m,1,2) = e(l,m,1,1)+e(l,m,3,1)
            echi(l,m,1,2) = echi(l,m,1,1)+echi(l,m,3,1)
            eeta(l,m,1,2) = eeta(l,m,1,1)+eeta(l,m,3,1)
            e(l,m,3,2) = 0.
            echi(l,m,3,2) = 0.
            eeta(l,m,3,2) = 0.
            e(l,m,2,2) = e(l,m,2,1)
            echi(l,m,2,2) = echi(l,m,2,1)
            eeta(l,m,2,2) = eeta(l,m,2,1)
            e(l,m,4,2) = e(l,m,4,1)
            echi(l,m,4,2) = echi(l,m,4,1)
            eeta(l,m,4,2) = eeta(l,m,4,1)
C     and now the pole at the back
            e(l,m,2,3) = e(l,m,2,1)+e(l,m,4,1)
            echi(l,m,2,3) = echi(l,m,2,1)+echi(l,m,4,1)
            eeta(l,m,2,3) = eeta(l,m,2,1)+eeta(l,m,4,1)
            e(l,m,4,3) = 0.
            echi(l,m,4,3) = 0.
            eeta(l,m,4,3) = 0.
            e(l,m,1,3) = e(l,m,1,1)
            echi(l,m,1,3) = echi(l,m,1,1)
            eeta(l,m,1,3) = eeta(l,m,1,1)
            e(l,m,3,3) = e(l,m,3,1)
            echi(l,m,3,3) = echi(l,m,3,1)
            eeta(l,m,3,3) = eeta(l,m,3,1)
 301     continue
C     
*     call print3('e(1)',e,ngpts,ngpts,4)
*     call print3('e(2)',e(1,1,1,2),ngpts,ngpts,4)
*     call print3('e(3)',e(1,1,1,3),ngpts,ngpts,4)
*     call print3('eeta(1)',eeta,ngpts,ngpts,4)
*     call print3('eeta(2)',eeta(1,1,1,2),ngpts,ngpts,4)
*     call print3('eeta(3)',eeta(1,1,1,3),ngpts,ngpts,4)
*     call print3('echi(1)',echi,ngpts,ngpts,4)
*     call print3('echi(2)',echi(1,1,1,2),ngpts,ngpts,4)
*     call print3('echi(3)',echi(1,1,1,3),ngpts,ngpts,4)
C
         do 1000 i=1,njp1-1
         do 1000 j=1,nk2p3-1
C
*     if ( i .ge. 6 .and. i .le. 7 .and. j .ge. 9 .and. j .le. 10 )
*     $         then
*     write (9,*) '******************************************'
*     write (9,*) '******************************************'
*     write (9,*) '******************************************'
*     write (9,*) '******************************************'
*     write (9,*) ' i   j  ',i,j
*     write (9,*) '******************************************'
*     write (9,*) 'u ',((x2ion(lli,llj),lli=i,i+1),llj=j,j+1)
*     write (9,*) 'v ',((y2ion(lli,llj),lli=i,i+1),llj=j,j+1)
*     endif
C

            do 400 ii=1,4
            do 400 jj=1,4
               cjwt(ii,jj) = 0.
               do 400 kl=1,4
                  cpwt(ii,jj,kl) = 0.
                  chwt(ii,jj,kl) = 0.
 400           continue
C
C     need to check the Jacobian matrix for ill condition
C     
               ueta = 0.5*((x2ion(i+1,j)+x2ion(i+1,j+1))-
     $              (x2ion(i,j)+x2ion(i,j+1)))
               veta = 0.5*((y2ion(i+1,j)+y2ion(i+1,j+1))-
     $              (y2ion(i,j)+y2ion(i,j+1)))
               uchi = 0.5*((x2ion(i,j+1)+x2ion(i+1,j+1))-
     $              (x2ion(i,j)+x2ion(i+1,j)))
               vchi = 0.5*((y2ion(i,j+1)+y2ion(i+1,j+1))-
     $              (y2ion(i,j)+y2ion(i+1,j)))
C     
               alpha1 = SIGN(abs(vchi),-vchi*veta)
     $              /sqrt(vchi**2+veta**2)
               beta1 = abs(veta)/sqrt(vchi**2+veta**2)
               condit1 = ((alpha1*uchi-beta1*ueta)
     $              /(alpha1*vchi-beta1*veta+1.e-15))**2
               alpha2 = sign(abs(ueta),-ueta*uchi)
     $              /sqrt(uchi**2+ueta**2)
               beta2 = abs(uchi)/sqrt(uchi**2+ueta**2)
               condit2 = ((alpha2*veta+beta2*vchi)
     $              /(alpha2*ueta+beta2*uchi+1.e-15))**2
               factor(i,j,1) = condit1
               factor(i,j,2) = condit2
               factor(i,j,3) = (ueta*uchi+veta*vchi)/sqrt(
     $              (ueta**2+veta**2)*(uchi**2+vchi**2))
               if ( condit1 .lt. condit2 ) then
                  if ( condit1 .gt. 0.5 ) then
                     write (9,*	)
     $                    'warning 1, possible ill condition at ',i,j,
     $                    '  value ',condit1
                  endif
                  alpha = alpha1
                  beta = beta1
               else
                  if ( condit2 .gt. 0.5 ) then
                     write (9,*)
     $                    'warning 2, possible ill condition at ',i,j,
     $                    '  value ',condit2
                  endif
                  alpha = alpha2
                  beta = beta2
               endif
C     
               ueta0 = ueta
               uchi0 = uchi
               veta0 = veta
               vchi0 = vchi
               ueta = (alpha*ueta0 + beta*uchi0)
               uchi = (alpha*uchi0 - beta*ueta0)
               veta = (alpha*veta0 + beta*vchi0)
               vchi = (alpha*vchi0 - beta*veta0)
C     
*     if ( i .ge. 6 .and. i .le. 7 .and. j .ge. 9 .and. j .le. 10 )
*     $         then
*     write (9,*) 'zero ueta,veta,uchi,vchi',
*     $         ueta0,veta0,uchi0,vchi0
*     write (9,*) 'ueta,veta,uchi,vchi',ueta,veta,uchi,vchi
*     write (9,*) 'alpha1,beta1,condit1',alpha1,beta1,condit1
*     write (9,*) 'alpha2,beta2,condit2',alpha2,beta2,condit2
*     write (9,*) 'alpha,beta',alpha,beta
*     write (9,*) 'factor',(factor(i,j,lm),lm=1,3)
*     endif

               if( i .eq. 1 .and.
     $              ( j .gt. 1 .and. j .lt. nk2p2 ))then
                  ix = 2
               elseif ( i .eq. nj .and.
     $                 ( j .gt. 1 .and. j .lt. nk2p2 )) then
                  ix = 3
               else
                  ix = 1
               endif
               do 600 l=1,ngpts
               do 600 m=1,ngpts
C     
                  u0 = fint(x2ion,l,m,i,j)
                  v0 = fint(y2ion,l,m,i,j)
                  ueta0 = feta(x2ion,l,m,i,j)
                  uchi0 = fchi(x2ion,l,m,i,j)
                  veta0 = feta(y2ion,l,m,i,j)
                  vchi0 = fchi(y2ion,l,m,i,j)
C     
                  u = alpha*u0 + beta*v0
                  v = -beta*u0 + alpha*v0
                  ueta = (alpha*ueta0 + beta*uchi0)
                  uchi = (alpha*uchi0 - beta*ueta0)
                  veta = (alpha*veta0 + beta*vchi0)
                  vchi = (alpha*vchi0 - beta*veta0)
C     
                  uvjac = ueta*vchi-uchi*veta
                  geo1 = e(l,m,1,ix)*gxx(i,j) + e(l,m,2,ix)*gxx(i+1,j) +
     $                 e(l,m,3,ix)*gxx(i,j+1) + e(l,m,4,ix)*gxx(i+1,j+1)
                  geo2 = e(l,m,1,ix)*gyy(i,j) + e(l,m,2,ix)*gyy(i+1,j) +
     $                 e(l,m,3,ix)*gyy(i,j+1) + e(l,m,4,ix)*gyy(i+1,j+1)
                  geo3 = e(l,m,1,ix)*gxy(i,j) + e(l,m,2,ix)*gxy(i+1,j) +
     $                 e(l,m,3,ix)*gxy(i,j+1) + e(l,m,4,ix)*gxy(i+1,j+1)
                  geoh = e(l,m,1,ix)*gh(i,j) + e(l,m,2,ix)*gh(i+1,j) +
     $                 e(l,m,3,ix)*gh(i,j+1) + e(l,m,4,ix)*gh(i+1,j+1)
                  geoj = e(l,m,1,ix)*gj(i,j) + e(l,m,2,ix)*gj(i+1,j) +
     $                 e(l,m,3,ix)*gj(i,j+1) + e(l,m,4,ix)*gj(i+1,j+1)
C     
C     
                  ceta2 = ( geo1*vchi**2 + geo2*uchi**2
     $                 - 2.*uchi*vchi*geo3 )/ uvjac
                  cchi2 = ( geo1*veta**2 + geo2*ueta**2
     $                 - 2.*ueta*veta*geo3 )/ uvjac
                  cetachi = ( -geo1*veta*vchi -geo2*ueta*uchi
     $                 + uvjac*geo3 )/ uvjac
                  cj = uvjac*geoj
                  ch = geoh
*     if ( i .ge. 6 .and. i .le. 7 .and. j .ge. 9 .and. j .le. 10 )
*     $         then
*     write (9,*) 'u,v,ueta,veta,uchi,vchi',u,v,ueta,veta,uchi,vchi
*     write (9,*) 'uvjac,theta,sinth,sinover',uvjac,thet,sinth,sinover
*     write (9,*) 'geo1,geo2,geo3',geo1,geo2,geo3
*     write (9,*) 'ceta2,cchi2,cetachi,cj',ceta2,cchi2,cetachi,cj
*     write (9,*) 'ix  ',ix
*     endif
C     
C     loops are over
C     ii - elements
C     jj - function values
C     kl - conductivity values
C     
             do 500 ii=1,4
             do 500 jj=1,4
                cjwt(ii,jj) = cj*wtij(l,m)*e(l,m,ii,ix)*e(l,m,jj,ix)
     $               + cjwt(ii,jj)
             do 500 kl=1,4
                cpwt(ii,jj,kl)= cpwt(ii,jj,kl)
     $               + wtij(l,m)*e(l,m,kl,ix)*(
     $               eeta(l,m,ii,ix)*eeta(l,m,jj,ix)*ceta2  +
     $               echi(l,m,ii,ix)*echi(l,m,jj,ix)*cchi2  +
     $               (eeta(l,m,ii,ix)*echi(l,m,jj,ix)+
     $               echi(l,m,ii,ix)*eeta(l,m,jj,ix))*cetachi )
                chwt(ii,jj,kl) 
     $               = chwt(ii,jj,kl) + wtij(l,m)*e(l,m,kl,ix)*ch*
     $               (echi(l,m,ii,ix)*eeta(l,m,jj,ix)-
     $               eeta(l,m,ii,ix)*echi(l,m,jj,ix))
 500         continue
C     
 600      continue
C     
C     now put the coefficients together
C     
          do 800 ii=1,4         ! sum over elements
          ielem = mod(ii-1,2)
          jelem = (ii-1)/2
          do 800 jj=1,4      ! sum over function values
             aj(i+ielem,j+jelem,kstencl(ii,jj)) =
     $            aj(i+ielem,j+jelem,kstencl(ii,jj)) + cjwt(ii,jj)
          do 800 kl=1,4      ! sum over conductivity values
             apeder(i+ielem,j+jelem,kstencl(ii,jj),kstencl(ii,kl)) =
     $        apeder(i+ielem,j+jelem,kstencl(ii,jj),kstencl(ii,kl)) -
     $            cpwt(ii,jj,kl)
             ahall(i+ielem,j+jelem,kstencl(ii,jj),kstencl(ii,kl))  =
     $         ahall(i+ielem,j+jelem,kstencl(ii,jj),kstencl(ii,kl)) -
     $               chwt(ii,jj,kl)
 800      continue
C     
 1000  continue
C     
       do 1100 i=1,njp1
       do 1100 j=1,nk2p3
       do 1100 kl=1,9
       apot1(i,j,kl) = 0.
       do 1100 ll=1,9
          apot1(i,j,kl) = apot1(i,j,kl) + apeder(i,j,kl,ll)
 1100  continue
*     do 1110 i=1,njp1
*     do 1110 j=1,nk2p3
*     write (9,*) '*******************  ',i,j,'  **************'
*     sumit = 0.
*     do 1111 kl=1,9
*     sumitll(kl) = 0.
*     do 1111 ll = 1,9
*     sumitll(kl) = sumitll(kl) + apeder(i,j,kl,ll)
*     sumit = sumit + apeder(i,j,kl,ll)
*     1111  continue
*     write (9,*) 'total sum    ',sumit
*     write (9,*) 'sub sums   ',sumitll
*     write (9,*) 'summed stencil'
*     write (9,*) apot1(i,j,7),apot1(i,j,6),apot1(i,j,5)
*     write (9,*) apot1(i,j,8),apot1(i,j,9),apot1(i,j,4)
*     write (9,*) apot1(i,j,1),apot1(i,j,2),apot1(i,j,3)
*     write (9,*) 'total subarray'
*     write (9,1112) ((apeder(i,J,kl,ll),kl=1,9),ll=1,9)
*     1112    format(9(3x,1pe12.4))
*     1110  continue

       do 1200 i=1,njp1
       do 1200 j=1,nk2p3
          dom(i,j) = abs(apot1(i,j,9))
          do 1199 kl=1,8
             dom(i,j) = dom(i,j)-abs(apot1(i,j,kl))
 1199     continue
          dom(i,j) = dom(i,j)/apot1(i,j,9)
 1200  continue
!
! ... End SUBROUTINE COEF1..............................................
!
       RETURN
       END

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!> CORD()
      SUBROUTINE CORD(XA1,YA1,ZA1,xahf,yahf,zahf,XA2,YA2,ZA2)
!
! ... Local variables ..................................................
!
#include "ion-param.inc"
      REAL LFM_PLANET_RION,LFM_PLANET_FAC
#include "ion_P++_pointers.inc"
#include "ion_fortran_P++_arrays.inc"
#include "iongeo.inc"
#include "curint.inc"
      REAL S2THET(NJP1,NK2P3),ROUTR(NJP1,NK2P3)
      REAL ShfTHET(NJ,NK2),ROUThf(NJ,NK2)
!
! ... Parameter variables ..............................................
!
      REAL XA1(NJP1),YA1(NJP1,NK2P1),ZA1(NJP1,NK2P1)
!> Description of XA1,YA1,ZA1
      REAL XA2(NJP1,NK2P3),YA2(NJP1,NK2P3),ZA2(NJP1,NK2P3)
!> Description of XA2,YA2,ZA2
      REAL XAhf(nj,nk2),YAhf(nj,nk2),ZAhf(nj,nk2)
!> Description of xahf,yahf,zahf
!
! ... Begin ............................................................
!
#ifdef DEBUG_MODE_ON
      write(*,*) "DEBUG:  in ionosphere.F::CORD(...)"
#endif      
      R0 = LFM_PLANET_RION()
      RINR = XA1(1)**2+YA1(1,1)**2+ZA1(1,1)**2
      RINR = SQRT(RINR)
      DO 2 J=2,nk2p2
      DO 2 I=1,Njp1
      ROUTR(I,J) = XA2(I,j)**2+YA2(I,J)**2+ZA2(I,J)**2
    2 ROUTR(I,J) = SQRT(ROUTR(I,J))
      DO 4 J=1,NK2P1
      DO 4 I=1,NJP1
C Opt-in geometry repair: preserve the actual source wall. The dipole
C invariant uses each node's radius, with horizontal components evaluated
C directly to avoid subtracting nearly equal numbers at magnetic poles.
      RWALL=SQRT(XA1(I)**2+YA1(I,J)**2+ZA1(I,J)**2)
      S2THET(I,J)=(R0/RWALL)*((XA1(I)/RWALL)**2+
     $ (YA1(I,J)/RWALL)**2)
      TETA1(I,J) = ASIN(SQRT(S2THET(I,J)))
    4 PHI1(I,J) = ATAN2(YA1(I,J),XA1(I))
        do 5 j=1,nk2p1
        do 5 i=1,njp1
         x1ion(i,j) = teta1(i,j)*cos(phi1(i,j))
         y1ion(i,j) = teta1(i,j)*sin(phi1(i,j))
    5   continue
*
*  standard mapping of grid
*
      DO 12 J=2,Nk2p2
      DO 12 I=1,Njp1
      S2THET(I,J) = R0/ROUTR(I,J)*(1.-(ZA2(I,J)/
     1              ROUTR(I,J))**2)
       TETA2(I,J) = ASIN(SQRT(S2THET(I,J)))
   12  PHI2(I,J) = ATAN2(YA2(I,J),XA2(I,j))
        do 15 j=2,nk2p2
        do 15 i=1,njp1
        x2ion(i,j) = teta2(i,j)*cos(phi2(i,j))
        y2ion(i,j) = teta2(i,j)*sin(phi2(i,j))
   15   continue
        do 16 i=1,njp1
           x2ion(i,1) = 2*x2ion(i,2) - x2ion(i,3)
           y2ion(i,1) = 2*y2ion(i,2) - y2ion(i,3)
           x2ion(i,nk2p3) = 2*x2ion(i,nk2p2) - x2ion(i,nk2p1)
           y2ion(i,nk2p3) = 2*y2ion(i,nk2p2) - y2ion(i,nk2p1)
 16     continue
        do j=1,nk2p3
           do i =1,njp1
              y2ion_old(i,j) = y2ion(i,j)
           enddo
        enddo
#ifdef Y2INTERP
*
*  constant line across r = x20*sec(phi2), even spacing in delta y
*
       do i=1,njp1
          deltay = (y2ion(i,nk2p3)-y2ion(i,1))/float(nk2p2)
          y2_0 = y2ion(i,1)
          do j=1,nk2p3
             y2ion(i,j) = y2_0 +deltay*float(j-1)
          enddo
       enddo
*             
* end of straight line mapping
*
c
c  get interpolation coefficients
c
        call cur_coefs
c
c
#endif
c       
      DO 20 J=1,nk2
      DO 20 I=1,Nj
      ROUTHF(I,J) = XAhf(I,j)**2+YAhf(I,J)**2+ZAhf(I,J)**2
   20 ROUTHF(I,J) = SQRT(ROUTHF(I,J))
      DO 22 J=1,Nk2
      DO 22 I=1,Nj
      ShfTHET(I,J) = R0/ROUTHF(I,J)*(1.-(ZAhf(I,J)/
     1              ROUTHF(I,J))**2)
       TETAhf(I,J) = ASIN(SQRT(ShfTHET(I,J)))
   22  PHIhf(I,J) = ATAN2(YAhf(I,J),XAhf(I,j))
        do 25 j=1,nk2
        do 25 i=1,nj
        xhfion(i,j) = tetahf(i,j)*cos(phihf(i,j))
        yhfion(i,j) = tetahf(i,j)*sin(phihf(i,j))
   25   continue
!
! ... End SUBROUTINE CORD(XA1,YA1,ZA1,xahf,yahf,zahf,XA2,YA2,ZA2).......
!
      RETURN
      END

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!> pardy - Compute the Hardy addition
!!    This routine computes the Hardy addition to the ionospheric
!!    conductivity due to the diffuse aurora.
!!
!!    The two parameters are:
!!    b0 = 1 for discrete; 0 otherwise
!!    c0 = scale factor for everything
!!    h0 = dummy
      subroutine pardy(dsigmap,dsigmah,dflux,denergy,b0,c0,h0)
!
! ... Local variables ..................................................
!
#include "ion-param.inc"
      REAL LFM_PLANET_RION,LFM_PLANET_FAC
#include "ion_P++_pointers.inc"
#include "ion_fortran_P++_arrays.inc"
#include "iongeo.inc"
      real discrete(njp1,nk2p3),fct1(njp1,nk2p3)
      real fct2(njp1,nk2p3)      
      common /sigpeuv/sigpeuv(njp1,nk2p3)
      common /ionrho/rhoion(njp1,nk2p3),cion(njp1,nk2p3),
     $     bjion(njp1,nk2p3)
      data ntimes/0/
      save ntimes,fct1,fct2
!
! ... Parameter variables ..............................................
!
      real dsigmap(njp1,nk2p3)
!> dsigmap description
      real dsigmah(njp1,nk2p3)
!> dsigmah description
      real dflux(njp1,nk2p3)
!> dflux description
      real denergy(njp1,nk2p3)
!> denergy description
!
! ... Begin ............................................................
!
#ifdef DEBUG_MODE_ON
      write(*,*) "DEBUG:  in ionosphere.F::pardy(...)"
#endif      
!      fpoop(w) = (1./(1.+(w)**4))
      ntimes = ntimes + 1
      if( ntimes .eq. 1 ) then
C damp
      redge = sqrt( x2ion(nj/2,1)**2 + y2ion(nj/2,1)**2)*57.3
      rmin = 0.85*redge
      rmin2 = 0.45*redge  ! this needs work for disturbed conditions
      rmax = 1.0*redge
      rfac = 0.1
      rmn = rmin/57.3
      rmn2 = rmin2/57.3
      rmx = rmax/57.3
      do k=1,nk2p3
         do j=1,njp1
          rnow = sqrt( x2ion(j,k)**2 + y2ion(j,k)**2)
          snow = rnow
          rnow = amin1( amax1(rnow,rmn),rmx )
          snow2 = amin1( amax1(snow,rmn2),rmx)
          fct = 1. + (rfac-1.)*(rnow-rmn)/(rmx-rmn)
          fct1(j,k) = fct
          fctt = (snow2-rmn2)/(rmx-rmn2)
          fct2(j,k) = 1. - 0.33*fctt*y2ion(j,k)/(snow+1.e-10)
         enddo
      enddo
         else
       endif
C 
C zero the parallel drop for NKILL grid points
C
      NKILL = 4
      do j=1,njp1
       do k=1,nk2p3
C dsigmah is E0
C fct2 gives a twist to the energy to account for grad B electron drift
C    at low latitude
c         dsigmah(j,k) = 12.52e-16*cion(j,k)**2*fct2(j,k)
*         dsigmah(j,k) = h0*12.52e-16*cion(j,k)**2*fct2(j,k)
c 4/18/95 put in new constant to correct to final Viking Model (L) h0=1
          dsigmah(j,k) = h0*9.39e-16*cion(j,k)**2*fct2(j,k)
c
C dsigmap is PHI0
*         dsigmap(j,k) = 3.26e20*rhoion(j,k)*sqrt(dsigmah(j,k)*1000.)
          dsigmap(j,k) = 0.815e20*rhoion(j,k)*sqrt(dsigmah(j,k)*1000.)
C bjion is delta phi
C enhance the discrete pot. by 8 or 1.6 (anom. resist.)
C reduce this by 1/2 for Galileo
         bsgn = sign(1.,bjion(j,k))
*         ares = (4.*(1.+bsgn) + 0.8*(1.-bsgn))*fct1(j,k)*0.5
c 4/18/95 put in new constant to corres to final Viking Model (L)
         ares = (6.*(1.+bsgn) + 1.2*(1.-bsgn))*fct1(j,k)
         fkill = float(1-int(j/(NJ-NKILL)))
C 4/18/95 put in new constant to corres to final VIKING model (L) (b0=1)
         rhomin = 3.3e-24*0.5*sigpeuv(j,k)
         bjion(j,k) = b0*bjion(j,k)/amax1(rhoion(j,k),rhomin)
     @      *15.73e-19*sqrt(dsigmah(j,k))*ares*fkill
C     @      *7.865e-19*sqrt(dsigmah(j,k))*ares*fkill
         bjion(j,k) = amin1(bjion(j,k),20.)
       enddo
      enddo
      do j=1,njp1
       do k=1,nk2p3
C dsigmap is corrected phi
         if(bjion(j,k).gt.0.) then
             crct = (8. -7.*exp(-(bjion(j,k)/(7.*dsigmah(j,k)))))
           else
             crct = exp(bjion(j,k)/dsigmah(j,k))
         endif
         dsigmap(j,k) = crct*dsigmap(j,k)
	 dflux(j,k) = dsigmap(j,k)
C dsigmah is corrected E0
         dsigmah(j,k) = dsigmah(j,k) + bjion(j,k)
         dsigmah(j,k) = amax1( dsigmah(j,k), 1.e-8)
	 denergy(j,k) = dsigmah(j,k)
       enddo
      enddo
C
#ifdef DEBUGHARDY
         nch = 1
C         call bilinion
         call iondiah('RHO and C',rhoion,cion,lstep,time,nch)
C temp
C      write(47,*) ' PHI corr '
C      write(47,847) NJP1,NK2P3
C      write(47,848) (( dsigmap(j,k),j=1,njp1 ), k=1,nk2p3 )
C      write(47,*) ' E0 corr '
C      write(47,847) NJP1,NK2P3
C      write(47,848) (( dsigmah(j,k),j=1,njp1 ), k=1,nk2p3 )
C      write(47,*) ' delta phi '
C      write(47,847) NJP1,NK2P3
C      write(47,848) (( bjion(j,k),j=1,njp1 ), k=1,nk2p3 )
#endif

c
      do k=1,nk2p3
         do j=1,njp1
C Hardy and discrete enhancement to the Pedersen cond.
c            dsigmap(j,k) = c0*2.5*dsigmah(j,k)**1.5*sqrt(dsigmap(j,k))
C 4/18/95 in the Viking parameters c0 was 2; now use c0 = 1 by mult by 2.5
            dsigmap(j,k) = c0*5.0*dsigmah(j,k)**1.5*sqrt(dsigmap(j,k))
     $                   /(1.+(0.25*dsigmah(j,k))**2)*fct1(j,k)
C put the cutoff ramp here
         enddo
      enddo

      do k=1,nk2p3
         do j=1,njp1
          h00 = dsigmah(j,k)**0.85/(1.+(dsigmah(j,k)*0.05)**2)
            dsigmah(j,k) =  -0.45*h00*dsigmap(j,k)
         enddo
      enddo
#ifdef DEBUGHARDY
c     ncar
          nch = 2
        call iondiah('Hardy correction',dsigmap,dsigmah,lstep,time,nch)
      write(9,*) ' PED/HAR '
      write(9,89) NJP1,NK2P3
89    format(x,2i5)
848   format(x,1p11e11.3)
      write(9,848) (( dsigmap(j,k),j=1,njp1 ), k=1,nk2p3 )
      write(9,*) ' HAL/HAR '
      write(9,89) NJP1,NK2P3
      write(9,848) (( dsigmah(j,k),j=1,njp1 ), k=1,nk2p3 )
#endif
!
! ... End subroutine pardy(dsigmap,dsigmah,dflux,denergy,b0,c0,h0)......
!
      return
      end
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!> intrp - Interpolate the potential onto the interior surface mesh
!! (long description)
      subroutine intrp
!
! ... Local variables ..................................................
!
#include "ion-param.inc"
      REAL LFM_PLANET_RION,LFM_PLANET_FAC
#include "volts.inc"
#include "iongeo.inc"
!
! ... Parameter variables ..............................................
!
! ... Begin ............................................................
!
#ifdef DEBUG_MODE_ON
      write(*,*) "DEBUG:  in ionosphere.F::intrp(...)"
#endif
        do 10 k=1,nk2p1
        do 10 j=1,njp1
	if ( eta1(j,k) .lt. 1.0 .and. chi1(j,k) .lt. 1.0) then
	 psii(j,k) = 0.0
	else
         jj = chi1(j,k)
         ii = eta1(j,k)
         eta11 = eta1(j,k)-float(ii)
         chi11 = chi1(j,k)-float(jj)
         psii(j,k) = psi(ii,jj)*(1.-eta11)*(1.-chi11) +
     1          psi(ii+1,jj)*eta11*(1.-chi11) +
     2          psi(ii,jj+1)*(1.-eta11)*chi11 +
     3          psi(ii+1,jj+1)*eta11*chi11
	endif
   10   continue
        do 20 j=1,nk2p1
        psii(1,j) = 0.
        psii(njp1,j) = 0.
  20    continue
        do 21 i=2,nj
        psii(i,1) = 0.
        psii(i,nk2p1) = 0.
  21    continue
C
        do 50 k=1,nk2
        do 50 j=1,nj
           if ( etahf(j,k) .lt. 1.0 .and. chihf(j,k) .lt. 1.0) then
              psieta(j,k) = 0.0
              psichi(j,k) = 0.0
           else
          jj=  chihf(j,k)
          ii = etahf(j,k)
          chihf1 = chihf(j,k)-float(jj)
          etahf1 = etahf(j,k) - float(ii)
            psieta(j,k) = (psi(ii+1,jj)-psi(ii,jj))*(1.-chihf1)+
     $               (psi(ii+1,jj+1)-psi(ii,jj+1))*chihf1
            psichi(j,k) = (psi(ii,jj+1)-psi(ii,jj))*(1.-etahf1)+
     $           (psi(ii+1,jj+1)-psi(ii+1,jj))*etahf1
            endif
   50   continue
C
        do 60 k=1,nk2
        do 60 j=1,nj
          vxint(j,k) = psieta(j,k)*vxeta(j,k) + psichi(j,k)*vxchi(j,k)
          vyint(j,k) = psieta(j,k)*vyeta(j,k) + psichi(j,k)*vychi(j,k)
          vzint(j,k) = psieta(j,k)*vzeta(j,k) + psichi(j,k)*vzchi(j,k)
   60   continue
!
! ... End subroutine intrp..............................................
!
      return
      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> iondiah
!! Ionosphere diagnostics????
	subroutine iondiah(title,one,two,lstep,time,nch)
!
! ... Local variables ..................................................
!
        real*8 time
#include "ion-param.inc"
      REAL LFM_PLANET_RION,LFM_PLANET_FAC
	character*80 titl
	dimension onea(njp3,nk2p3),twoa(njp3,nk2p3)	
!
! ... Parameter variables ..............................................
!
	character*(*) title
!> title description
	dimension one(njp1,nk2p3)
!> one description
        dimension two(njp1,nk2p3)	
!> two description
!
! ... Begin ............................................................
!
#ifdef DEBUG_MODE_ON
      write(*,*) "DEBUG:  in ionosphere.F::iondiah(...)"
#endif
c
	titl = title
	do k=1,nk2p3
	do j=1,njp1
	  onea(j,k) = one(j,k)
	  twoa(j,k) = two(j,k)
	enddo
	enddo
c
#ifdef DIAGPLOTS
C	call bilinion
C	call iondiag(titl,onea,twoa,lstep,time)
#endif
c
!
! ... End subroutine iondiah(title,one,two,lstep,time,nch).............
!
	return
	end
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!> ioninit
!!
      subroutine ioninit
!
! ... Local variables ..................................................
!
#include "ion-param.inc"
      REAL LFM_PLANET_RION,LFM_PLANET_FAC
#include "ion-var.inc"
#include "ion_P++_pointers.inc"
#include "ion_fortran_P++_arrays.inc"
#include "global_dims.inc"
c
      common /ionscrach/
     $     xa1(njp1),ya1(njp1,nk2p1),za1(njp1,nk2p1),
     $     xa2(njp1,nk2p3),ya2(njp1,nk2p3),za2(njp1,nk2p3),
     $     xahf(nj,nk2),yahf(nj,nk2),zahf(nj,nk2)
!
! ... Parameter variables ..............................................
!
! ... Begin ............................................................
!
#ifdef DEBUG_MODE_ON
      write(*,*) "DEBUG:  in ionosphere.F::ioninit(...)"
#endif
c
      mion_i = nion_i
      mion_j = nion_j
C     note the 1.e-2 to change from centimeters to meters for the
C     ionospheric part of the code
C
      do 200 j=1,njp1
          xa1(j) = x(1,j,1)*1.e-2
C
          do 190 k=1,nk2P1
          ya1(j,k) = y(1,j,k)*1.e-2
          za1(j,k) = z(1,j,k)*1.e-2
          xa2(j,k+1) = x(2,j,1)*1.e-2
          ya2(j,k+1) = y(2,j,k)*1.e-2
          za2(j,k+1) = z(2,j,k)*1.e-2
  190     continue
C
C
  200   continue
        do 201 k=1,nk2p1
        xa2(1,k+1) = x(2,1,nk2/2+1)*1.e-2
        ya2(1,k+1) = y(2,1,nk2/2+1)*1.e-2
        za2(1,k+1) = z(2,1,nk2/2+1)*1.e-2
        xa2(njp1,k+1) = x(2,njp1,nk2/2+1)*1.e-2
        ya2(njp1,k+1) = y(2,njp1,nk2/2+1)*1.e-2
        za2(njp1,k+1) = z(2,njp1,nk2/2+1)*1.e-2
  201   continue
C
*        do 240 j=1,njp1
*        xa2(j,1) = x(1,j,1)*1.e-2
*        xa2(j,nk2p3) = x(1,j,nk2p1)*1.e-2
*        ya2(j,1) = y(1,j,1)*1.e-2
*        ya2(j,nk2p3) = y(1,j,nk2p1)*1.e-2
*        za2(j,1) = z(1,j,1)*1.e-2
*        za2(j,nk2p3) = z(1,j,nk2p1)*1.e-2
*  240   continue
C
C
        do 300 k=1,nk2
        do 300 j=1,nj
        xahf(j,k) = 1.25e-3*( x(1,j,k) + x(2,j,k)
     $          + x(1,j+1,k) + x(2,j+1,k)
     $          + x(1,j+1,k+1) + x(2,j+1,k+1)
     $          + x(1,j,k+1) + x(2,j,k+1))
        yahf(j,k) = 1.25e-3*( y(1,j,k) + y(2,j,k)
     $          + y(1,j+1,k) + y(2,j+1,k)
     $          + y(1,j+1,k+1) + y(2,j+1,k+1)
     $          + y(1,j,k+1) + y(2,j,k+1))
        zahf(j,k) = 1.25e-3*( z(1,j,k) + z(2,j,k)
     $          + z(1,j+1,k) + z(2,j+1,k)
     $          + z(1,j+1,k+1) + z(2,j+1,k+1)
     $          + z(1,j,k+1) + z(2,j,k+1))
  300   continue
C
C
        call cord(xa1,ya1,za1,xahf,yahf,zahf,xa2,ya2,za2)
C
        call subwts
        call vcoef(xahf,yahf,zahf)
C
C  note coef1 uses /ionscrach/ so we need to be done with it before
C       it is called
C
        call coef1
        call jsetup
!
! ... End subroutine ioninit ...........................................
!
        return
        end
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!> ionosf()
!! 
      subroutine ionosf
!
! ... Local variables ..................................................
!
#include "ion-param.inc"
      REAL LFM_PLANET_RION,LFM_PLANET_FAC
#include "ion_P++_pointers.inc"
#include "ion_fortran_P++_arrays.inc"
#include "jpara.inc"
#include "ion-meter.inc"
#include "volts.inc"
#include "ion-var.inc"
#include "ionC90.inc"
#include "run-time.inc"
      real*4 time4
#include "boundx.inc"
***#include "iongeo.inc"
#include "curint.inc"
      common /ionrho/rhoion(njp1,nk2p3),cion(njp1,nk2p3),
     $     bjion(njp1,nk2p3)
      dimension cdum(njp1,nk2p3),rhodum(njp1,nk2p3),
     $     ciondum(njp1,nk2p3)
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>  JGL
      integer counter
      data counter /-1/

      dimension vxint_n(nj,nk2),vyint_n(nj,nk2),vzint_n(nj,nk2)
      dimension vxint_s(nj,nk2),vyint_s(nj,nk2),vzint_s(nj,nk2)
      dimension psii_n(njp1,nk2p1),psii_s(njp1,nk2p1)
      integer modulus
      save counter
      save vxint_n,vyint_n,vzint_n
      save vxint_s,vyint_s,vzint_s
      save psii_n,psii_s
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
c
#ifdef SMOOTHION
      dimension curd(njp1,nkp1)
#endif
*      dimension psiv(2,njp1,nkp1)
#ifdef DEBUG
      dimension taudum(njp1,nk),curdum(NJP1,NK),CURDUM2(NJP1,nk)
C /* #define PSICOMPARE 1 Needs conraq_ */
!
! ... Parameter variables...............................................
!
! ... Begin ............................................................
!
#ifdef DEBUG_MODE_ON
      write(*,*) "DEBUG:  in ionosphere.F::ionosf(...)"
#endif
c
      write (6,*) 'in ionosf, modulus =', modulus
#endif /*DEBUG*/
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>   JGL
c
      time4=real(time,4)
       if (table) then
        ind = min(itable,nint((time4-tzero)/deltaT))
	if(ind.lt.1) ind = 1
                 if (SM) tilt_angle = Zangle(ind)
C                write(6,*) ind,'tilt_angle changed to ',tilt_angle
c       
           endif
      counter = counter + 1
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
C
#ifdef DEBUG
      write (6,*) counter, 'in full branch'
      call pictur(bx,'BX',1,2,1,njp1,1,nkp1,0,2,3,1,
     $     nip1,njp1,nkp1,nip1*njp1)
      call pictur(by,'BY',1,2,1,njp1,1,nkp1,0,2,3,1,
     $     nip1,njp1,nkp1,njp1*nip1)
      call pictur(bz,'BZ',1,2,1,njp1,1,nkp1,0,2,3,1,
     $     nip1,njp1,nkp1,njp1*nip1)
#endif /*DEBUG*/
      call jpara(bx,by,bz)
#ifdef DEBUG
      write(6,*) ' IONOSP: current(1,1)=',current(1,1)
      write(9,*) ' IONOSP: current(1,1)=',current(1,1)
      call print2d(current(1,1),'jpar',njp1,nk,njp1,nk)
#endif /*DEBUG*/
      dtconv = dt*1.e8
#ifdef DEBUG
      do 1100 k=1,nk
      do 1100 j=1,njp1
 1100    curdum(j,k) = current(j,k)
#endif /*DEBUG*/
c
c
#ifdef SMOOTHION
         do k=2,nk-1
            do j=2,nj
               curd(j,k) = 0.25*current(j,k) + 0.125*(
     $              current(j+1,k)+current(j,k+1)+current(j-1,k)+
     $              current(j,k-1)) + 0.0625*(current(j-1,k-1)+
     $              current(j-1,k+1)+current(j+1,k-1)+
     $              current(j+1,k+1))
            enddo
         enddo
         do k=2,nk-1
            do j=2,nj
               current(j,k) = curd(j,k)
            enddo
         enddo
#endif
         do k=1,nkp1
            DO J=1,NJP1
               rhoion0(j,k) = 0.
               cion0(j,k) = 0.
            enddo
         enddo
         do k=2,nk
            klow = mod(k-2+nk,nk) + 1
            khi = mod(k-1+nk,nk) + 1
            do j=2,nj
            rhoion0(j,k) = 0.125*(rho(2,j,klow)+rho(2,j-1,khi)+
     $           rho(2,j,klow)+rho(2,j-1,klow)+
     $           rho(1,j,khi)+rho(1,j-1,khi)+
     $           rho(1,j,klow)+rho(1,j-1,klow))
            cion0(j,k) = 0.125*(c(2,j,khi)+c(2,j-1,khi)+
     $           c(2,j,klow)+c(2,j-1,klow)+
     $           c(1,j,khi)+c(1,j-1,khi)+
     $           c(1,j,klow)+c(1,j-1,klow))
            enddo
         enddo
c               
C     contains a conversion from volts to statvolts and also
C     a c for Faraday's law
C
C   do the northern hemisphere
C
      do 200 k=2,nk2
         do 200 j=2,nj
            curnt(j,k+1) = -current(j,k)
            rhoion(j,k+1) = rhoion0(j,k)
            cion(j,k+1) = cion0(j,k)
 200  continue
#ifdef Y2INTERP
      do k=2,nk2p1
         do j=2,nj
            loff = lmin(j,k)-1
            cdum(j,k) = 0.
            rhodum(j,k) = 0.
            ciondum(j,k) = 0.
            do lint = lmin(j,k),lmax(j,k)
               cdum(j,k) = cdum(j,k) + curcoef(j,k,lint-loff)*
     $              curnt(j,lint)
               rhodum(j,k) = rhodum(j,k) +
     $              curcoef(j,k,lint-loff)*rhoion(j,lint)
               ciondum(j,k) = ciondum(j,k) +
     $              curcoef(j,k,lint-loff)*cion(j,lint)
            enddo
         enddo
      enddo
c
#ifdef IONINTERP_CHECK
      write (21,*) njp1,nk2p3
      write (21,901) ((y2ion(j,k),j=1,njp1),k=1,nk2p3)
      write (21,901) ((y2ion_old(j,k),j=1,njp1),k=1,nk2p3)
      write (21,901) ((curnt(j,k),j=1,njp1),k=1,nk2p3)
      write (21,901) ((cdum(j,k),j=1,njp1),k=1,nk2p3)
      write (21,901) ((rhoion(j,k),j=1,njp1),k=1,nk2p3)
      write (21,901) ((rhodum(j,k),j=1,njp1),k=1,nk2p3)
      write (21,901) ((cion(j,k),j=1,njp1),k=1,nk2p3)
      write (21,901) ((ciondum(j,k),j=1,njp1),k=1,nk2p3)
 901  format(10(1pe12.4))
#endif
      do k=2,nk2p1
         do j=2,nj
            curnt(j,k) = cdum(j,k)
            rhoion(j,k) = rhodum(j,k)
            cion(j,k) = ciondum(j,k)
         enddo
      enddo
#endif               
C
      do 220 j=1,njp1
         curnt(j,1) = 0.
         curnt(j,2) = 0.
         curnt(j,NK2+2) = 0.
         curnt(j,nk2p3) = 0.
         rhoion(j,1) = rhoion(j,3)
         rhoion(j,2) = rhoion(j,3)
         rhoion(j,nk2+2) = rhoion(j,nk2p1)
         rhoion(j,nk2p3) = rhoion(j,nk2p1)
         cion(j,1) = cion(j,3)
         cion(j,2) = cion(j,3)
         cion(j,nk2+2) = cion(j,nk2+1)
         cion(j,nk2p3) = cion(j,nk2+1)
 220  continue
      do 221 k=1,nk2p3
         curnt(1,k) = 0.
         curnt(njp1,k) = 0.
         rhoion(1,k) = rhoion(2,k)
         rhoion(njp1,k) = rhoion(nj,k)
         cion(1,k) = cion(2,k)
         cion(njp1,k) = cion(nj,k)
 221  continue
c
      do k=1,nk2p3
         do j=1,njp1
            bjion(j,k) = curnt(j,k)
         enddo
      enddo
C
C      call print2d(rhoion,'RHO N',NJP1,NK2p3,NJp1,NK2P3)
C      call print2d(cion,'C N',NJP1,NK2p3,NJp1,NK2P3)
C
c$$$      do 250 k=4,nk2
c$$$         do 250 j=3,njm1
c$$$            psi(j,k) = psinorth(j,k)
c$$$            curnorth(j,k) = curnt(j,k)
c$$$ 250  continue
      do 250 k=1,nk2p3
         do 250 j=1,njp1
            psi(j,k) = psinorth(j,k)
            curnorth(j,k) = curnt(j,k)
 250  continue

C
#ifdef DEBUG
      write(9,*) ' IONOSF: about to execute PSISLV for N'
#endif

      call psislv(err,.true.)

#ifdef DEBUG
      write(9,*) ' IONOSF: PSISLV for N err=',err
#endif
C
      psimax = psi(1,1)
      psimin = psi(1,1)
      cnabs = 0.
      cndiff = 0.
      do 251 k=1,nk2p3
         do 251 j=1,njp1
            sig_p_N(j,k) = sigmap(j,k)
            sig_h_N(j,k) = sigmah(j,k)
            cnabs = cnabs + abs(curnt(j,k))
            cndiff = cndiff + curnt(j,k)
            flux_N(j,k) = dflux(j,k)
            av_energy_N(j,k) = denergy(j,k)
            psimax = max(psimax,psi(j,k))
            psimin = min(psimin,psi(j,k))
 251  continue
C
*      write (17,*) 'IONOSPHERE-NORTH  ',psimax,psimin,cnabs,cndiff
      do 260 k=1,nk2p3
         do 260 j=1,njp1
            psinorth(j,k) = psi(j,k)
 260  continue
C
      call intrp
      CALL LFM_PLANET_MAP_COROTATION
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>  JGL
      do k=1,nk2
         do j=1,nj
            vxint_n(j,k) = vxint(j,k)
            vyint_n(j,k) = vyint(j,k)
            vzint_n(j,k) = vzint(j,k)
          enddo
       enddo
       do k=1,nk2p1
          do j=1,njp1
             psii_n(j,k) = psii(j,k)
          enddo
       enddo
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
C
      call psicomp(psi,psii)
#ifdef DEBUG
      call print2d(psii,'psii',njp1,nk2p1,njp1,nk2p1)
#endif
c
c  >>>>>>>>>>>>>>>>>>>
c      note the 1.e8 to convert from mks to cgs
c
      do k=1,nk2p1
         do j=1,nj
            ej_pot(1,j,k) = 1.e8*(psii(j+1,k) - psii(j,k))
         enddo
      enddo
      do k=1,nk2
         do j=1,njp1
            ek_pot(1,j,k) = 1.e8*(psii(j,k+1) - psii(j,k))
         enddo
      enddo
      do 320 k=1,nk2
         do 320 j=1,nj
            vx_pot(1,j,k) = vxint(j,k)
            vy_pot(1,j,k) = vyint(j,k)
            vz_pot(1,j,k) = vzint(j,k)
 320  continue
#ifdef DEBUG
      CALL print2d(vxint,'vxint',nj,nk2,nj,nk2)
      CALL print2d(vyint,'vyint',nj,nk2,nj,nk2)
      CALL print2d(vzint,'vzint',nj,nk2,nj,nk2)
#endif
C
C     now the southern hemisphere
C
C     note the inversion in the k index; done to make the two
C     hemispheres look indentical when viewed from above the
C     respective pole, also note the sign inversion done on the
C     potential
C
      do 500 k=2,nk2
         do 500 j=2,nj
            curnt(j,k+1) = current(j,nkp2-k)
            rhoion(j,k+1) = rhoion0(j,nkp2-k)
            cion(j,k) = cion0(j,nkp2-k)
 500  continue
C
#ifdef Y2INTERP
      do k=2,nk2p1
         do j=2,nj
            loff = lmin(j,k)-1
            cdum(j,k) = 0.
            rhodum(j,k) = 0.
            ciondum(j,k) = 0.
            do lint = lmin(j,k),lmax(j,k)
               cdum(j,k) = cdum(j,k) + curcoef(j,k,lint-loff)
     $              *curnt(j,lint)
               rhodum(j,k) = rhodum(j,k) +
     $              curcoef(j,k,lint-loff)*rhoion(j,lint)
               ciondum(j,k) = ciondum(j,k) +
     $              curcoef(j,k,lint-loff)*cion(j,lint)
            enddo
         enddo
      enddo
c
#ifdef IONINTERP_CHECK
      write (21,*) njp1,nk2p3
      write (21,901) ((y2ion(j,k),j=1,njp1),k=1,nk2p3)
      write (21,901) ((y2ion_old(j,k),j=1,njp1),k=1,nk2p3)
      write (21,901) ((curnt(j,k),j=1,njp1),k=1,nk2p3)
      write (21,901) ((cdum(j,k),j=1,njp1),k=1,nk2p3)
      write (21,901) ((rhoion(j,k),j=1,njp1),k=1,nk2p3)
      write (21,901) ((rhodum(j,k),j=1,njp1),k=1,nk2p3)
      write (21,901) ((cion(j,k),j=1,njp1),k=1,nk2p3)
      write (21,901) ((ciondum(j,k),j=1,njp1),k=1,nk2p3)
#endif
      do k=2,nk2p1
         do j=2,nj
            curnt(j,k) = cdum(j,k)
            rhoion(j,k) = rhodum(j,k)
            cion(j,k) = ciondum(j,k)
         enddo
      enddo
#endif
c
      do 520 j=1,njp1
         curnt(j,1) = 0.
         curnt(j,2) = 0.
         curnt(j,nk2+2) = 0.
         curnt(j,nk2p3) = 0.
         rhoion(j,1) = rhoion(j,3)
         rhoion(j,2) = rhoion(j,3)
         rhoion(j,nk2+2) = rhoion(j,nk2+1)
         rhoion(j,nk2p3) = rhoion(j,nk2+1)
         cion(j,1) = cion(j,3)
         cion(j,2) = cion(j,3)
         cion(j,nk2+2) = cion(j,nk2+1)
         cion(j,nk2p3) = cion(j,nk2+1)
 520  continue
      do 521 k=1,nk2p3
         curnt(1,k) = 0.
         curnt(njp1,k) = 0.
         rhoion(1,k) = rhoion(2,k)
         rhoion(njp1,k) = rhoion(nj,k)
         cion(1,k) = cion(2,k)
         cion(njp1,k) = cion(nj,k)
 521  continue
c
      do k=1,nk2p3
         do j=1,njp1
            bjion(j,k) = curnt(j,k)
         enddo
      enddo
C
C      call print2d(rhoion,'RHO S',NJP1,NK2p3,NJp1,NK2P3)
C      call print2d(cion,'C S',NJP1,NK2p3,NJp1,NK2P3)
C
c$$$      do 550 k=4,nk2
c$$$         do 550 j=3,njm1
c$$$            psi(j,k) = psisouth(j,k)
c$$$            cursouth(j,k) = curnt(j,k)
c$$$ 550  continue
      do 550 k=1,nk2p3
         do 550 j=1,njp1
            psi(j,k) = psisouth(j,k)
            cursouth(j,k) = curnt(j,k)
 550  continue
C
#ifdef DEBUG
      write(9,*) ' IONOSF: about to execute PSISLV for S'
#endif

      call psislv(err,.false.)

#ifdef DEBUG
      write(9,*) ' IONOSF: PSISLV for S err=',err
#endif
C 
      psimax = psi(1,1)
      psimin = psi(1,1)
      csabs = 0.
      csdiff = 0.
      do 562 k=1,nk2p3
         do 562 j=1,njp1
            sig_p_S(j,k) = sigmap(j,k)
            sig_h_S(j,k) = sigmah(j,k)
            flux_S(j,k) = dflux(j,k)
            av_energy_S(j,k) = denergy(j,k)
            csabs = csabs + abs(curnt(j,k))
            csdiff = csdiff + curnt(j,k)
            psimax = max(psimax,psi(j,k))
            psimin = min(psimin,psi(j,k))
 562  continue
*      write (18,*) 'IONOSPHERE-SOUTH  ',psimax,psimin,csabs,csdiff    
C

      do 570 k=1,nk2p3
         do 570 j=1,njp1
            psisouth(j,k) = psi(j,k)
 570  continue
C
*      write(9,*) 'north/south currents, abs(j) and net j  ',cnabs,
*     $     cndiff,csabs,csdiff
      call intrp
      CALL LFM_PLANET_MAP_COROTATION
      call psicomp(psi,psii)
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>  JGL
      do k=1,nk2
         do j=1,nj
            vxint_s(j,k) = vxint(j,k)
            vyint_s(j,k) = vyint(j,k)
            vzint_s(j,k) = vzint(j,k)
          enddo
       enddo
       do k=1,nk2p1
          do j=1,njp1
             psii_s(j,k) = psii(j,k)
          enddo
       enddo
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
C
C     note the k index inversion again
C
c  >>>>>>>>>>>>>>>>>>>
c      note the 1.e8 to convert from mks to cgs
c
      do 600 k=1,nk2p1
         do 600 j=1,nj
            ej_pot(1,j,nkp2-k) = 1.e8*(psii(j+1,k)-psii(j,k))
 600  continue
c
c  note the inversion here for the k index
c
      do 610 k=1,nk2
         do 610 j=1,njp1
            ek_pot(1,j,nkp1-k) = 1.e8*(psii(j,k) - psii(j,k+1))
 610  continue
      do 620 k=1,nk2
         do 620 j=1,nj
            vx_pot(1,j,nkp1-k) = vxint(j,k)
            vy_pot(1,j,nkp1-k) = vyint(j,k)
            vz_pot(1,j,nkp1-k) = -vzint(j,k)
 620  continue
C     note the signs here, only a half flip because of the interchange
C     
C     velocity boundary condition
!
! ... End subroutine ionosf.............................................
!
      return
      end
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!> jpara(bxx,byy,bzz)
!!
        Subroutine jpara(bxx,byy,bzz)
!
! ... Local variables ..................................................
!
#include "ion-param.inc"
      REAL LFM_PLANET_RION,LFM_PLANET_FAC
#include "ion_P++_pointers.inc"
#include "ion_fortran_P++_arrays.inc"
#include "jpara.inc"
#include "ion-meter.inc"
      SAVE
      common /ionscrach/
     $     xc(2,nj,nkp1),yc(2,nj,nkp1),zc(2,nj,nkp1),
     $     ji(2,njm1,nk),jj(1,nj,nk),jk(1,njm1,nkp1),
     $     jx(1,njp1,nkp1),jy(1,njp1,nkp1),jz(1,njp1,nkp1),
     $     xaver(2,nj,nkp1),yaver(2,nj,nkp1),zaver(2,nj,nkp1),
     $     binti(1,nj,nkp1),bintj(2,njm1,nkp1),bintk(2,nj,nk),
     $     volinv(njp1),fdum(1,njm1,nk)
#ifdef DEBUG
       dimension dumx(nj,nk,2),dumy(nj,nk,2),dumz(nj,nk,2)
#endif
      real nose,tail,ji,jj,jk,jx,jy,jz
!
! ... Parameter variables ..............................................
! 
      dimension bxx(nip1,njp1,nkp1)
      dimension byy(nip1,njp1,nkp1)
      dimension bzz(nip1,njp1,nkp1)
!
! ... Begin ............................................................
!
C
C  find the cell centered coordinates
C
#ifdef DEBUG_MODE_ON
      write(*,*) "DEBUG:  in ionosphere.F::jpara(...)"
#endif
#ifdef DEBUG
	write(9,*) ' JPARA: x2c(1,1),x2c(10,10)=',
     $    x2c(1,1),x2c(10,10)
	write(9,*) ' JPARA: y2c(1,1),y2c(10,10)=',
     $    y2c(1,1),y2c(10,10)
#endif
        do 200 i=1,2
        do 200 k=1,nk
        do 200 j=1,nj
          xc(i,j,k) = x2c(i,j)
          yc(i,j,k) = y2c(i,j)*cosarc2(k)*cos2phi(k)
          zc(i,j,k) = y2c(i,j)*cosarc2(k)*sin2phi(k)
  200   continue
C
        do 230 i=1,2
        do 230 j=1,nj
          xc(i,j,nkp1) = xc(i,j,1)
          yc(i,j,nkp1) = yc(i,j,1)
          zc(i,j,nkp1) = zc(i,j,1)
  230   continue
Cc
C       calculate the line integrals of b between cell centers
C
C
#ifdef DEBUG
       do 270 k=1,nk
       do 270 j=1,nj
       do 270 i=1,2
       dumx(j,k,i) = bxx(i,j,k)
       dumy(j,k,i) = byy(i,j,k)
       dumz(j,k,i) = bzz(i,j,k)
  270  continue
       call print2d(dumx,'bxx(1)',nj,nk,nj,nk)
       call print2d(dumx(1,1,2),'bxx(2)',nj,nk,nj,nk)
       call print2d(dumy,'byy(1)',nj,nk,nj,nk)
       call print2d(dumy(1,1,2),'byy(2)',nj,nk,nj,nk)
       call print2d(dumz,'bzz(1)',nj,nk,nj,nk)
       call print2d(dumz(1,1,2),'bzz(2)',nj,nk,nj,nk)
#endif
C
        do 300 i=1,1
        do 300 k=1,nk
        do 300 j=1,nj
         binti(i,j,k) = 0.5*(
     $          (bxx(i,j,k)+bxx(i+1,j,k))*(xc(i+1,j,k)-xc(i,j,k))
     $       +  (byy(i,j,k)+byy(i+1,j,k))*(yc(i+1,j,k)-yc(i,j,k))
     $       +  (bzz(i,j,k)+bzz(i+1,j,k))*(zc(i+1,j,k)-zc(i,j,k)))
  300      continue
C
        do 305 j=1,nj
         binti(1,j,nkp1) = binti(1,j,1)
  305   continue
C
        do 320 i=1,2
        do 320 k=1,nk
        do 320 j=1,njm1
         bintj(i,j,k) = 0.5*(
     $          (bxx(i,j,k)+bxx(i,j+1,k))*(xc(i,j+1,k)-xc(i,j,k))
     $       +  (byy(i,j,k)+byy(i,j+1,k))*(yc(i,j+1,k)-yc(i,j,k))
     $       +  (bzz(i,j,k)+bzz(i,j+1,k))*(zc(i,j+1,k)-zc(i,j,k)))
  320   continue
C
        do 325 i=1,2
        do 325 j=1,njm1
         bintj(i,j,nkp1) = bintj(i,j,1)
  325   continue
C
        do 340 i=1,2
        do 340 k=1,nkm1
        do 340 j=1,nj
         bintk(i,j,k) = 0.5*(
     $          (bxx(i,j,k)+bxx(i,j,k+1))*(xc(i,j,k+1)-xc(i,j,k))
     $       +  (byy(i,j,k)+byy(i,j,k+1))*(yc(i,j,k+1)-yc(i,j,k))
     $       +  (bzz(i,j,k)+bzz(i,j,k+1))*(zc(i,j,k+1)-zc(i,j,k)))
  340        continue
C
        do 342 i=1,2
        do 342 j=1,nj
         bintk(i,j,nk) = 0.5*(
     $          (bxx(i,j,nk)+bxx(i,j,1))*(xc(i,j,1)-xc(i,j,nk))
     $       +  (byy(i,j,nk)+byy(i,j,1))*(yc(i,j,1)-yc(i,j,nk))
     $       +  (bzz(i,j,nk)+bzz(i,j,1))*(zc(i,j,1)-zc(i,j,nk)))
  342        continue
C
C    find area integrals of j
C
#ifdef DEBUG
      call pictur(binti,'BNI',1,1,1,nj,1,nkp1,0,2,3,1,
     $     1,nj,nkp1,nj)
      call pictur(bintj,'BNj',1,2,1,njm1,1,nkp1,0,2,3,1,
     $     2,njm1,nkp1,2*njm1)
      call pictur(binti,'BNI',1,2,1,nj,1,nk,0,2,3,1,
     $     2,nj,nk,2*nj)
#endif

        do 400 i=1,1
        do 400 k=1,nkp1
        do 400 j=1,njm1
        jk(i,j,k) = binti(i,j,k)+bintj(i+1,j,k)-
     $          binti(i,j+1,k)-bintj(i,j,k)
  400   continue
C
        do 420 i=1,2
        do 420 k=1,nk
        do 420 j=1,njm1
        ji(i,j,k) = bintj(i,j,k)+bintk(i,j+1,k)-
     $          bintj(i,j,k+1)-bintk(i,j,k)
  420   continue
C
        do 440 i=1,1
        do 440 k=1,nk
        do 440 j=1,nj
        jj(i,j,k) = bintk(i,j,k)+binti(i,j,k+1)-
     $          bintk(i+1,j,k)-binti(i,j,k)
  440   continue
C
C    find the value of j over cell volume
C
        do 500 i=1,1
        do 500 k=1,nk
        do 500 j=1,njm1
  500   fdum(i,j,k) = ji(i,j,k)+ji(i+1,j,k)
C
        DO 510 I=1,2
        DO 510 K=1,NK
        DO 510 J=1,NJm1
        XAVER(I,J,K) = Xc(I,J,K)+Xc(I,J+1,K)+Xc(I,J,K+1)+Xc(I,J+1,K+1)
        YAVER(I,J,K) = Yc(I,J,K)+Yc(I,J+1,K)+Yc(I,J,K+1)+Yc(I,J+1,K+1)
  510   ZAVER(I,J,K) = Zc(I,J,K)+Zc(I,J+1,K)+Zc(I,J,K+1)+Zc(I,J+1,K+1)
C
        DO 520 I=1,1
        DO 520 K=1,NK
        DO 520 J=1,NJm1
        JX(I,J+1,k+1) = FDUM(I,J,K)*(XAVER(I+1,J,K)-XAVER(I,J,K))
        JY(I,J+1,k+1) = FDUM(I,J,K)*(YAVER(I+1,J,K)-YAVER(I,J,K))
  520   JZ(I,J+1,k+1) = FDUM(I,J,K)*(ZAVER(I+1,J,K)-ZAVER(I,J,K))
C
#ifdef DEBUG
      call pictur(jx,'jx1',1,1,1,njp1,1,nkp1,0,2,3,1,
     $     1,njp1,nkp1,njp1)
      call pictur(jy,'jy1',1,1,1,njp1,1,nkp1,0,2,3,1,
     $     1,njp1,nkp1,njp1)
      call pictur(jz,'jz1',1,1,1,njp1,1,nkp1,0,2,3,1,
     $     1,njp1,nkp1,njp1)
#endif
C
        DO 600 I=1,1
        DO 600 K=1,NK
        DO 600 J=1,NJm1
  600   FDUM(I,J,K)= jj(i,j,k)+jj(i,j+1,k)
C
        DO 610 I=1,1
        DO 610 K=1,NK
        DO 610 J=1,NJ
        XAVER(I,J,K) = Xc(I,J,K)+Xc(I+1,J,K)+Xc(I,J,K+1)+Xc(I+1,J,K+1)
        YAVER(I,J,K) = Yc(I,J,K)+Yc(I+1,J,K)+Yc(I,J,K+1)+Yc(I+1,J,K+1)
  610   ZAVER(I,J,K) = Zc(I,J,K)+Zc(I+1,J,K)+Zc(I,J,K+1)+Zc(I+1,J,K+1)
C
        DO 620 I=1,1
        DO 620 K=1,NK
        DO 620 J=1,NJm1
        JX(I,J+1,k+1)= JX(I,J+1,k+1)+FDUM(I,J,K)*
     $          (XAVER(I,J+1,K)-XAVER(I,J,K))
        JY(I,J+1,k+1)= JY(I,J+1,k+1)+FDUM(I,J,K)*
     $          (YAVER(I,J+1,K)-YAVER(I,J,K))
  620   JZ(I,J+1,k+1)= JZ(I,J+1,k+1)+FDUM(I,J,K)*
     $          (ZAVER(I,J+1,K)-ZAVER(I,J,K))
C
#ifdef DEBUG
      call pictur(jx,'jx2',1,1,1,njp1,1,nkp1,0,2,3,1,
     $     1,njp1,nkp1,njp1)
      call pictur(jy,'jy2',1,1,1,njp1,1,nkp1,0,2,3,1,
     $     1,njp1,nkp1,njp1)
      call pictur(jz,'jz2',1,1,1,njp1,1,nkp1,0,2,3,1,
     $     1,njp1,nkp1,njp1)
#endif
C
        DO 700 I=1,1
        DO 700 K=1,NK
        DO 700 J=1,NJm1
  700   FDUM(I,J,K)= jk(i,j ,k)+jk(i,j,k+1)
C
        DO 710 I=1,1
        DO 710 K=1,NKP1
        DO 710 J=1,NJm1
        XAVER(I,J,K) = Xc(I,J,K)+Xc(I+1,J,K)+Xc(I,J+1,K)+Xc(I+1,J+1,K)
        YAVER(I,J,K) = Yc(I,J,K)+Yc(I+1,J,K)+Yc(I,J+1,K)+Yc(I+1,J+1,K)
  710   ZAVER(I,J,K) = Zc(I,J,K)+Zc(I+1,J,K)+Zc(I,J+1,K)+Zc(I+1,J+1,K)
C
        DO 720 I=1,1
        DO 720 K=1,NK
        DO 720 J=1,NJm1
        volinv(j) = 0.125/(volq2(i,j)*cosarc2(k)*darc2(k))
        JX(I,J+1,k+1)= (JX(I,J+1,k+1)+FDUM(I,J,K)*
     $          (XAVER(I,J,K+1)-XAVER(I,J,K)))*volinv(j)
        JY(I,J+1,k+1)= (JY(I,J+1,k+1)+FDUM(I,J,K)*
     $          (YAVER(I,J,K+1)-YAVER(I,J,K)))*volinv(j)
  720   JZ(I,J+1,k+1)= (JZ(I,J+1,k+1)+FDUM(I,J,K)*
     $          (ZAVER(I,J,K+1)-ZAVER(I,J,K)))*volinv(j)
C
#ifdef DEBUG
        write (9,*) 'metric stuff in jpara'
        write (9,*) (k,cosarc2(k),darc2(k),k=1,nk)
        write (9,*) (j,volq2(1,j),j=1,njm1)
      call pictur(jx,'jx2',1,1,1,njp1,1,nkp1,0,2,3,1,
     $     1,njp1,nkp1,njp1)
      call pictur(jy,'jy2',1,1,1,njp1,1,nkp1,0,2,3,1,
     $     1,njp1,nkp1,njp1)
      call pictur(jz,'jz2',1,1,1,njp1,1,nkp1,0,2,3,1,
     $     1,njp1,nkp1,njp1)
#endif
C
        do 740 i=1,1
        do 740 k=1,nkp1
        jz(i,1,k) = 0.
        jx(i,1,k) = 0.
        jy(i,1,k) = 0.
        jz(i,njp1,k) = 0.
        jx(i,njp1,k) = 0.
        jy(i,njp1,k) = 0.
  740   continue
C
        do 760 i=1,1
        do 760 j=1,njp1
        jx(i,j,1) = jx(i,j,nkp1)
        jy(i,j,1) = jy(i,j,nkp1)
        jz(i,j,1) = jz(i,j,nkp1)
  760   continue
C
      CALL LFM_MI_CAPTURE_CURL(JX,JY,JZ)
C
C  now dot the j's with the right constant  to get j_parallel for the
C         field solve
        do 900 k=1,nkp1
        do 900 j=1,njp1
         current(j,k) = bxdot(j,k)*jx(1,j,k) + bydot(j,k)*jy(1,j,k) +
     $          bzdot(j,k)*jz(1,j,k)
  900   continue
C
C
#ifdef DEBUG
      write(9,*) ' JPARA: current(1,1)=',current(1,1)
#endif
!
! ... End jpara(bxx,byy,bzz) ...........................................
!
        return
        end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> jsetup
!!
      SUBROUTINE jsetup
!
! ... Local variables ..................................................
!
#include "ion-param.inc"
      REAL LFM_PLANET_RION,LFM_PLANET_FAC
#include "ion-meter.inc"
#include "ion_P++_pointers.inc"
#include "ion_fortran_P++_arrays.inc"
#include "jpara.inc"
c /* #include "iongeo.inc" */
c
      COMMON /IONSCRACH/
     $    DX(NI,NJP1),DY(NI,NJP1),DELTAX(NIP1,NJ),DELTAY(NIP1,NJ),
     $    DELDX(NI,NJ),DELDY(NI,NJ),
     $     J0(NI,NJ),JPSI(NI,NJ),JETA(NI,NJ),
     $    bsq(njp1)
#ifndef T3E
      double precision bsq, bxdott,bydott,bzdott
#endif
      REAL JXY
      REAL J0,JPSI,JETA
!
! ... Parameter variables ..............................................
!
! ... Begin ............................................................
!

#ifdef DEBUG_MODE_ON
      write(*,*) "DEBUG:  in ionosphere.F::jsetup(...)"
#endif
      mip1 = min(nip1,3)
      mi = mip1 - 1
      mim1 = mi - 1
      pi = atan(1.)*4.0
      do j=1,njp1
         do i=1,mip1
            X2(I,J) = x(i,j,1)
            y2(i,j) = sqrt(y(i,j,1)**2+z(i,j,1)**2)
         enddo
      enddo
      do k=1,nkp1
         cosphi(k) = y(1,2,k)/sqrt(y(1,2,k)**2+z(1,2,k)**2)
         sinphi(k) = z(1,2,k)/sqrt(y(1,2,k)**2+z(1,2,k)**2)
         phi(k) = atan2(z(1,2,k),y(1,2,k))
         if ( phi(k) .lt. 0. ) phi(k) = phi(k) + 2.*pi
      enddo
         phi(nkp1) = phi(1) + 2.*pi
         
            
        do 120 j=1,nj
        do 120 i=1,mi
        x2c(i,j) = 0.25*(x2(i,j)+x2(i+1,j)+x2(i,j+1)+x2(i+1,j+1))
        y2c(i,j) = 0.25*(y2(i,j)+y2(i+1,j)+y2(i,j+1)+y2(i+1,j+1))
  120   continue
#ifdef DEBUG
	write(9,*) ' JSETUP: x2c(1,1),x2c(10,10)=',
     $    x2c(1,1),x2c(10,10),((i,j,x2c(i,j),j=1,nj),i=1,mi)
	write(9,*) ' JSETUP: y2c(1,1),y2c(10,10)=',
     $    y2c(1,1),y2c(10,10),((i,j,y2c(i,j),j=1,nj),i=1,mi)
#endif
C
        DO 150 J=1,NJ
        DO 150 I=1,mIm1
        DX(I,J) = X2c(I+1,J)-X2c(I,J)
  150   DY(I,J) = Y2c(I+1,J)-Y2c(I,J)
C
C
C
        DO 200 J=1,NJm1
        DO 200 I=1,mI
        DELTAX(I,J) = X2c(I,J+1)-X2c(I,J)
  200   DELTAY(I,J) = Y2c(I,J+1)-Y2c(I,J)
        DO 210 J=1,NJm1
        DO 210 I=1,mIm1
        DELDX(I,J) = DX(I,J+1)-DX(I,J)
        DELDY(I,J) = DY(I,J+1)-DY(I,J)
        J0(I,J) = DX(I,J)*DELTAY(I,J) - DY(I,J)*DELTAX(I,J)
        JETA(I,J)= DX(I,J)*DELDY(I,J) - DY(I,J)*DELDX(I,J)
        JPSI(I,J) = DELTAY(I,J)*DELDX(I,J) - DELTAX(I,J)*DELDY(I,J)
           adum = J0(I,J)*(Y2c(I,J)+0.5*(DY(I,J)+DELTAY(I,J))+
     $                          0.25*DELDY(I,J)) 
          bdum =  JETA(I,J)*(0.5*Y2c(I,J)+ 0.333333*DY(I,J) +
     $       0.25*DELTAY(I,J)
     $                 + 0.16666666*DELDY(I,J)) 
          cdum =  JPSI(I,J)*(0.5*Y2c(I,J) + 0.25*DY(I,J)+
     $         0.333333*DELTAY(I,J)
     $                  + 0.16666666*DELDY(I,J))
        VOLq2(I,J) = adum + bdum + cdum
#ifdef DEBUG
        write (9,*) 'vlq2 ',i,j,adum,bdum,cdum,volq2(i,j)
#endif
  210   CONTINUE
#ifdef DEBUG
        call print2d(DX,'dx',mi,njp1,ni,njp1)
        call print2d(Dy,'dy',mi,njp1,ni,njp1)
        call print2d(DELTAX,'dltx',mip1,nj,nip1,nj)
        call print2d(DELTAy,'dlty',mip1,nj,nip1,nj)
        call print2d(DELDX,'ddx',mi,nj,ni,nj)
        call print2d(DELDY,'ddy',mi,nj,ni,nj)
        call print2d(J0,'j0',mi,nj,ni,nj)
        call print2d(Jpsi,'jpsi',mi,nj,ni,nj)
        call print2d(Jeta,'jeta',mi,nj,ni,nj)
        call print2d(volq2,'vlq2',nim1,njm1,nim1,njm1)
#endif
C
C
C
        do 450 k=1,nk
         phiq2(k) = 0.5*(phi(k)+phi(k+1))
  450   continue
         phiq2(nkp1) = phiq2(1) + 6.283185307
C
        DO 500 K=1,NK
        COS2PHI(K) = COS(PHIQ2(K))
        SIN2PHI(K) = SIN(PHIQ2(K))
        DARC2(K)   = 2.0*SIN(0.5*(PHIQ2(K+1)-PHIQ2(K)))
        COSARC2(K) = COS( 0.5*(PHIQ2(K+1)-PHIQ2(K)))
  500   CONTINUE
C
C
        do 600 k=1,nkp1
        do 600 j=1,njp1
          bxdot(j,k) = bxqq0(x(2,j,k),y(2,j,k),z(2,j,k))
          bydot(j,k) = byqq0(x(2,j,k),y(2,j,k),z(2,j,k))
          bzdot(j,k) = bzqq0(x(2,j,k),y(2,j,k),z(2,j,k))
  600   continue
        do 610 k=1,nkp1
        do 610 j=1,njp1
C  the constant gets the computed j to an ionospheric mks value
C    it is (5.5e-5tesla)*(100 cm/m)/mu_0
#ifndef CRAY /* CMM 16-aug-90 */
C     This version of the code uses Gaussian units for the magnetic field.
C     We have a floating point underflow here.
C     The current mips operating system just returns silent NaNs.
c      bxdott = max(sqrt(toosmall),abs(bxdot(j,k)))
c      bydott = max(sqrt(toosmall),abs(bydot(j,k)))
c      bzdott = max(sqrt(toosmall),abs(bzdot(j,k)))
c      bsq(j) = LFM_PLANET_FAC()/(bxdott**2+bydott**2+bzdott**2)
c
      bxdott = bxdot(j,k)
      bydott = bydot(j,k)
      bzdott = bzdot(j,k)
      bsq(j) = LFM_PLANET_FAC()/(bxdott**2+bydott**2+bzdott**2)
#else      
        bsq(j) = LFM_PLANET_FAC()/(bxdot(j,k)**2+bydot(j,k)**2+bzdot(j,k)**2)
#endif
        bxdot(j,k) = bxdot(j,k)*bsq(j)
        bydot(j,k) = bydot(j,k)*bsq(j)
        bzdot(j,k) = bzdot(j,k)*bsq(j)
  610   continue
!
! ... End SUBROUTINE jsetup.............................................
!
        RETURN
        END
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> map2ion(p,an)
!! @FIXME - optimize this code!
      subroutine map2ion(p,an)
!
! ... Local variables ..................................................
!
#include "ion-param.inc"
      REAL LFM_PLANET_RION,LFM_PLANET_FAC
c     For scalar quantities above the ionosphere.
!
! ... Parameter variables ..............................................
!
      real p(ni,nj,nk)
!> Description of p
      real an(njp1,nk2p3)
!> Description of "an"
!
! ... Begin ............................................................
!
#ifdef DEBUG_MODE_ON
      write(*,*) "DEBUG:  in ionosphere.F::map2ion(...)"
#endif
      do k=1,nk2
         km1=1 + mod(k+nk2-2,nk2)
         m=nkp1-k
         mm1=1 + mod(m+nk2-2,nk2)
         do j=2,nj
c           Bozo interpolation.
            an(j,k+1) = 0.0625*(
     $           p(2,j,k)+p(2,j-1,k)+p(2,j,km1)+p(2,j-1,km1)
     $          +p(2,j,m)+p(2,j-1,m)+p(2,j,mm1)+p(2,j-1,mm1)
     $          +p(1,j,k)+p(1,j-1,k)+p(1,j,km1)+p(1,j-1,km1)
     $          +p(1,j,m)+p(1,j-1,m)+p(1,j,mm1)+p(1,j-1,mm1))
         enddo
c     Bozo boundary conditions.
         j=1
            an(j,k+1) = 0.125*(
     $           p(2,j,k)+p(2,j,km1)+p(2,j,m)+p(2,j,mm1)+
     $           p(1,j,k)+p(1,j,km1)+p(1,j,m)+p(1,j,mm1))
         j=njp1
            an(j,k+1) = 0.125*(
     $           p(2,j-1,k)+p(2,j-1,km1)+p(2,j-1,m)+p(2,j-1,mm1)+
     $           p(1,j-1,k)+p(1,j-1,km1)+p(1,j-1,m)+p(1,j-1,mm1))
      enddo
      do j=1,njp1
         an(j,1) = 0.0
         an(j,2) = an(j,3)
         an(j,NK2p2) = an(j,nk2p1)
         an(j,nk2p3) = 0.0
      enddo
!
! ... End subroutine map2ion(p,an) .....................................
!
      return
      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> psicomp(psi, psii)
!!  This subroutine is only defined if the preprocessor flag
!!  "-DPSICOMPARE" is passed.
       subroutine psicomp(psi,psii)
!
! ... Local variables ..................................................
!
#include "ion-param.inc"
      REAL LFM_PLANET_RION,LFM_PLANET_FAC
#include "ion_P++_pointers.inc"
#include "ion_fortran_P++_arrays.inc"
#include "iongeo.inc"
       dimension zplot(10000)
       dimension xplot(10000)
       dimension yplot(10000)
       dimension dumqq(10000),idumqq(10000)
!
! ... Parameter variables ..............................................
!
       dimension psi(njp1,nk2p3)
       dimension psii(njp1,nk2p1)
!
! ... Begin ............................................................
!
#ifdef DEBUG_MODE_ON
      write(*,*) "DEBUG:  in ionosphere.F::psicomp(...)"
#endif
#ifdef PSICOMPARE
       do 200 k = 1,nk2p1
       do 200 j = 2, nj
         ind = j + (k-1)*(nj-1)
         xplot(ind) = x1ion(j,k)
         yplot(ind) = y1ion(j,k)
         zplot(ind) = psii(j,k)
  200  continue
       xplot(1) = x1ion(1,2)
       yplot(1) = y1ion(1,2)
       zplot(1) = psii(1,2)
       indtot = njp1 + nk2*(nj-1)
       xplot(indtot) = x1ion(njp1,2)
       yplot(indtot) = y1ion(njp1,2)
       zplot(indtot) = psii(njp1,2)
       call conraq(xplot,yplot,zplot,indtot,dumqq,idumqq)
       call frame
       do 300 k = 1,nk2p3
       do 300 j = 2, nj
         ind = j + (k-1)*(nj-1)
         xplot(ind) = x2ion(j,k)
         yplot(ind) = y2ion(j,k)
         zplot(ind) = psi(j,k)
  300  continue
       xplot(1) = x2ion(1,2)
       yplot(1) = y2ion(1,2)
       zplot(1) = psi(1,2)
       indtot = njp1 + (nk2+2)*(nj-1)
       xplot(indtot) = x2ion(njp1,2)
       yplot(indtot) = y2ion(njp1,2)
       zplot(indtot) = psi(njp1,2)
       call conraq(xplot,yplot,zplot,indtot,dumqq,idumqq)
       call frame
c
#endif PSICOMPARE
!
! ... End subroutine psicomp(psi,psii) .................................
!
       return
       end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> PSISLV(err, north)
!!
      SUBROUTINE PSISLV(err,north)
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
               sigmap(j,k) = ped0
               sigmah(j,k) = 0.0
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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
