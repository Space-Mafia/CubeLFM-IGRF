
      SUBROUTINE SUBWTs
#include "ion-param.inc"
#include "ion_P++_pointers.inc"
#include "ion_fortran_P++_arrays.inc"
#include "iongeo.inc"
#include "volts.inc"
c
      SAVE
      common /ionscrsch1/
     $     dxi(2,nj,nk2p3),dxj(2,njp1,nk2p2)
      dimension xdum(2)
      logical minusok,plusok,BAD
C
#ifdef DEBUG_MODE_ON
      write(*,*) "DEBUG:  in subwts.F::SUBWTs(...)"
#endif
        do 200 k=1,nk2p3
        do 200 j=1,nj
          dxi(1,j,k) = x2ion(j+1,k)-x2ion(j,k)
          dxi(2,j,k) = y2ion(j+1,k)-y2ion(j,k)
  200   continue
        do 201 k=1,nk2p2
        do 201 j=1,njp1
          dxj(1,j,k) = x2ion(j,k+1)-x2ion(j,k)
          dxj(2,j,k) = y2ion(j,k+1)-y2ion(j,k)
  201   continue
        do 202 k=2,nk2p1
          dxj(1,1,k) = 0.
          dxj(2,1,k) = -1.e-6
          dxj(1,njp1,k) = 0.
          dxj(2,njp1,k) = -1.e-6
  202   continue
C
C   loop over inner sphere
C
        iint = 1
        jint = 1
        do 401 k=2,nk2
        if ( k .eq. nk2/2 + 1 ) then
          go to 401
        else if ( k .le. nk2/2 ) then
           koff = 1
        else
           koff = nk2p3
        endif
        do 400 j=2,nj
C
c check to see if inside the grid
c
           xd = x1ion(j,k)
           idj = 0
	if ( xd .le. x2ion(1,koff) .and. xd .gt. x2ion(njp1,koff)) then
           do 388 jjj=1,nj
              if ( x1ion(j,k) .le. x2ion(jjj,koff) .and.
     $             x1ion(j,k) .gt. x2ion(jjj+1,koff) )  then
                 idj = jjj
                 delt=(x1ion(j,k)-x2ion(jjj,koff))/(x2ion(jjj+1,koff)-
     $                x2ion(jjj,koff))
                 yd = y2ion(jjj,koff)*(1.-delt)+y2ion(jjj+1,koff)*delt
                 if ( ((k .le. nk2/2) .and. ( yd .le. y1ion(j,k))) .or.
     $               ((k .gt. nk2/2) .and. ( yd .ge. y1ion(j,k)))) then
                      idj = 0
                 endif
                 go to 389
              endif
 388          continue
 389          continue
         endif
c
         if ( idj .eq. 0 ) then
            eta1(j,k) =0.0
            chi1(j,k) = 0.0
         else
            iint = idj
            jint = k
c
        do 399 kounter=1,200
          iold = iint
          jold = jint
          xdum(1) = x1ion(j,k) - x2ion(iint,jint)
          xdum(2) = y1ion(j,k) - y2ion(iint,jint)
          si = sign(0.5,cross(dxi(1,iint,jint),xdum))
          sj = sign(0.5,cross(xdum,dxj(1,iint,jint)))
          if ( si + sj .gt. 0.75 ) then
             xdum(1) = x1ion(j,k) - x2ion(iint+1,jint+1)
             xdum(2) = y1ion(j,k) - y2ion(iint+1,jint+1)
             ssi = sign(0.5,cross(xdum,dxi(1,iint,jint+1)))
             ssj = sign(0.5,cross(dxj(1,iint+1,jint),xdum))
             if ( ssi + ssj .gt. 0.75 ) then
                eta1(j,k) = float(iint)+0.5
                chi1(j,k) = float(jint)+0.5
                go to 390   !success
             else
                if ( ssi .lt. 0. ) jint = max0(1,min0(nk2p2,jint+1))
                if ( ssj .lt. 0. ) then
                   if ( iint .eq. nj .and. jold .ne. 1 ) then
                        iint = nj
                        jint = 1
                   elseif ( iint .eq. nj .and. jold .eq. 1 ) then
                        iint = nj
                        jint = nk2p2
                   else
                        iint = max0(1,min0(nj,iint+1))
                   endif
                endif
            endif
          else
                if ( si .lt. 0. ) jint = max0(1,min0(nk2p2,jint-1))
                if ( sj .lt. 0. ) then
                   if (iint .eq. 1 .and. jold .ne. 1 ) then
                        iint = 1
                        jint = 1
                   elseif ( iint .eq. 1 .and. jold .eq. 1 ) then
                        iint = 1
                        jint = nk2p2
                   else
                        iint = max0(1,min0(nj,iint-1))
                   endif
                endif
          endif
          if ( iold .eq. iint .and. jold .eq. jint ) go to 397
  399   continue
C
          write(9,*)
     $       'tried maximum number of times but could not find'
     $          ,'  cell, j,k = ',j,k
        go to 395
  397   continue
             write (9,*)
     $       'not in cell, but no update, i,j = ',iint,jint
     $                          ,' i1,j1 =  ',j,k
  395   continue
             eta1(j,k) = float(iint)+0.5
             chi1(j,k) = float(jint)+0.5
 390    continue
        endif
  400   continue
  401   continue
C
C       figure out where in the cell we are
C
        do 451 k=2,nk2
        if ( k .eq. nk2/2 + 1 ) go to 451
        do 450 j=2,nj    ! only interior points needed
        if ( eta1(j,k) .ge. 1.0 ) then
          ii = eta1(j,k)
          jj = chi1(j,k)
C                       regular quadrilaterals
          dx = x1ion(j,k) - x2ion(ii,jj)
          dy = y1ion(j,k) - y2ion(ii,jj)
	  xeta = dxi(1,ii,jj)
	  xchi = dxj(1,ii,jj)
	  xetachi = dxj(1,ii+1,jj)-dxj(1,ii,jj)
	  yeta = dxi(2,ii,jj)
	  ychi = dxj(2,ii,jj)
	  yetachi = dxj(2,ii+1,jj)-dxj(2,ii,jj)
          aa = yeta*xetachi-xeta*yetachi
	  if ( aa .ne. 0.0 ) then
	    bb = xchi*yeta-xeta*ychi +yetachi*dx - xetachi*dy
	    cc = ychi*dx - xchi*dy
	    disc = sqrt(bb**2-4*aa*cc)/(2.*aa)
	    bt = -bb/(2.*aa)
	    etap = bt+disc
	    plusok = etap .ge. 0. .and. etap .le. 1.0
	    etam = bt-disc
	    minusok = etam .ge. 0. .and. etam .le. 1.0
	    if ( minusok .and. plusok ) then
		  etainc = etap
		  if ( abs(0.5-etap) .gt. abs(0.5-etam) )
     $			etainc = etam
		  chiinc = (dy -yeta*etainc)/(ychi+yetachi*etainc)
	          if ( chiinc .gt. 1.0 .or. chiinc .lt. 0.0 ) then
		     etainc=etap
		     chiinc = (dy -yeta*etainc)/(ychi+yetachi*etainc)
	          endif
	    elseif ( plusok ) then
		  etainc = etap
		  chiinc = (dy -yeta*etainc)/(ychi+yetachi*etainc)
	    else
		  etainc = etam
		  chiinc = (dy -yeta*etainc)/(ychi+yetachi*etainc)
	    endif
	    eta1(j,k) = float(ii) + etainc
	    chi1(j,k) = float(jj) + chiinc
	    if ( etainc .gt. 1. .or. etainc .lt. 0. .or.
     $		  chiinc .gt. 1. .or. chiinc .lt. 0. ) then
		  et1 = bt+disc
	    ch1 = (dy -yeta*et1)/(ychi+yetachi*et1)
		  et2 = bt-disc
	    ch2 = (dy -yeta*et2)/(ychi+yetachi*et2)
c   check if only slightly wrong
	     errmax = 3.e-5
	   if ( (et1 .le. 1.+errmax .and. et1 .ge. -errmax) .and.
     $		(ch1 .le. 1.+errmax .and. ch1 .ge. -errmax) ) then
	     eta1(j,k) = float(ii) + amax1(0.,amin1(1.,et1))
	     chi1(j,k) = float(jj) + amax1(0.,amin1(1.,ch1))
	     BAD = .false.
	   elseif ( (et2 .le. 1.+errmax .and. et2 .ge. -errmax) .and.
     $		(ch2 .le. 1.+errmax .and. ch2 .ge. -errmax) ) then
	     eta1(j,k) = float(ii) + amax1(0.,amin1(1.,et2))
	     chi1(j,k) = float(jj) + amax1(0.,amin1(1.,ch2))
	     BAD = .false.
	   else
	    BAD = .true.
	       write (6,*) 'bad solution for position within cell',
     $		  '   xion,yion,i2,j2,eta,chi  ',xion,yion,ii,jj,
     $	     	    etainc,chiinc
	       write (6,*) 'bt,disc,+ and -',bt,disc,et1,ch1,et2,ch2
	   endif
	     endif
	   else
	     BAD = .true.
*	     write (6,*) 'stuck in a bad place -- xion,yion,iint,jint',
*     $	     	xion,yion,iint,jint
	   endif
           endif
  450   continue
  451   continue
C
        do 455 j=1,njp1
        chi1(j,1) = 1.0000001
        chi1(j,nk2p1) = float(nk2p2) + 0.99999
        chi1(j,nk2/2+1) = float(nk2/2+2)+0.000001
  455   continue
        do 456 k=2,nk2
        eta1(1,k) = 1.000001
        chi1(1,k) = 1.000001
        eta1(njp1,k) = float(nj) + 0.99999
        chi1(njp1,k) = 1.000001
  456   continue
        do 457 k=1,nk2p1,nk2/2
          kp = k + 2*k/nk2
        eta1(1,k) = 1.000001
        eta1(njp1,k) = float(nj) + 0.99999
          do 458 j=2,nj
            do 459 jj=1,nj
              if ( x1ion(j,k) .ge. x2ion(jj+1,kp) .aND.
     $             x1ion(j,k) .lt. x2ion(jj,kp)) then
                eta1(j,k) = float(jj) + (x1ion(j,k)-x2ion(jj,kp))/
     $                                  (x2ion(jj+1,kp)-x2ion(jj,kp))
                go to 458
              endif
  459       continue
            if ( k .eq. nk2/2+1 .and.
     $          x1ion(j,k) .gt. x2ion(1,kp) .and.
     $          x1ion(j,k) .le. x2ion(1,1) ) then
              eta1(j,k) = 1.000001
              chi1(j,k) = 1. +
     $          (x1ion(j,k)-x2ion(1,1))/(x2ion(1,2)-x2ion(1,1))
            elseif ( k .eq. nk2/2+1 .and.
     $          x1ion(j,k) .gt. x2ion(njp1,1) .and.
     $          x1ion(j,k) .le. x2ion(njp1,kp) ) then
              eta1(j,k) = float(nj) + 0.99999
              chi1(j,k) = 1. +
     $          (x1ion(j,k)-x2ion(njp1,1))/(x2ion(njp1,2)-x2ion(njp1,1))
            else
              write (9,*) 'could not find intersection for j,k',j,k
              eta1(j,k) = float(nj)
            endif
  458     continue
  457   continue
C
C
C   loop over first computational cell
C
        iint = 1
        jint = 1
        do 601 k=1,nk2
        if ( k .le. nk2/2 ) then
           koff = 1
        else
           koff = nk2p3
        endif
        do 600 j=1,nj
C
c check to see if inside the grid
c
           xd = xhfion(j,k)
           idj = 0
        if (xd .le. x2ion(1,koff) .and. xd .gt. x2ion(njp1,koff))then
           do 588 jjj=1,nj
              if ( xhfion(j,k) .le. x2ion(jjj,koff) .and.
     $             xhfion(j,k) .gt. x2ion(jjj+1,koff) )  then
                 idj = jjj
                delt=(xhfion(j,k)-x2ion(jjj,koff))/(x2ion(jjj+1,koff)-
     $                x2ion(jjj,koff))
                 yd = y2ion(jjj,koff)*(1.-delt)+y2ion(jjj+1,koff)*delt
                 if ( ((k .le. nk2/2) .and. ( yd .le. yhfion(j,k))) .or.
     $               ((k .gt. nk2/2) .and. ( yd .ge. yhfion(j,k)))) then
                      idj = 0
                 endif
                 go to 589
              endif
 588          continue
 589          continue
         endif
c
         if ( idj .eq. 0 ) then
            eta1(j,k) =0.0
            chi1(j,k) = 0.0
         else
            iint = idj
            jint = k
c
C
        do 599 kounter=1,200
          iold = iint
          jold = jint
          xdum(1) = xhfion(j,k) - x2ion(iint,jint)
          xdum(2) = yhfion(j,k) - y2ion(iint,jint)
          si = sign(0.5,cross(dxi(1,iint,jint),xdum))
          sj = sign(0.5,cross(xdum,dxj(1,iint,jint)))
          if ( si + sj .gt. 0.75 ) then
             xdum(1) = xhfion(j,k) - x2ion(iint+1,jint+1)
             xdum(2) = yhfion(j,k) - y2ion(iint+1,jint+1)
             ssi = sign(0.5,cross(xdum,dxi(1,iint,jint+1)))
             ssj = sign(0.5,cross(dxj(1,iint+1,jint),xdum))
             if ( ssi + ssj .gt. 0.75 ) then
                etahf(j,k) = float(iint)+0.5
                chihf(j,k) = float(jint)+0.5
                go to 590   !success
             else
                if ( ssi .lt. 0. ) jint = max0(1,min0(nk2p2,jint+1))
                if ( ssj .lt. 0. ) then
                   if ( iint .eq. nj .and. jold .ne. 1 ) then
                        iint = nj
                        jint = 1
                   elseif ( iint .eq. nj .and. jold .eq. 1 ) then
                        iint = nj
                        jint = nk2p2
                   else
                        iint = max0(1,min0(nj,iint+1))
                   endif
                endif
           endif
          else
                if ( si .lt. 0. ) jint = max0(1,min0(nk2p2,jint-1))
                if ( sj .lt. 0. ) then
                   if (iint .eq. 1 .and. jold .ne. 1 ) then
                        iint = 1
                        jint = 1
                   elseif ( iint .eq. 1 .and. jold .eq. 1 ) then
                        iint = 1
                        jint = nk2p2
                   else
                        iint = max0(1,min0(nj,iint-1))
                   endif
                endif
           endif
          if ( iold .eq. iint .and. jold .eq. jint ) go to 597
  599   continue
C
          write (9,*)
     $       'tried maximum number of times but could not find',
     $          '  cell, j,k = ',j,k
        go to 595
  597   continue
             write (9,*) 'not in cell, but no update, i,j = ',iint,jint
     $                          ,' ihf,jhf =  ',j,k
  595   continue
             etahf(j,k) = float(iint)+0.5
             chihf(j,k) = float(jint)+0.5
 590    continue
        endif
  600   continue
  601   continue
C
C       figure out where in the cell we are
C
*	call print2d(x2ion,'x2ion',njp1,nk2p3,njp1,nk2p3)
*	call print2d(y2ion,'y2ion',njp1,nk2p3,njp1,nk2p3)
*	call print2d(xhfion,'xhfion',nj,nk2,nj,nk2)
*	call print2d(yhfion,'yhfion',nj,nk2,nj,nk2)
        do 650 k=1,nk2
        do 650 j=1,nj    ! only interior points needed
        if ( etahf(j,k) .ge. 1.0 ) then
          ii = etahf(j,k)
          jj = chihf(j,k)
C                       regular quadrilaterals
          dx = xhfion(j,k) - x2ion(ii,jj)
          dy = yhfion(j,k) - y2ion(ii,jj)
	  xeta = dxi(1,ii,jj)
	  xchi = dxj(1,ii,jj)
	  xetachi = dxj(1,ii+1,jj)-dxj(1,ii,jj)
	  yeta = dxi(2,ii,jj)
	  ychi = dxj(2,ii,jj)
	  yetachi = dxj(2,ii+1,jj)-dxj(2,ii,jj)
*	write (6,*) dx,dy,xeta,xchi,xetachi,yeta,ychi,yetachi
          aa = yeta*xetachi-xeta*yetachi
	  if ( aa .ne. 0.0 ) then
	    bb = xchi*yeta-xeta*ychi +yetachi*dx - xetachi*dy
	    cc = ychi*dx - xchi*dy
	    disc = sqrt(bb**2-4*aa*cc)/(2.*aa)
	    bt = -bb/(2.*aa)
	    etap = bt+disc
	    plusok = etap .ge. 0. .and. etap .le. 1.0
	    etam = bt-disc
	    minusok = etam .ge. 0. .and. etam .le. 1.0
	    if ( minusok .and. plusok ) then
		  etainc = etap
		  if ( abs(0.5-etap) .gt. abs(0.5-etam) )
     $			etainc = etam
		  chiinc = (dy -yeta*etainc)/(ychi+yetachi*etainc)
	          if ( chiinc .gt. 1.0 .or. chiinc .lt. 0.0 ) then
		     etainc=etap
		     chiinc = (dy -yeta*etainc)/(ychi+yetachi*etainc)
	          endif
	    elseif ( plusok ) then
		  etainc = etap
		  chiinc = (dy -yeta*etainc)/(ychi+yetachi*etainc)
	    else
		  etainc = etam
		  chiinc = (dy -yeta*etainc)/(ychi+yetachi*etainc)
	    endif
	    etahf(j,k) = float(ii) + etainc
	    chihf(j,k) = float(jj) + chiinc
	    if ( etainc .gt. 1. .or. etainc .lt. 0. .or.
     $		  chiinc .gt. 1. .or. chiinc .lt. 0. ) then
		  et1 = bt+disc
	    ch1 = (dy -yeta*et1)/(ychi+yetachi*et1)
		  et2 = bt-disc
	    ch2 = (dy -yeta*et2)/(ychi+yetachi*et2)
c   check if only slightly wrong
	     errmax = 3.e-5
	   if ( (et1 .le. 1.+errmax .and. et1 .ge. -errmax) .and.
     $		(ch1 .le. 1.+errmax .and. ch1 .ge. -errmax) ) then
	     etahf(j,k) = float(ii) + amax1(0.,amin1(1.,et1))
	     chihf(j,k) = float(jj) + amax1(0.,amin1(1.,ch1))
	     BAD = .false.
	   elseif ( (et2 .le. 1.+errmax .and. et2 .ge. -errmax) .and.
     $		(ch2 .le. 1.+errmax .and. ch2 .ge. -errmax) ) then
	     etahf(j,k) = float(ii) + amax1(0.,amin1(1.,et2))
	     chihf(j,k) = float(jj) + amax1(0.,amin1(1.,ch2))
	     BAD = .false.
	   else
	    BAD = .true.
	       write (6,*) 'bad solution for position within cell',
     $		  '   xion,yion,i2,j2,eta,chi  ',xion,yion,ii,jj,
     $	     	    etainc,chiinc
	       write (6,*) 'bt,disc,+ and -',bt,disc,et1,ch1,et2,ch2
	   endif
	     endif
	   else
	     BAD = .true.
*	     write (6,*) 'stuck in a bad place -- xion,yion,iint,jint',
*     $	     	xion,yion,iint,jint
	   endif
           endif
  650   continue
C
c   now let's set those things that are outside the boundary equal to
c     something fairly safe like a point on the edge = (1,1)
        do 700 k=1,nk2p1
           do 700 j=1,njp1
              if ( eta1(j,k) .lt. 0.5 ) then
                 eta1(j,k) = 1.000001
                 chi1(j,k) = 1.000001
              endif
 700    continue
C
        do 710 k=1,nk2
           do 710 j=1,nj
              if ( etahf(j,k) .lt. 0.5 ) then
                 etahf(j,k) = 1.000001
                 chihf(j,k) = 1.000001
              endif
 710    continue
C
        RETURN
        END

        subroutine vcoef(xahf,yahf,zahf)
#include "ion-param.inc"
#include "ion_P++_pointers.inc"
#include "ion_fortran_P++_arrays.inc"
#include "iongeo.inc"
#include "volts.inc"
#include "dipole.inc"
#include "ionconst.inc"
C     We need geoqmu.
C
        dimension xahf(nj,nk2),yahf(nj,nk2),zahf(nj,nk2)
*        dimension bxhf(nj,nk2),byhf(nj,nk2),bzhf(nj,nk2)
*        dimension psiteta(nj,nk2),psiphi(nj,nk2)
*        dimension psix(nj,nk2),psiy(nj,nk2),vdotb(nj,nk2)
*        dimension vr(nj,nk2),vtet(nj,nk2),vphi(nj,nk2)
*        dimension ex(nj,nk2),ey(nj,nk2),ez(nj,nk2)
*        dimension psihf(nj,nk2),dpsi(nj,nk2)

C
      SAVE
      common /ionscrsch1/
     $     rpr(nj,nk2),rsq(nj,nk2),costhpr(nj,nk2),
     $     sinthpr(nj,nk2),cosph(nj,nk2),sinph(nj,nk2),
     $     rot(3,3,nj,nk2),eul(3,2,nj,nk2),der(2,2,nj,nk2),
     $     vspher(3,2,nj,nk2)
        geomu = geoqmu*1.e-10
C
#ifdef DEBUG_MODE_ON
      write(*,*) "DEBUG:  in subwts.F::vcoef(...)"
#endif
        do 200 k=1,nk2
        do 200 j=1,nj
          rpr(j,k) = sqrt(xahf(j,k)**2+yahf(j,k)**2+zahf(j,k)**2)
          rsq(j,k) = rpr(j,k)**2
          costhpr(j,k) = zahf(j,k)/rpr(j,k)
          sinthpr(j,k) = sqrt(1. - costhpr(j,k)**2)
          cosph(j,k) = cos(phihf(j,k))
          sinph(j,k) = sin(phihf(j,k))
  200   continue
C
        do 500 k=1,nk2
        do 500 j=1,nj
          rot(1,1,j,k) = sinthpr(j,k)*cosph(j,k)
          rot(1,2,j,k) = costhpr(j,k)*cosph(j,k)
          rot(1,3,j,k) = -sinph(j,k)
          rot(2,1,j,k) = sinthpr(j,k)*sinph(j,k)
          rot(2,2,j,k) = costhpr(j,k)*sinph(j,k)
          rot(2,3,j,k) = cosph(j,k)
          rot(3,1,j,k) = costhpr(j,k)
          rot(3,2,j,k) = -sinthpr(j,k)
          rot(3,3,j,k) = 0.
C
C   now the Euler potential derivatives
C
          denom = 1./(geomu*(1.+3.*costhpr(j,k)**2))
          Eul(1,1,j,k) = yhfion(j,k)*rsq(j,k)*denom
          Eul(1,2,j,k) = -xhfion(j,k)*rsq(j,k)*denom
          Eul(2,1,j,k) = -2.*Eul(1,1,j,k)*costhpr(j,k)/sinthpr(j,k)
          Eul(2,2,j,k) = -2.*Eul(1,2,j,k)*costhpr(j,k)/sinthpr(j,k)
          denom = rpr(j,k)*sinthpr(J,k)/
     $          (geomu/rionm*tetahf(j,k)*sin(2.*tetahf(j,k)))
          Eul(3,1,j,k) = -xhfion(j,k)*denom
          Eul(3,2,j,k) = -yhfion(j,k)*denom
C
C   now the loCAL DERIVATIVES
C
	if ( etahf(j,k) .lt. 1.002 ) then
	  denom = 1.0
	  xeta = 0.0
	  yeta = 0.0
	  xchi = 0.0
	  ychi = 0.0
	else
          ii = etahf(j,k)
          jj = chihf(j,k)
          eta = etahf(j,k)-float(ii)
          chi = chihf(j,k) - float(jj)
          xchi = (x2ion(ii,jj+1)-x2ion(ii,jj))*(1.-eta) +
     $           (x2ion(ii+1,jj+1)-x2ion(ii+1,jj))*eta
          ychi = (y2ion(ii,jj+1)-y2ion(ii,jj))*(1.-eta) +
     $           (y2ion(ii+1,jj+1)-y2ion(ii+1,jj))*eta
          xeta = (x2ion(ii+1,jj)-x2ion(ii,jj))*(1.-chi) +
     $           (x2ion(ii+1,jj+1)-x2ion(ii,jj+1))*chi
          yeta = (y2ion(ii+1,jj)-y2ion(ii,jj))*(1.-chi) +
     $           (y2ion(ii+1,jj+1)-y2ion(ii,jj+1))*chi
          denom = 1./(ychi*xeta-yeta*xchi)
	endif
          Der(1,1,j,k) = ychi*denom
          Der(1,2,j,k) = -yeta*denom
          Der(2,1,j,k) = -xchi*denom
          Der(2,2,j,k) = xeta*denom
C
C   multiply them all together
C
          Vspher(1,1,j,k) = Eul(1,1,j,k)*Der(1,1,j,k) +
     $           Eul(1,2,j,k)*Der(2,1,j,k)
          Vspher(1,2,j,k) = Eul(1,1,j,k)*Der(1,2,j,k) +
     $           Eul(1,2,j,k)*Der(2,2,j,k)
          Vspher(2,1,j,k) = Eul(2,1,j,k)*Der(1,1,j,k) +
     $           Eul(2,2,j,k)*Der(2,1,j,k)
          Vspher(2,2,j,k) = Eul(2,1,j,k)*Der(1,2,j,k) +
     $           Eul(2,2,j,k)*Der(2,2,j,k)
          Vspher(3,1,j,k) = Eul(3,1,j,k)*Der(1,1,j,k) +
     $           Eul(3,2,j,k)*Der(2,1,j,k)
          Vspher(3,2,j,k) = Eul(3,1,j,k)*Der(1,2,j,k) +
     $           Eul(3,2,j,k)*Der(2,2,j,k)
C
          vxeta(j,k) = ( rot(1,1,j,k)*vspher(1,1,j,k) +
     $           rot(1,2,j,k)*vspher(2,1,j,k) +
     $          rot(1,3,j,k)*vspher(3,1,j,k))*100.
          vyeta(j,k) = ( rot(2,1,j,k)*vspher(1,1,j,k) +
     $           rot(2,2,j,k)*vspher(2,1,j,k) +
     $          rot(2,3,j,k)*vspher(3,1,j,k))*100.
          vzeta(j,k) = ( rot(3,1,j,k)*vspher(1,1,j,k) +
     $           rot(3,2,j,k)*vspher(2,1,j,k) +
     $          rot(3,3,j,k)*vspher(3,1,j,k))*100.
          vxchi(j,k) = ( rot(1,1,j,k)*vspher(1,2,j,k) +
     $           rot(1,2,j,k)*vspher(2,2,j,k) +
     $          rot(1,3,j,k)*vspher(3,2,j,k))*100.
          vychi(j,k) = ( rot(2,1,j,k)*vspher(1,2,j,k) +
     $           rot(2,2,j,k)*vspher(2,2,j,k) +
     $          rot(2,3,j,k)*vspher(3,2,j,k))*100.
          vzchi(j,k) = ( rot(3,1,j,k)*vspher(1,2,j,k) +
     $           rot(3,2,j,k)*vspher(2,2,j,k) +
     $          rot(3,3,j,k)*vspher(3,2,j,k))*100.
C
C
  500   continue
*       WRITE (9,*) 'xahf,yahf,zahf,tetahf,phihf(3,4)',
*     $         xahf(3,4),yahf(3,4),zahf(3,4),tetahf(3,4),phihf(3,4)
*       write (9,*) 'vspher'
*       write (9,*) vspher(1,1,3,4),vspher(1,2,3,4)
*       write (9,*) vspher(2,1,3,4),vspher(2,2,3,4)
*       write (9,*) vspher(3,1,3,4),vspher(3,2,3,4)
*       write (9,*) ((rot(i,j,3,4),i=1,3),j=1,3)
*       write (9,*) ((eul(i,j,3,4),i=1,2),j=1,2)
*c
*c    test the sucker out
*c
*       do 590 k=1,nk2
*       do 590 j=1,nj
*         bxhf(j,k) = bxqq0(100.*xahf(j,k),100.*yahf(j,k)
*     $         ,100.*zahf(j,k))
*         byhf(j,k) = byqq0(100.*xahf(j,k),100.*yahf(j,k)
*     $         ,100.*zahf(j,k))
*         bzhf(j,k) = bzqq0(100.*xahf(j,k),100.*yahf(j,k)
*     $         ,100.*zahf(j,k))
*  590  continue
*c
*       do 600 k=1,nk2p3
*       do 600 j=1,njp1
**        psi(j,k) = sin(teta2(j,k))
*         psi(j,k) = sin(phi2(j,k))
*  600  continue
*c
*       call intrp
*c
*       do 610 k=1,nk2
*       do 610 j=1,nj
*         psix(j,k) = psieta(j,k)*der(1,1,j,k)+psichi(j,k)*
*     $         der(1,2,j,k)
*         psiy(j,k) = psieta(j,k)*der(2,1,j,k)+psichi(j,k)*
*     $         der(2,2,j,k)
*         psiteta(j,k) = (xhfion(j,k)*psix(j,k)+yhfion(j,k)*psiy(j,k))/
*     $         (tetahf(j,k)*cos(tetahf(j,k)))
*         psiphi(j,k) = -(yhfion(j,k)*psix(j,k)-xhfion(j,k)*psiy(j,k))
*         vr(j,k) = vspher(1,1,j,k)*psieta(j,k) +
*     $             vspher(1,2,j,k)*psichi(j,k)
*         vtet(j,k) = vspher(2,1,j,k)*psieta(j,k) +
*     $             vspher(2,2,j,k)*psichi(j,k)
*         vphi(j,k) = vspher(3,1,j,k)*psieta(j,k) +
*     $             vspher(3,2,j,k)*psichi(j,k)
*  610  continue
*c
*       call print2d(psix,'px',nj,nk2,nj,nk2)
*       call print2d(psiy,'py',nj,nk2,nj,nk2)
*       call print2d(psiteta,'ptet',nj,nk2,nj,nk2)
*       call print2d(psiphi,'pphi',nj,nk2,nj,nk2)
*       call print2d(vr,'vr',nj,nk2,nj,nk2)
*       call print2d(vtet,'vtet',nj,nk2,nj,nk2)
*       call print2d(vphi,'vphi',nj,nk2,nj,nk2)
*       call print2d(vxint,'vxint',nj,nk2,nj,nk2)
*       call print2d(vyint,'vyint',nj,nk2,nj,nk2)
*       call print2d(vzint,'vzint',nj,nk2,nj,nk2)
*c
*       do 650 k=1,nk2
*       do 650 j=1,nj
*         vdotb(j,k) = (vxint(j,k)*bxhf(j,k)+vyint(j,k)*byhf(j,k)
*     $         +vzint(j,k)*bzhf(j,k))/sqrt(
*     $     (vxint(j,k)**2+vyint(j,k)**2+vzint(j,k)**2)*
*     $     (bxhf(j,k)**2+byhf(j,k)**2+bzhf(j,k)**2) )
*  650  continue
*c
*       call print2d(vdotb,'vdotb',nj,nk2,nj,nk2)
*c
*       do 660 k=1,nk2
*       do 660 j=1,nj
*         ii = etahf(j,k)
*         jj = chihf(j,k)
*         eta = etahf(j,k)-float(ii)
*         chi = chihf(j,k) - float(jj)
*         psihf(j,k) = psi(ii,jj)*(1.-eta)*(1.-chi) +
*     $                 psi(ii+1,jj)*eta*(1.-chi)   +
*     $                 psi(ii,jj+1)*(1.-eta)*chi   +
*     $                 psi(ii+1,jj+1)*eta*chi
*         Ex(j,k) = vzint(j,k)*byhf(j,k)-vyint(j,k)*bzhf(j,k)
*         ey(j,k) = vxint(j,k)*bzhf(j,k)-vzint(j,k)*bxhf(j,k)
*         ez(j,k) = vyint(j,k)*bxhf(j,k)-vxint(j,k)*byhf(j,k)
*  660  continue
*c
*       do 670 k=1,nk2
*       do 670 j=1,nj-1
*        dpsi(j,k) = ( 0.5*(xahf(j+1,k)-xahf(j,k))*(ex(j,k)+ex(j+1,k)) +
*     $         0.5*(yahf(j+1,k)-yahf(j,k))*(ey(j,k)+ey(j+1,k)) +
*     $         0.5*(zahf(j+1,k)-zahf(j,k))*(ez(j,k)+ez(j+1,k)) )/
*     $         (psihf(j+1,k)-psihf(j,k))
*  670  continue
*c
*       call print2d(dpsi,'dpsi',nj,nk2,nj,nk2)
*c
        return
        end

      function fchi(f,l,m,i,j)
#include "ion-param.inc"
#include "gcoef.inc"
      dimension f(njp1,nk2p3)
C
#ifdef DEBUG_MODE_ON
      write(*,*) "DEBUG:  in subwts.F::fchi(...)"
#endif
      fchi = f(i,j)*echi(l,m,1,1)+f(i+1,j)*echi(l,m,2,1)
     $     +  f(i,j+1)*echI(l,m,3,1)+f(i+1,j+1)*echi(l,m,4,1)
C
      return
      end
      function feta(f,l,m,i,j)
#include "ion-param.inc"
#include "gcoef.inc"
      dimension f(njp1,nk2p3)
#ifdef DEBUG_MODE_ON
      write(*,*) "DEBUG:  in subwts.F::feta(...)"
#endif
      feta = f(i,j)*eeta(l,m,1,1)+f(i+1,j)*eeta(l,m,2,1)
     $     +  f(i,j+1)*eeta(l,m,3,1)+f(i+1,j+1)*eeta(l,m,4,1)
C
      return
      end
      function fint(f,l,m,i,j)
#include "ion-param.inc"
#include "gcoef.inc"
      dimension f(njp1,nk2p3)
C
#ifdef DEBUG_MODE_ON
      write(*,*) "DEBUG:  in subwts.F::fint(...)"
#endif
C
      fint = f(i,j)*e(l,m,1,1) + f(i+1,j)*e(l,m,2,1) +
     $     f(i,j+1)*e(l,m,3,1) + f(i+1,j+1)*e(l,m,4,1)
      return
      end
