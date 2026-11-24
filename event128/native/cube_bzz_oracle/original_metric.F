      subroutine metric
#include "global_dims.inc"
#include "fortran_P++_arrays.inc"
#include "mhd_P++_pointers.inc"
#include "param.inc"
#include "help.inc"
#ifdef SPHERE
#include "meter.inc"
#endif
      dimension volrel(ni,nj,nk),fi_rel(ni,nj,nk),fj_rel(ni,nj,nk)
      parameter ( icent = ni+no/2, jcent=nj+no/2, kcent=nk+no/2)
      common /scrach/
     $   xcent_i(llow:icent,llow:jcent,llow:kcent),
     $   ycent_i(llow:icent,llow:jcent,llow:kcent),
     $   zcent_i(llow:icent,llow:jcent,llow:kcent),
     $   xcent_j(llow:icent,llow:jcent,llow:kcent),
     $   ycent_j(llow:icent,llow:jcent,llow:kcent),
     $   zcent_j(llow:icent,llow:jcent,llow:kcent),
     $   xcent_k(llow:icent,llow:jcent,llow:kcent),
     $   ycent_k(llow:icent,llow:jcent,llow:kcent),
     $   zcent_k(llow:icent,llow:jcent,llow:kcent)
*      external malloc,free
#ifdef PT8BYTE
      integer*8 fmalloc
#else
      integer*4 fmalloc
#endif
      data small/1.e-10/
      data big /1.e20/
c
c       temporary scratch space
c
***      pointer (pxc_i,xcent_i), (pyc_i,ycent_i), (pzc_i,zcent_i)
***      pointer (pxc_j,xcent_j), (pyc_j,ycent_j), (pzc_j,zcent_j)
***      pointer (pxc_k,xcent_k), (pyc_k,ycent_k), (pzc_k,zcent_k)
***      real*4 xcent_i(ilow:ihigh,jlow:jhigh,klow:khigh)
***      real*4 ycent_i(ilow:ihigh,jlow:jhigh,klow:khigh)
***      real*4 zcent_i(ilow:ihigh,jlow:jhigh,klow:khigh)
***      real*4 xcent_j(ilow:ihigh,jlow:jhigh,klow:khigh)
***      real*4 ycent_j(ilow:ihigh,jlow:jhigh,klow:khigh)
***      real*4 zcent_j(ilow:ihigh,jlow:jhigh,klow:khigh)
***      real*4 xcent_k(ilow:ihigh,jlow:jhigh,klow:khigh)
***      real*4 ycent_k(ilow:ihigh,jlow:jhigh,klow:khigh)
***      real*4 zcent_k(ilow:ihigh,jlow:jhigh,klow:khigh)
***c
***      length = 4*(ihigh-ilow+1)*(jhigh-jlow+1)*(khigh-klow+1)
***#ifdef T3E
***      pxc_i = fmalloc(length)
***      pyc_i = fmalloc(length)
***      pzc_i = fmalloc(length)
***      pxc_j = fmalloc(length)
***      pyc_j = fmalloc(length)
***      pzc_j = fmalloc(length)
***      pxc_k = fmalloc(length)
***      pyc_k = fmalloc(length)
***      pzc_k = fmalloc(length)
***#else
***      pxc_i = malloc(length)
***      pyc_i = malloc(length)
***      pzc_i = malloc(length)
***      pxc_j = malloc(length)
***      pyc_j = malloc(length)
***      pzc_j = malloc(length)
***      pxc_k = malloc(length)
***      pyc_k = malloc(length)
***      pzc_k = malloc(length)
***#endif
      integer i,j,k
#ifdef DEBUG_MODE_ON
      write(*,*) "DEBUG:  in metric.F::metric(...)"
#endif
      if ( .not. MHD_PE ) return
	do i = llow,ihigh
          do j = llow,jhigh
            do k = llow,khigh
              xcent_i(i,j,k) = 0.0
              ycent_i(i,j,k) = 0.0
              zcent_i(i,j,k) = 0.0
              xcent_j(i,j,k) = 0.0
              ycent_j(i,j,k) = 0.0
              zcent_j(i,j,k) = 0.0
              xcent_k(i,j,k) = 0.0
              ycent_k(i,j,k) = 0.0
              zcent_k(i,j,k) = 0.0
            enddo
          enddo
        enddo
c
c
c          i  direction stuff
c
c
#ifdef DEBUG
        write (6,*) 'in metric, before i face  ', mype
#endif
#ifdef RENORM
c  code to keep overflows from happening on small dynamic range machines
      a05 = 0.5e-10
      a1inv = 1.e20
#else
      a05 = 0.5
      a1inv =1.0
#endif
      do k=klow,khigh-1
         do j=jlow,jhigh-1
            do i=ilow,ihigh
               xcent_i(i,j,k) = 0.25*(x(i,j,k)+x(i,j+1,k)+
     $              x(i,j,k+1)+x(i,j+1,k+1))
               ycent_i(i,j,k) = 0.25*(y(i,j,k)+y(i,j+1,k)+
     $              y(i,j,k+1)+y(i,j+1,k+1))
               zcent_i(i,j,k) = 0.25*(z(i,j,k)+z(i,j+1,k)+
     $              z(i,j,k+1)+z(i,j+1,k+1))
               dx_0 = a05*(x(i,j+1,k)+x(i,j+1,k+1)-x(i,j,k)-x(i,j,k+1))
               dy_0 = a05*(y(i,j+1,k)+y(i,j+1,k+1)-y(i,j,k)-y(i,j,k+1))
               dz_0 = a05*(z(i,j+1,k)+z(i,j+1,k+1)-z(i,j,k)-z(i,j,k+1))
               dx_1 = a05*(x(i,j,k+1)+x(i,j+1,k+1)-x(i,j,k)-x(i,j+1,k))
               dy_1 = a05*(y(i,j,k+1)+y(i,j+1,k+1)-y(i,j,k)-y(i,j+1,k))
               dz_1 = a05*(z(i,j,k+1)+z(i,j+1,k+1)-z(i,j,k)-z(i,j+1,k))
               dx_a = dy_0*dz_1-dz_0*dy_1
               dy_a = dz_0*dx_1-dx_0*dz_1
               dz_a = dx_0*dy_1-dy_0*dx_1
               Face_I(i,j,k) = sqrt(dx_a**2+dy_a**2+dz_a**2)
               dinv = 1./max(Face_i(i,j,k),small)
               tran_i(i,j,k,1) = dinv*dx_a
               tran_i(i,j,k,2) = dinv*dy_a
               tran_i(i,j,k,3) = dinv*dz_a
               dinv = 1./max(small,sqrt(dx_0**2+dy_0**2+dz_0**2))
               tran_i(i,j,k,4) = dinv*dx_0
               tran_i(i,j,k,5) = dinv*dy_0
               tran_i(i,j,k,6) = dinv*dz_0
               dinv = 1./max(small,sqrt(dx_1**2+dy_1**2+dz_1**2))
               tran_i(i,j,k,7) = dinv*dx_1
               tran_i(i,j,k,8) = dinv*dy_1
               tran_i(i,j,k,9) = dinv*dz_1
               Face_I(i,j,k) = a1inv*Face_i(i,j,k)
            enddo
         enddo
       enddo
#ifdef AXIS_FIX
c  modifications for quasi-spherical grids
	if ( jlow_bound) then
	do ltran=1,9
	do k=klow,khigh-1
	do j=jlow,0
	do i=ilow,ihigh
	   tran_i(i,j,k,ltran) = - tran_i(i,j,k,ltran)
	enddo
	enddo
	enddo
	enddo
	endif
c
	if (jhigh_bound) then
	do ltran=1,9
	do k=klow,khigh-1
	do j=jhigh-no2,jhigh-1
	do i=ilow,ihigh
	   tran_i(i,j,k,ltran) = - tran_i(i,j,k,ltran)
	enddo
	enddo
	enddo
	enddo
	endif
	
#endif
#ifdef DEBUG
        write (6,*) 'in metric, before j face  ', mype
#endif
c
c        j-face stuff
c
      do k=klow,khigh-1
         do j=jlow,jhigh
            do i=ilow,ihigh-1
               xcent_j(i,j,k) = 0.25*(x(i,j,k)+x(i+1,j,k)+
     $              x(i,j,k+1)+x(i+1,j,k+1))
               ycent_j(i,j,k) = 0.25*(y(i,j,k)+y(i+1,j,k)+
     $              y(i,j,k+1)+y(i+1,j,k+1))
               zcent_j(i,j,k) = 0.25*(z(i,j,k)+z(i+1,j,k)+
     $              z(i,j,k+1)+z(i+1,j,k+1))
               dx_1 = a05*(x(i+1,j,k)+x(i+1,j,k+1)-x(i,j,k)-x(i,j,k+1))
               dy_1 = a05*(y(i+1,j,k)+y(i+1,j,k+1)-y(i,j,k)-y(i,j,k+1))
               dz_1 = a05*(z(i+1,j,k)+z(i+1,j,k+1)-z(i,j,k)-z(i,j,k+1))
               dx_0 = a05*(x(i,j,k+1)+x(i+1,j,k+1)-x(i,j,k)-x(i+1,j,k))
               dy_0 = a05*(y(i,j,k+1)+y(i+1,j,k+1)-y(i,j,k)-y(i+1,j,k))
               dz_0 = a05*(z(i,j,k+1)+z(i+1,j,k+1)-z(i,j,k)-z(i+1,j,k))
               dx_a = dy_0*dz_1-dz_0*dy_1
               dy_a = dz_0*dx_1-dx_0*dz_1
               dz_a = dx_0*dy_1-dy_0*dx_1
               Face_j(i,j,k) = sqrt(dx_a**2+dy_a**2+dz_a**2)
               dinv = 1.0/max(small,Face_j(i,j,k))
               tran_j(i,j,k,1) = dinv*dx_a
               tran_j(i,j,k,2) = dinv*dy_a
               tran_j(i,j,k,3) = dinv*dz_a
               dinv = 1.0/max(small,sqrt(dx_0**2+dy_0**2+dz_0**2))
               tran_j(i,j,k,4) = dinv*dx_0
               tran_j(i,j,k,5) = dinv*dy_0
               tran_j(i,j,k,6) = dinv*dz_0
               dinv = 1.0/max(small,sqrt(dx_1**2+dy_1**2+dz_1**2))
               tran_j(i,j,k,7) = dinv*dx_1
               tran_j(i,j,k,8) = dinv*dy_1
               tran_j(i,j,k,9) = dinv*dz_1
               Face_J(i,j,k) = a1inv*Face_J(i,j,k)
            enddo
         enddo
       enddo
c
#ifdef AXIS_FIX
c  modifications for quasi-spherical grids
	if ( jlow_bound) then
	do ltran=1,9
	do k=klow,khigh-1
	do j=jlow,0
	do i=ilow,ihigh-1
	   tran_j(i,j,k,ltran) = - tran_j(i,j,k,ltran)
	enddo
	enddo
	enddo
	enddo
	endif
c
	if (jhigh_bound) then
	do ltran=1,9
	do k=klow,khigh-1
	do j=jhigh-no2+1,jhigh
	do i=ilow,ihigh-1
	   tran_j(i,j,k,ltran) = - tran_j(i,j,k,ltran)
	enddo
	enddo
	enddo
	enddo
	endif
	
#endif
c
c
c       k-face stuff
c
c
        iunit = mype+20
#ifdef DEBUG
        write (6,*) 'in metric, before k face  ', mype
        call forflush
        if (mype .eq. 3) then
        kup = khigh-klow+1
        jup = jhigh-jlow+1
        iup = ihigh-ilow+1
        write (6,*) 'ranges ',ilow,ihigh,jlow,jhigh,klow,khigh
        call pictur(x,'XXX',1,iup,1,jup,1,kup,0,1,2,3,
     $       iup,jup,kup,iup*jup)
        call pictur(y,'YYY',1,iup,1,jup,1,kup,0,1,2,3,
     $       iup,jup,kup,iup*jup)
        call pictur(z,'ZZZ',1,iup,1,jup,1,kup,0,1,2,3,
     $       iup,jup,kup,iup*jup)
        endif
#endif
      do k=klow,khigh
         do j=jlow,jhigh-1
            do i=ilow,ihigh-1
#ifdef DEBUG
               write (iunit,*) '>>>>>>>>>>>>>>>>  ',i,j,k
               call forflush
#endif
               xcent_k(i,j,k) = 0.25*(x(i,j,k)+x(i,j+1,k)+
     $              x(i+1,j,k)+x(i+1,j+1,k))
               ycent_k(i,j,k) = 0.25*(y(i,j,k)+y(i,j+1,k)+
     $              y(i+1,j,k)+y(i+1,j+1,k))
               zcent_k(i,j,k) = 0.25*(z(i,j,k)+z(i,j+1,k)+
     $              z(i+1,j,k)+z(i+1,j+1,k))
               dx_1 = a05*(x(i,j+1,k)+x(i+1,j+1,k)-x(i,j,k)-x(i+1,j,k))
               dy_1 = a05*(y(i,j+1,k)+y(i+1,j+1,k)-y(i,j,k)-y(i+1,j,k))
               dz_1 = a05*(z(i,j+1,k)+z(i+1,j+1,k)-z(i,j,k)-z(i+1,j,k))
               dx_0 = a05*(x(i+1,j,k)+x(i+1,j+1,k)-x(i,j,k)-x(i,j+1,k))
               dy_0 = a05*(y(i+1,j,k)+y(i+1,j+1,k)-y(i,j,k)-y(i,j+1,k))
               dz_0 = a05*(z(i+1,j,k)+z(i+1,j+1,k)-z(i,j,k)-z(i,j+1,k))
               dx_a = dy_0*dz_1-dz_0*dy_1
               dy_a = dz_0*dx_1-dx_0*dz_1
               dz_a = dx_0*dy_1-dy_0*dx_1
#ifdef DEBUG
               write(iunit,*) xcent_k(i,j,k),ycent_k(i,j,k),
     $              zcent_k(i,j,k)
               call forflush
               write(iunit,*) dx_1,dy_1,dz_1,dx_0,dy_0,dz_0
               call forflush
#endif
               Face_k(i,j,k) = sqrt(dx_a**2+dy_a**2+dz_a**2)
               dinv = 1.0/max(small,Face_k(i,j,k))
               tran_k(i,j,k,1) = dinv*dx_a
               tran_k(i,j,k,2) = dinv*dy_a
               tran_k(i,j,k,3) = dinv*dz_a
               dinv = 1.0/max(small,sqrt(dx_0**2+dy_0**2+dz_0**2))
               tran_k(i,j,k,4) = dinv*dx_0
               tran_k(i,j,k,5) = dinv*dy_0
               tran_k(i,j,k,6) = dinv*dz_0
               dinv = 1.0/max(small,sqrt(dx_1**2+dy_1**2+dz_1**2))
               tran_k(i,j,k,7) = dinv*dx_1
               tran_k(i,j,k,8) = dinv*dy_1
               tran_k(i,j,k,9) = dinv*dz_1
               Face_K(i,j,k) = a1inv*Face_K(i,j,k)
            enddo
         enddo
       enddo
c
c      cell volumes
c
#ifdef DEBUG
        write (6,*) 'in metric, before volumes  ', mype
#endif
       do k=klow,khigh-1
          do j=jlow,jhigh-1
             do i=ilow,ihigh-1
                dx_i = xcent_i(i+1,j,k)-xcent_i(i,j,k)
                dy_i = ycent_i(i+1,j,k)-ycent_i(i,j,k)
                dz_i = zcent_i(i+1,j,k)-zcent_i(i,j,k)
                dx_j = xcent_j(i,j+1,k)-xcent_j(i,j,k)
                dy_j = ycent_j(i,j+1,k)-ycent_j(i,j,k)
                dz_j = zcent_j(i,j+1,k)-zcent_j(i,j,k)
                dx_k = xcent_k(i,j,k+1)-xcent_k(i,j,k)
                dy_k = ycent_k(i,j,k+1)-ycent_k(i,j,k)
                dz_k = zcent_k(i,j,k+1)-zcent_k(i,j,k)
                volumetmp = abs(dx_i*(dy_j*dz_k-dy_k*dz_j) +
     $                          dy_i*(dz_j*dx_k-dx_j*dz_k) +
     $                          dz_i*(dx_j*dy_k-dy_j*dx_k))
                volume(i,j,k) = abs(dx_i*(dy_j*dz_k-dy_k*dz_j) +
     $                          dy_i*(dz_j*dx_k-dx_j*dz_k) +
     $                          dz_i*(dx_j*dy_k-dy_j*dx_k))
                volinv(i,j,k) = 1./max(small,volume(i,j,k))
                if ( i .eq. 10 .and. j .eq. 10 ) then
                   klev00 = k
                endif
             enddo
          enddo
       enddo
c
c
c      edge stuff for magnetic field push
c
c      (7-9) vector in edge direction (Z, say)
c      (1-3) vector in direction of sweep (x, say)
c      (4-6) vector right handed perp to both
c
c         
c                         k -edge
c
#ifdef DEBUG
        write (6,*) 'in metric, before k edge  ', mype
#endif
       do k=klow,khigh-1
          do j=jlow+1,jhigh-1
             do i=ilow+1,ihigh-1
                dx_e = x(i,j,k+1)-x(i,j,k)
                dy_e = y(i,j,k+1)-y(i,j,k)
                dz_e = z(i,j,k+1)-z(i,j,k)
                edge_k(i,j,k) = sqrt(dx_e**2+dy_e**2+dz_e**2)
                dinv = 1./max(small,edge_k(i,j,k))
                dx_e = dx_e*dinv
                dy_e = dy_e*dinv
                dz_e = dz_e*dinv
                vec_k(i,j,k,7) = dx_e
                vec_k(i,j,k,8) = dy_e
                vec_k(i,j,k,9) = dz_e
                dxd = 0.25*(x(i,j+1,k+1)+x(i,j+1,k)-
     $                           x(i,j-1,k+1)-x(i,j-1,k))
                dyd = 0.25*(y(i,j+1,k+1)+y(i,j+1,k)-
     $                           y(i,j-1,k+1)-y(i,j-1,k))
                dzd = 0.25*(z(i,j+1,k+1)+z(i,j+1,k)-
     $                           z(i,j-1,k+1)-z(i,j-1,k))
                dx_0 = dyd*dz_e-dzd*dy_e
                dy_0 = dzd*dx_e-dxd*dz_e
                dz_0 = dxd*dy_e-dyd*dx_e
                dinv = 1./max(small,sqrt(dx_0**2+dy_0**2+dz_0**2))
                vec_k(i,j,k,1) = dx_0*dinv
                vec_k(i,j,k,2) = dy_0*dinv
                vec_k(i,j,k,3) = dz_0*dinv
                dx_1 = -dy_0*dz_e+dz_0*dy_e
                dy_1 = -dz_0*dx_e+dx_0*dz_e
                dz_1 = -dx_0*dy_e+dy_0*dx_e
                dinv = 1./max(small,sqrt(dx_1**2+dy_1**2+dz_1**2))
                vec_k(i,j,k,4) = dx_1*dinv
                vec_k(i,j,k,5) = dy_1*dinv
                vec_k(i,j,k,6) = dz_1*dinv
             enddo
          enddo
       enddo
c
#ifdef DEBUG
        write (6,*) 'in metric, before i edge  ', mype
#endif
c         
c                         i -edge
c
       do k=klow+1,khigh-1
          do j=jlow+1,jhigh-1
             do i=ilow,ihigh-1
                dx_e = x(i+1,j,k)-x(i,j,k)
                dy_e = y(i+1,j,k)-y(i,j,k)
                dz_e = z(i+1,j,k)-z(i,j,k)
                edge_i(i,j,k) = sqrt(dx_e**2+dy_e**2+dz_e**2)
                dinv = 1./max(small,edge_i(i,j,k))
                dx_e = dx_e*dinv
                dy_e = dy_e*dinv
                dz_e = dz_e*dinv
                vec_i(i,j,k,7) = dx_e
                vec_i(i,j,k,8) = dy_e
                vec_i(i,j,k,9) = dz_e
                dxd = 0.25*(x(i+1,j,k+1)+x(i,j,k+1)-
     $                           x(i+1,j,k-1)-x(i,j,k-1))
                dyd = 0.25*(y(i+1,j,k+1)+y(i,j,k+1)-
     $                           y(i+1,j,k-1)-y(i,j,k-1))
                dzd = 0.25*(z(i+1,j,k+1)+z(i,j,k+1)-
     $                           z(i+1,j,k-1)-z(i,j,k-1))
                dx_0 = dyd*dz_e-dzd*dy_e
                dy_0 = dzd*dx_e-dxd*dz_e
                dz_0 = dxd*dy_e-dyd*dx_e  
                dinv = 1./max(small,sqrt(dx_0**2+dy_0**2+dz_0**2))
                vec_i(i,j,k,1) = dx_0*dinv
                vec_i(i,j,k,2) = dy_0*dinv
                vec_i(i,j,k,3) = dz_0*dinv
                dx_1 = -dy_0*dz_e+dz_0*dy_e
                dy_1 = -dz_0*dx_e+dx_0*dz_e
                dz_1 = -dx_0*dy_e+dy_0*dx_e
                dinv = 1./max(small,sqrt(dx_1**2+dy_1**2+dz_1**2))
                vec_i(i,j,k,4) = dx_1*dinv
                vec_i(i,j,k,5) = dy_1*dinv
                vec_i(i,j,k,6) = dz_1*dinv
             enddo
          enddo
       enddo
c
c         
c                         j -edge
c
#ifdef DEBUG
        write (6,*) 'in metric, before j edge  ', mype
#endif
       do k=klow+1,khigh-1
          do j=jlow+1,jhigh-1
             do i=ilow+1,ihigh-1
                dx_e = x(i,j+1,k)-x(i,j,k)
                dy_e = y(i,j+1,k)-y(i,j,k)
                dz_e = z(i,j+1,k)-z(i,j,k)
                edge_j(i,j,k) = sqrt(dx_e**2+dy_e**2+dz_e**2)
                dinv = 1./max(small,edge_j(i,j,k))
                dx_e = dx_e*dinv
                dy_e = dy_e*dinv
                dz_e = dz_e*dinv
                vec_j(i,j,k,7) = dx_e
                vec_j(i,j,k,8) = dy_e
                vec_j(i,j,k,9) = dz_e
                dxd = 0.25*(x(i+1,j+1,k)+x(i+1,j,k)-
     $                           x(i-1,j+1,k)-x(i-1,j,k))
                dyd = 0.25*(y(i+1,j+1,k)+y(i+1,j,k)-
     $                           y(i-1,j+1,k)-y(i-1,j,k))
                dzd = 0.25*(z(i+1,j+1,k)+z(i+1,j,k)-
     $                           z(i-1,j+1,k)-z(i-1,j,k))
                dx_0 = dyd*dz_e-dzd*dy_e
                dy_0 = dzd*dx_e-dxd*dz_e
                dz_0 = dxd*dy_e-dyd*dx_e
                dinv = 1./max(small,sqrt(dx_0**2+dy_0**2+dz_0**2))
                vec_j(i,j,k,1) = dx_0*dinv
                vec_j(i,j,k,2) = dy_0*dinv
                vec_j(i,j,k,3) = dz_0*dinv
                dx_1 = -dy_0*dz_e+dz_0*dy_e
                dy_1 = -dz_0*dx_e+dx_0*dz_e
                dz_1 = -dx_0*dy_e+dy_0*dx_e
                dinv = 1./max(small,sqrt(dx_1**2+dy_1**2+dz_1**2))
                vec_j(i,j,k,4) = dx_1*dinv
                vec_j(i,j,k,5) = dy_1*dinv
                vec_j(i,j,k,6) = dz_1*dinv
             enddo
          enddo
       enddo
c
c     for courant calculation
c
       do k=klow,khigh
          do j=jlow,jhigh
             do i=ilow,ihigh
                cell_length(i,j,k) = BIG
             enddo
          enddo
       enddo
       do k=1,nl_k
          do j=1,nl_j
             avit = 1.0
#ifdef SPHERE
             jglob = j + j_global_off
             if ( jglob .gt. nj_global/2 ) jglob = nj_global+1 - jglob
             if ( jglob .le. nrings) then
               avit = 1.0/naver(jglob)
             endif
#endif
             do i=1,nl_i
                cell_length(i,j,k) = volume(i,j,k)/
     $               max( face_i(i,j,k),face_i(i+1,j,k),
     $                    face_j(i,j,k),face_j(i,j+1,k),
     $                    avit*face_k(i,j,k))
             enddo
          enddo
       enddo
c
c
c
c
c   free up temps
c
***       call free(pxc_i)
***       call free(pyc_i)
***       call free(pzc_i)
***       call free(pxc_j)
***       call free(pyc_j)
***       call free(pzc_j)
***       call free(pxc_k)
***       call free(pyc_k)
***       call free(pzc_k)
c
c
c   legacy stuff to make old code run right
c
#ifdef DEBUG
        write (6,*) 'in metric, before old 2d  ', mype
#endif
#ifdef SPHERE
       pi = 4.0*atan(1.0)
      do k=1,nl_k+1
         phi(k) = atan2(z(2,2,k),y(2,2,k))
         if (phi(k) .lt. 0.0 ) phi(k) = phi(k) + 2.*pi
      enddo
         phi(nl_k+1) = phi(1) + 2.*pi
      do k=1,nl_k+1
        COSPHI(K) = COS(PHI(K))
        SINPHI(K) = SIN(PHI(K))
      enddo
      DO 500 K=1,nl_k
        CPHIqM(K) = COS(0.5*(PHI(K)+PHI(K+1)))
        SPHIqM(K) = SIN(0.5*(PHI(K)+PHI(K+1)))
        DARC(K) = 2.0*SIN(0.5*(PHI(K+1)-PHI(K)))
        COSARC(K) = COS( 0.5*(PHI(K+1)-PHI(K)))
  500   CONTINUE
#endif
C
#ifdef METRICBUG
        iunit = mype + 20
        write (iunit,*) 'metric stuff'
 912    format (8(1pe12.4))
        write (iunit,*) 'face_i'
        write (iunit,912) (face_i(5,10,kf),kf=1,nl_k)
        write (iunit,*) 'face_j'
        write (iunit,912) (face_j(5,10,kf),kf=1,nl_k)
        write (iunit,*) 'face_k'
        write (iunit,912) (face_k(5,10,kf),kf=1,nl_k)
        write (iunit,*) 'edge_i'
        write (iunit,912) (edge_i(5,10,kf),kf=1,nl_k)
        write (iunit,*) 'edge_j'
        write (iunit,912) (edge_j(5,10,kf),kf=1,nl_k)
        write (iunit,*) 'edge_k'
        write (iunit,912) (edge_k(5,10,kf),kf=1,nl_k)
        do klev = 1,nl_k
        write (iunit,*) '>>>>>>>>>>>>>>>> K value',klev
 911    format (i3,2x,9(1pe11.3))
        write (iunit,*) 'tran_i and vec_i' 
        dettran = tran_i(5,10,klev,1)*(
     $       tran_i(5,10,klev,5)*tran_i(5,10,klev,9) -
     $         tran_i(5,10,klev,6)*tran_i(5,10,klev,8)) +
     $     tran_i(5,10,klev,2)*(
     $       tran_i(5,10,klev,6)*tran_i(5,10,klev,7) -
     $       tran_i(5,10,klev,4)*tran_i(5,10,klev,9))  +
     $     tran_i(5,10,klev,3)*(
     $       tran_i(5,10,klev,4)*tran_i(5,10,klev,8)-
     $       tran_i(5,10,klev,5)*tran_i(5,10,klev,7))
        detvec = vec_i(5,10,klev,1)*(
     $       vec_i(5,10,klev,5)*vec_i(5,10,klev,9) -
     $         vec_i(5,10,klev,6)*vec_i(5,10,klev,8)) +
     $     vec_i(5,10,klev,2)*(
     $       vec_i(5,10,klev,6)*vec_i(5,10,klev,7) -
     $       vec_i(5,10,klev,4)*vec_i(5,10,klev,9))  +
     $     vec_i(5,10,klev,3)*(
     $       vec_i(5,10,klev,4)*vec_i(5,10,klev,8)-
     $       vec_i(5,10,klev,5)*vec_i(5,10,klev,7))
        write (iunit,*) klev,dettran,detvec
           
        write (iunit,911) klev,(tran_i(5,10,klev,l),l=1,9)
        write (iunit,911) klev,(vec_i(5,10,klev,l),l=1,9)
        write (iunit,*) 'tran_j and vec_j'
        dettran = tran_j(5,10,klev,1)*(
     $       tran_j(5,10,klev,5)*tran_j(5,10,klev,9) -
     $         tran_j(5,10,klev,6)*tran_j(5,10,klev,8)) +
     $     tran_j(5,10,klev,2)*(
     $       tran_j(5,10,klev,6)*tran_j(5,10,klev,7) -
     $       tran_j(5,10,klev,4)*tran_j(5,10,klev,9))  +
     $     tran_j(5,10,klev,3)*(
     $       tran_j(5,10,klev,4)*tran_j(5,10,klev,8)-
     $       tran_j(5,10,klev,5)*tran_j(5,10,klev,7))
        detvec = vec_j(5,10,klev,1)*(
     $       vec_j(5,10,klev,5)*vec_j(5,10,klev,9) -
     $         vec_j(5,10,klev,6)*vec_j(5,10,klev,8)) +
     $     vec_j(5,10,klev,2)*(
     $       vec_j(5,10,klev,6)*vec_j(5,10,klev,7) -
     $       vec_j(5,10,klev,4)*vec_j(5,10,klev,9))  +
     $     vec_j(5,10,klev,3)*(
     $       vec_j(5,10,klev,4)*vec_j(5,10,klev,8)-
     $       vec_j(5,10,klev,5)*vec_j(5,10,klev,7))
        write (iunit,*) klev,dettran,detvec
        write (iunit,911) klev,(tran_j(5,10,klev,l),l=1,9)
        write (iunit,911) klev,(vec_j(5,10,klev,l),l=1,9)
        write (iunit,*) 'tran_k and vec_k'
        dettran = tran_k(5,10,klev,1)*(
     $       tran_k(5,10,klev,5)*tran_k(5,10,klev,9) -
     $         tran_k(5,10,klev,6)*tran_k(5,10,klev,8)) +
     $     tran_k(5,10,klev,2)*(
     $       tran_k(5,10,klev,6)*tran_k(5,10,klev,7) -
     $       tran_k(5,10,klev,4)*tran_k(5,10,klev,9))  +
     $     tran_k(5,10,klev,3)*(
     $       tran_k(5,10,klev,4)*tran_k(5,10,klev,8)-
     $       tran_k(5,10,klev,5)*tran_k(5,10,klev,7))
        detvec = vec_k(5,10,klev,1)*(
     $       vec_k(5,10,klev,5)*vec_k(5,10,klev,9) -
     $         vec_k(5,10,klev,6)*vec_k(5,10,klev,8)) +
     $     vec_k(5,10,klev,2)*(
     $       vec_k(5,10,klev,6)*vec_k(5,10,klev,7) -
     $       vec_k(5,10,klev,4)*vec_k(5,10,klev,9))  +
     $     vec_k(5,10,klev,3)*(
     $       vec_k(5,10,klev,4)*vec_k(5,10,klev,8)-
     $       vec_k(5,10,klev,5)*vec_k(5,10,klev,7))
        write (iunit,*) klev,dettran,detvec
        write (iunit,911) klev,(tran_k(5,10,klev,l),l=1,9)
        write (iunit,911) klev,(vec_k(5,10,klev,l),l=1,9)
        write (iunit,*) 'pressure test'
        fxi = -tran_i(5,10,klev,1)*face_i(5,10,klev) +
     $       tran_i(6,10,klev,1)*face_i(6,10,klev)
        fxj = -tran_j(5,10,klev,1)*face_j(5,10,klev) +
     $       tran_j(5,11,klev,1)*face_j(5,11,klev) 
        fxk =   -tran_k(5,10,klev,1)*face_k(5,10,klev) +
     $       tran_k(5,10,klev+1,1)*face_k(5,10,klev+1)
        fxx = fxi + fxj + fxk
        fyi = -tran_i(5,10,klev,2)*face_i(5,10,klev) +
     $       tran_i(6,10,klev,2)*face_i(6,10,klev)
        fyj =   -tran_j(5,10,klev,2)*face_j(5,10,klev) +
     $       tran_j(5,11,klev,2)*face_j(5,11,klev) 
        fyk =   - tran_k(5,10,klev,2)*face_k(5,10,klev) +
     $       tran_k(5,10,klev+1,2)*face_k(5,10,klev+1)
        fyy = fyi + fyj + fyk
        fzi = -tran_i(5,10,klev,3)*face_i(5,10,klev) +
     $       tran_i(6,10,klev,3)*face_i(6,10,klev)
        fzj =  -tran_j(5,10,klev,3)*face_j(5,10,klev) +
     $       tran_j(5,11,klev,3)*face_j(5,11,klev) 
        fzk = -tran_k(5,10,klev,3)*face_k(5,10,klev) +
     $       tran_k(5,10,klev+1,3)*face_k(5,10,klev+1)
        fzz = fzi + fzj + fzk
        write (iunit,*) fxx,fyy,fzz
        write (iunit,*) fxx,fxi,fxj,fxk
        write (iunit,*) fyy,fyi,fyj,fyk
        write (iunit,*) fzz,fzi,fzj,fzk
        enddo
#endif
       return
       end
