!> curinterp.F
!! long description of what's here
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> cur_coefs
!!
      subroutine cur_coefs
!
! ... Local variables ..................................................
!
#include "ion-param.inc"
#include "ion_P++_pointers.inc"
#include "ion_fortran_P++_arrays.inc"
#include "iongeo.inc"
#include "curint.inc"
      dimension dy2(njp1)
      dimension y2norm(njp1,nk2p3)
      dimension summit(njp1,nk2p3)
      logical stop_in_int
!
! ... Begin ............................................................
!
#ifdef DEBUG_MODE_ON
      write(*,*) "DEBUG:  in curinterp.F::cur_coefs()"
#endif
      do j=2,nj
         dy2(j) = (y2ion(j,nk2p3)-y2ion(j,1))/float(nk2p2)
      enddo
c
      do k=1,nk2p3
         do j=2,nj
            y2norm(j,k) = (y2ion_old(j,k)-y2ion(j,1))/dy2(j) + 1.
         enddo
      enddo
c
c   check the number of cells over which interpolation is done
c
      stop_in_int = .false.
      do j=2,nj
         do k=2,nk2p2
            lmin(j,k) = nk2p3
            lmax(j,k) = 1
            do lint = 2,nk2p2
               if ( y2norm(j,lint) .le. k+1 ) lmax(j,k)=lint+1
               if ( y2norm(j,nk2p3-lint+1) .ge. k-1 )
     $              lmin(j,k)=nk2p3 - lint 
            enddo
           lcint(j,k) = lmax(j,k)-lmin(j,k)+1
           if (lcint(j,k) .gt. ncint) then
              write (6,*) 'interpolation covers too many original cells'
              write (6,*) ' recompile with a bigger ncint '
              write (6,*) 'j,k,lmax,lmin'
              stop_in _int = .true.
           endif
         enddo
      enddo
c
c  get the coefficients: we're doing an integration using chapeau functions
c
c
c  first upward slant
c
      do j=2,nj
         do k=2,nk2p2
            do lint = 1,nk2p2
             if (       .not.
     $        (( y2norm(j,lint) .ge. k .and. y2norm(j,lint+1) .gt. k )
     $                              .or.
     $        ((y2norm(j,lint) .lt. k-1 .and. y2norm(j,lint+1) .le. k-1)
     $              )))   then
                  dyinv = 1./(y2norm(j,lint+1)-y2norm(j,lint))
                  llow = lint - lmin(j,k) + 1
                  lhi = llow + 1
                  ynplus = y2norm(j,lint+1)-float(k-1)
                  ynminus = y2norm(j,lint) -float(k-1)
                  ynlow = max(float(k-1),y2norm(j,lint))-float(k-1)
                  ynhi = min(float(k),y2norm(j,lint+1))-float(k-1)
                  curcoef(j,k,lhi) = dyinv*(
     $                 (ynhi**3-ynlow**3)/3.0 -
     $                 0.5*ynminus*(ynhi**2-ynlow**2))
     $                 +curcoef(j,k,lhi)
                  curcoef(j,k,llow) = dyinv*(
     $                 0.5*ynplus*(ynhi**2-ynlow**2) -
     $                 (ynhi**3-ynlow**3)/3.0) +curcoef(j,k,llow)
              endif
            enddo
          enddo
       enddo

c
c  second downward slant
c
      do j=2,nj
         do k=2,nk2p2
            do lint = 2,nk2p2
             if (       .not.
     $        (( y2norm(j,lint) .ge. k+1
     $              .and. y2norm(j,lint+1) .gt. k+1 )
     $                              .or.
     $        (y2norm(j,lint) .lt. k .and. y2norm(j,lint+1) .le. k)
     $              ))   then
                  dyinv = 1./(y2norm(j,lint+1)-y2norm(j,lint))
                  llow = lint - lmin(j,k) + 1
                  lhi = llow + 1
                  ynplus = y2norm(j,lint+1)-float(k)
                  ynminus = y2norm(j,lint) -float(k)
                  ynlow = max(float(k),y2norm(j,lint))-float(k)
                  ynhi = min(float(k+1),y2norm(j,lint+1))-float(k)
                  curcoef(j,k,llow) = dyinv*(
     $                 ynplus*(ynhi-ynlow) -
     $                 (1.+ynplus)*0.5*(ynhi**2-ynlow**2)
     $                 +(ynhi**3-ynlow**3)/3.0)
     $                 +curcoef(j,k,llow)
                  curcoef(j,k,lhi) = dyinv*(
     $                 -ynminus*(ynhi-ynlow)+
     $                 (ynminus+1.)*0.5*(ynhi**2-ynlow**2) -
     $                 (ynhi**3-ynlow**3)/3.0) +curcoef(j,k,lhi)
              endif
            enddo
          enddo
       enddo
c
!jgl  debug code
       do k=2,nk2p2
          do j=2,nj
             summit(j,k) = 0.
             do lint=1,ncint
                summit(j,k) = summit(j,k) + curcoef(j,k,lint)
             enddo
          enddo
       enddo
c
       summax = summit(2,2)
       summin = summit(2,2)
       do k=2,nk2p2
       do j=2,nj
          smaxold = summax
          summax = max(summax,summit(j,k))
          if ( summax .gt. smaxold ) then
             jsmax = j
             ksmax = k
          endif
          sminold = summin
          summin = min(summin,summit(j,k))
          if ( summin .gt. sminold ) then
             jsmin = j
             ksmin = k
          endif
       enddo
       enddo
#ifdef DEBUG_MODE_ON
       write (6,*) 'Curinterp info - max: ',summax,jsmax,ksmax
       write (6,*) 'Curinterp info - min: ',summin,jsmin,ksmin
#endif
!
! ... End subroutine cur_coefs .........................................
!
       return
       end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> flop_y
!! Long description
      subroutine flop_y
!
! ... Local variables ..................................................
!
#include "ion-param.inc"
#include "ion_P++_pointers.inc"
#include "ion_fortran_P++_arrays.inc"
#include "iongeo.inc"
#include "curint.inc"
#include "global_dims.inc"
!
! ... Begin ............................................................
!
#ifdef DEBUG_MODE_ON
      write(*,*) "DEBUG:  in curinterp.F::flop_y()"
#endif
      do k=1,nk2p3
         do j=1,njp1
            xinterp(j,k) = x2ion(j,k)
            yinterp(j,k) = y2ion(j,k)
         enddo
      enddo
c
#ifdef Y2INTERP
      do k=1,nk2p3
         do j=1,njp1
            y2ion(j,k) = y2ion_old(j,k)
         enddo
      enddo
#endif
!
! ... End subroutine flop_y ............................................
      return
      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> unflop_y
!!
      subroutine unflop_y
!
! ... Local variables ..................................................
!
#include "ion-param.inc"
#include "ion_P++_pointers.inc"
#include "ion_fortran_P++_arrays.inc"
#include "iongeo.inc"
#include "curint.inc"
#include "global_dims.inc"
!
! ... Begin ............................................................
!
#ifdef DEBUG_MODE_ON
      write(*,*) "DEBUG:  in curinterp.F::unflop_y()"
#endif

#ifdef Y2INTERP
      do k=1,nk2p3
         do j=1,njp1
            y2ion(j,k) = yinterp(j,k)
         enddo
      enddo
#endif
! ... End subroutine unflop_y...........................................
      return
      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
