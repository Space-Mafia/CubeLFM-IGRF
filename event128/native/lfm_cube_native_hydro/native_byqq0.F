        function byqq0(x,y,z)
!
! ... Local Variables ..................................................
!
#include "dipole.inc"
!
! ... Begin ............................................................
!
#ifdef DEBUG_MODE_ON
      write(*,*) "DEBUG:  in dipole.F::byqq0(x,y,z)"
#endif
        rsq = x**2 + y**2 + z**2
        r = sqrt(rsq)
        if ( r .le. rinner ) then
          byqq0 = -(geoqmu/rsq)*(3.*y*z/rsq)/r
        else if ( r .ge. router ) then
          byqq0 = 0.
        else
          rprim = amax1(rinner,amin1(router,r))
          byqq0 = -(geoqmu/rsq)*(3.*y*z/rsq)/r*(a1 + b1*rprim)
        endif
!
! ... End function byqq(x,y,z) .........................................
!
        return
        end
