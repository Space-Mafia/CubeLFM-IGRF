        function bxqq0(x,y,z)
!
! ... Local variables ..................................................
!
#include "dipole.inc"
!
! ... Begin ............................................................
!
#ifdef DEBUG_MODE_ON
      write(*,*) "DEBUG:  in dipole.F::bxqq0(x,y,z)"
#endif
        rsq = x**2 + y**2 + z**2
        r = sqrt(rsq)
        if ( r .le. rinner ) then
          bxqq0 = -(geoqmu/rsq)*(3.*x*z/rsq)/r
        else if ( r .ge. router ) then
          bxqq0 = 0.
        else
          rprim = amax1(rinner,amin1(router,r))
          bxqq0 = -(geoqmu/rsq)*(3.*x*z/rsq)/r*(a1 + b1*rprim)
        endif
!
! ... End function bxqq0(x,y,z).........................................
!
        return
        end
