      function bzqq0(x,y,z)
!
! ... Local variables...................................................
!
#include "dipole.inc"
! 
! ... Begin ............................................................
!
#ifdef DEBUG_MODE_ON
      write(*,*) "DEBUG:  in dipole.F::bzqq0(x,y,z)"
#endif

      rsq = x**2 + y**2 + z**2
      r = sqrt(rsq)
      if ( r .le. rinner ) then
         bzqq0 = -(geoqmu/rsq)*((2. - 3.*(x**2+y**2)/rsq)/r)
      else if ( r .ge. router ) then
         bzqq0 = 0.
      else
         rprim = amax1(rinner,amin1(router,r))
         bzqq0 = -(geoqmu/rsq)*(2.-3.*(x**2+y**2)/rsq)
     $        /(r)*(a1 + b1*rprim)
      endif
!
! ... End function bzqq0(x,y,z).........................................
!   
      return
      end
