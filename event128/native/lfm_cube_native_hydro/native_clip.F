      subroutine clip(fhf,f0)
!
! ... Local variables ..................................................
!
#include "param.inc"
#include "mdims.inc"
!
! ... Parameter variables ..............................................
!
      dimension fhf(LIp1,LJ)
!> Description of fhf
      dimension f0(LIp4,LJ)
!> Description of f0
!
! ... Begin ............................................................
!
#ifdef DEBUG_MODE_ON
      write(*,*) "DEBUG:  in clips.F::clip(fhf,f0)"
#endif
        do 200 j=1,MJ
        do 200 i=1,MIp1
  200   fhf(i,j) = amax1( amin1(f0(i+2,j),f0(i+1,j)),
     $          amin1( fhf(i,j), amax1(f0(i+2,j),f0(i+1,j))))
C
!
! ... End subroutine clip(fhf,f0) ......................................
!
      return
      end
