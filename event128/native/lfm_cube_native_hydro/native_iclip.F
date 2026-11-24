       subroutine Iclip(fhf,f0)
!
! ... Local variables ..................................................
!
#include "param.inc"
#include "mdims.inc"
!
! ... Parameter variables ..............................................
!
      dimension fhf(lIp1,lJP1)
!> description of fhf
      dimension f0(lIp4,lJP1)
!> description of f0
!
! ... Begin ............................................................
!
#ifdef DEBUG_MODE_ON
      write(*,*) "DEBUG:  in clips.F::Iclip(fhf,f0)"
#endif
        do 200 j=1,MJP1
        do 200 i=1,MIp1
  200   fhf(i,j) = amax1( amin1(f0(i+2,j),f0(i+1,j)),
     $          amin1( fhf(i,j), amax1(f0(i+2,j),f0(i+1,j))))
!
! ... End subroutine Iclip(fhf,f0) .....................................
!
        return
        end
