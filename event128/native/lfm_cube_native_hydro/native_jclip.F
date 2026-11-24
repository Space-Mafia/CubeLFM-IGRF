      subroutine jclip(fhf,f0)
!
! ... Local variables ..................................................
!
#include "param.inc"
#include "mdims.inc"
!
! ... Parameter variables ..............................................
!
      dimension fhf(lIP1,lJp1)
!> description of fhf
      dimension f0(lIP1,lJp4)
!> description of f0
!
! ... Begin ............................................................
!
#ifdef DEBUG_MODE_ON
      write(*,*) "DEBUG:  in clips.F::jclip(fhf,f0)"
#endif
        do 200 j=1,MJp1
        do 200 i=1,MIP1
  200   fhf(i,j) = amax1( amin1(f0(i,j+2),f0(i,j+1)),
     $          amin1( fhf(i,j), amax1(f0(i,j+2),f0(i,j+1))))
!
! ... End subroutine jclip(fhf,f0) .....................................
!
        return
        end
