        subroutine ave2flux
! ... Local variables ..................................................
#include "global_dims.inc"
#include "mhd_P++_pointers.inc"
#include "fortran_P++_arrays.inc"
#include "param.inc"
! ... Begin ............................................................
#ifdef DEBUG_MODE_ON
      write(*,*) "DEBUG:  in therest.F::ave2flux()"
#endif
C
      if ( .not. MHD_PE ) return
C 
        DO 340 K=1,NL_K
        DO 340 J=1,NL_J
        DO 340 I=1,NL_I+1
  340   BI_g(I,J,K) = BI_g(I,J,K)*Face_I(I,J,K)
C
        DO 350 K=1,NL_K
        DO 350 J=1,NL_J+1
        DO 350 I=1,NL_I
  350   BJ_g(I,J,K) = BJ_g(I,J,K)*Face_J(I,J,K)
C
        DO 360 K=1,NL_K+1
        DO 360 J=1,NL_J
        DO 360 I=1,NL_I
  360   BK_g(I,J,K) = BK_g(I,J,K)*Face_K(I,J,K)
!
! ... End subroutine ave2flux ...........................................
        return
        end
