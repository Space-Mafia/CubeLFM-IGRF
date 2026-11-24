      SUBROUTINE HIGH(F,FHF,LVAR)
#include "param.inc"
#include "mdims.inc"
#include "help.inc"
C
      DIMENSION F(llow:lIhigh,lJ),FHF(lIP1,lJ)
      COMMON /SCRCH1/
     $     FDUM(lIqq,lJ)
#ifdef DEBUG_MODE_ON
      write(*,*) "DEBUG:  in highs.F::HIGH(...)"
#endif
C
        DO 200 J=1,MJ
        DO 200 I=-ord2+1,MI+ord2
  200   FDUM(I+ORD2,J) = F(I,J)
C
*        CALL BRDHI(FDUM,LVAR)
C
        DO 299 J=1,MJ
        DO 299 I=1,MIP1
  299   FHF(I,J) = 0.
        DO 300 K=1,ORD2
        DO 300 J=1,MJ
        DO 300 I=1,MIP1
  300   FHF(I,J) = FHF(I,J) + FCOEF(K)*
     $          ( FDUM(I+K-1,J) + FDUM(I+ORDER-K,J) )
C 
C
*
*        if ( isweep .eq. 2 .and. jhigh_bound ) then
*        do j=1,mj
*        fhf(mip1,j) = 0.5*(fdum(mi+ord2,j)+fdum(mip1+ord2,j))
*        enddo
c
*
*        endif
C
        RETURN
        END
