        SUBROUTINE YBHI(F,FHF,LVAR)
#include "param.inc"
#include "global_dims.inc"
#include "mdims.inc"
#include "help.inc"
C
        DIMENSION F(Llow:lihigh,llow:LJhigh),FHF(LIP1,LJP1)
        DIMENSION FDUM(LIP1,LJqq)
C
#ifdef DEBUG_MODE_ON
      write(*,*) "DEBUG:  in highs.F::YBHI(...)"
#endif
        DO 200 J=-no2+1,mj+no2
        DO 200 I=1,MIP1
  200   FDUM(I,J+ORD2) = F(I,J)
C
C
*        call yhi2(fdum,lvar)
C
        DO 299 J=1,MJP1
        DO 299 I=1,MIP1
  299   FHF(I,J) = 0.
        do 300 k=1,ORD2
        DO 300 J=1,MJP1
        DO 300 I=1,MIP1
  300   FHF(I,J) = FHF(I,J) + FCOEF(K)*
     $          ( FDUM(I,J+K-1) + FDUM(I,J+ORDER-K) )
C
        RETURN
        END
