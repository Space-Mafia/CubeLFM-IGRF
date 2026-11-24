        SUBROUTINE HIGH2(F,FqCORN,LVAR)
#include "param.inc"
#include "global_dims.inc"
#include "mdims.inc"
#include "help.inc"
C
        DIMENSION F(llow:LIhigh,llow:LJhigh),FqCORN(LIP1,LJP1)
        DIMENSION FHF(LIP1,llow:LJhigh)
        DIMENSION FXDUM(llow:liqq,llow:ljqq),
     $       FYDUM(Llow:liqq,llow:ljqq)
        EQUIVALENCE (FXDUM,FYDUM)
C
#ifdef DEBUG_MODE_ON
      write(*,*) "DEBUG:  in highs.F::HIGH2(...)"
#endif
        DO 200 J=-no2+1,MJ+no2
        DO 200 I=-no2+1,MI+no2
  200   FXDUM(I+ord2,J) = F(I,J)
C
*        CALL IHI2(FXDUM,LVAR)
C
        DO 299 J=-no2+1,MJ+no2
        DO 299 I=1,MIP1
  299   FHF(I,J) = 0.
        DO 300 K=1,ORD2
        DO 300 J=-no2+1,MJ+no2
        DO 300 I=1,MIP1
  300   FHF(I,J) = FHF(I,J) + FCOEF(K)*
     $          ( FXDUM(I+K-1,J) + FXDUM(I+ORDER-K,J) )
C
C
        DO 400 J=-no2+1,MJ+no2
        DO 400 I=1,MIP1
  400   FYDUM(I,J+ord2) = FHF(I,J)
C
*        CALL JHI2(FYDUM,LVAR)
C
        DO 499 J=1,MJP1
        DO 499 I=1,MIP1
  499   FqCORN(I,J) = 0.
        do 500 k=1,ord2
        DO 500 J=1,MJP1
        DO 500 I=1,MIP1
  500   FqCORN(I,J) = FqCORN(I,J) + FCOEF(K)*
     $          ( FYDUM(I,J+K-1) + FYDUM(I,J+ORDER-K) )
C
        RETURN
        END
