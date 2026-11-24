      SUBROUTINE XBHI(F,FHF,LVAR)
#include "param.inc"
#include "mdims.inc"
#include "help.inc"
C
      DIMENSION F(llow:lihigh,Llow:ljhigh),FHF(LIP1,LJP1)
      DIMENSION FDUM(LIqq,LJP1)
C
#ifdef DEBUG_MODE_ON
      write(*,*) "DEBUG:  in highs.F::XBHI(...)"
#endif
      DO 200 J=1,MJP1
      DO 200 I=-no2+1,MI+no2
      FDUM(I+ord2,J) = F(I,J)
  200 continue
C
C
*      call xhi2(fdum,lvar)
C
C
      DO 299 J=1,MJP1
      DO 299 I=1,MIP1
      FHF(I,J) = 0.
 299  continue
c
      DO 300 K=1,ORD2
      DO 300 J=1,MJP1
      DO 300 I=1,MIP1
      FHF(I,J) = FHF(I,J) + FCOEF(K)*
     $    ( FDUM(I+K-1,J) + FDUM(I+ORDER-K,J) )
 300    continue
      RETURN
      END
