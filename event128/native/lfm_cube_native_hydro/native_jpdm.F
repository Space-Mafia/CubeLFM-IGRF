      SUBROUTINE JPDM(F,F0,FRT,FLFT)
#include "param.inc"
#include "mdims.inc"
#include "help.inc"
#include "run-time.inc"
C
C
      DIMENSION F(LIP1,LJP1),F0(LIP1,LJP4),FRT(LIP1,LJP1),
     $          FLFT(LIP1,LJP1)
C
      COMMON /SCRCH1/
     $     DF(LIP1,LJP3),S(LIP1,LJP3),Q(LIP1,LJP2),
     $     DFLFT(LIP1,LJP1),DFRT(LIP1,LJP1)
C
C
#ifdef DEBUG_MODE_ON
      write(*,*) "DEBUG:  in pdms.F::JPDM(...)"
#endif
      DO 200 J = 1,MJP3
      DO 200 I = 1,MIP1
      DF(I,J) = PDMB*(F0(I,J+1)-F0(I,J))
      S(I,J) = SIGN(1.,DF(I,J))
  200 DF(I,J) = ABS(DF(I,J))
C
      DO 300 J = 1,MJP2
      DO 300 I = 1,MIP1
  300 Q(I,J) = ABS(S(I,J) +S(I,J+1))
C
      DO 400 J = 1,MJP1
      DO 400 I = 1,MIP1
      DFLFT(I,J) = F(I,J) - F0(I,J+1)
      DFRT(I,J) = F0(I,J+2) - F(I,J)
      FLFT(I,J) = F(I,J) - S(I,J+1)*AMAX1(0.,
     $     ABS(DFLFT(I,J))- Q(I,J)*DF(I,J))
  400 FRT(I,J) = F(I,J) + S(I,J+1)*AMAX1(0.,
     $     ABS(DFRT(I,J))-Q(I,J+1)*DF(I,J+2))
C
      RETURN
      END
