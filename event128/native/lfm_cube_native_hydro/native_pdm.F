      SUBROUTINE PDM(F,F0,FRT,FLFT)
#include "param.inc"
#include "mdims.inc"
#include "help.inc"
#include "run-time.inc"
C
C
      DIMENSION F(LIP1,LJ),F0(LIP4,LJ),FRT(LIP1,LJ),
     $          FLFT(LIP1,LJ)
C
      COMMON /SCRCH1/ 
     $     DF(LIP3,LJ),S(LIP3,LJ),Q(LIP2,LJ),
     $     DFLFT(LIP1,LJ),DFRT(LIP1,LJ)
C
#ifdef DEBUG_MODE_ON
      write(*,*) "DEBUG:  in pdms.F::PDM(...)"
#endif
      DO 200 J = 1,MJ
      DO 200 I = 1,MIP3
      DF(I,J) = PDMB*(F0(I+1,J)-F0(I,J))
      S(I,J) = SIGN(1.,DF(I,J))
  200 DF(I,J) = ABS(DF(I,J))
C
      DO 300 J = 1,MJ
      DO 300 I = 1,MIP2
  300 Q(I,J) = ABS(S(I,J) +S(I+1,J))
C
      DO 400 J = 1,MJ
      DO 400 I = 1,MIP1
      DFLFT(I,J) = F(I,J) - F0(I+1,J)
      DFRT(I,J) = F0(I+2,J) - F(I,J)
      FLFT(I,J) = F(I,J) - S(I+1,J)*AMAX1(0.,
     $     ABS(DFLFT(I,J))- Q(I,J)*DF(I,J))
  400 FRT(I,J) = F(I,J) + S(I+1,J)*AMAX1(0.,
     $     ABS(DFRT(I,J))-Q(I+1,J)*DF(I+2,J))
C
      RETURN
      END
