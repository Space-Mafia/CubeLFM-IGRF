      SUBROUTINE PDM2(F,F0,Fq00,Fq01,Fq10,Fq11)
#include "param.inc"
#include "mdims.inc"
#include "help.inc"
#include "run-time.inc"
C
C
        DIMENSION F(LIP1,LJP1),F0(LIP4,LJP4),Fq00(LIP1,LJP1),
     $          Fq01(LIP1,LJP1),Fq10(LIP1,LJP1),Fq11(LIP1,LJP1)
C
      COMMON /SCRCH1/
     $     DF(LIP3,LJP3),DF0(LIP3,LJP3),S(LIP3,LJP3),
     $     q(LIP2,LJP2),DFq00(LIP1,LJP1),DFq01(LIP1,LJP1),
     $     DFq10(LIP1,LJP1),DFq11(LIP1,LJP1),S0(LIP3,LJP3)
C
#ifdef DEBUG_MODE_ON
      write(*,*) "DEBUG:  in pdms.F::PDM2(...)"
#endif
C
C
      DO 200 J = 1,MJp2
      DO 200 I = 1,MIP3
        DF0(I,J) = F0(I+1,J+1)-F0(I,J+1)
        S(I,J) = SIGN(1.0,DF0(I,J))
        DF0(I,J) = ABS(DF0(I,J))
        DF(I,J) = PDMB*DF0(I,J)
        S0(I,J) = 0.5*S(I,J)
  200   CONTINUE
C
      DO 300 J = 1,MJP2
      DO 300 I = 1,MIP2
  300 q(I,J) = ABS(S(I,J) +S(I+1,J))
C
      DO 400 J = 1,MJP1
      DO 400 I = 1,MIP1
        DFq00(I,J) = F(I,J) - F0(I+1,J+1)
        DFq01(I,J) = F(i,J) - F0(I+1,J+2)
        DFq10(I,J) = F0(I+2,J+1) - F(I,J)
        DFq11(I,J) = F0(I+2,J+2) - F(I,J)
        Fq00(I,J) = F0(I+1,J+1) + S(I+1,J)*AMIN1( DF0(I+1,J),
     $          q(I,J)*DF(I,J),
     $          ABS(DFq00(I,J)*(SIGN(0.5,DFq00(I,J))+S0(I+1,J))))
        Fq01(I,J) = F0(I+1,J+2) + S(I+1,J+1)*AMIN1( DF0(I+1,J+1),
     $          q(I,J+1)*DF(I,J+1),
     $          ABS(DFq01(I,J)*(SIGN(0.5,DFq01(I,J))+S0(I+1,J+1))))
        Fq10(I,J) = F0(I+2,J+1) - S(I+1,J)*AMIN1( DF0(I+1,J),
     $          q(I+1,J)*DF(I+2,J),
     $          ABS(DFq10(I,J)*(SIGN(0.5,DFq10(I,J))+S0(I+1,J))))
        Fq11(I,J) = F0(I+2,J+2) - S(I+1,J+1)*AMIN1( DF0(I+1,J+1),
     $          q(I+1,J+1)*DF(I+2,J+1),
     $          ABS(DFq11(I,J)*(SIGN(0.5,DFq11(I,J))+S0(I+1,J+1))))
  400   CONTINUE
C
C
C
      DO 500 J = 1,MJP3
      DO 500 I = 1,MIP2
        DF0(I,J) = F0(I+1,J+1)-F0(I+1,J)
        S(I,J) = SIGN(1.0,DF0(I,J))
        DF0(I,J) = ABS(DF0(I,J))
        S0(I,J) = 0.5*S(I,J)
        DF(I,J) = PDMB*DF0(I,J)
  500   CONTINUE
C
      DO 600 J = 1,MJP2
      DO 600 I = 1,MIP2
  600 q(I,J) = ABS(S(I,J) +S(I,J+1))
C
C
      DO 700 J = 1,MJP1
      DO 700 I = 1,MIP1
        DFq00(I,J) = F(I,J) - Fq00(I,J)
        DFq10(I,J) = F(I,J) - Fq10(I,J)
        DFq01(I,J) = Fq01(I,J) - F(I,J)
        DFq11(I,J) = Fq11(I,J) - F(I,J)
        Fq00(I,J) = Fq00(I,J) + S(I,J+1)*AMIN1( DF0(I,J+1),
     $          q(I,J)*DF(I,J),
     $          ABS(DFq00(I,J)*(SIGN(0.5,DFq00(I,J))+S0(I,J+1))))
        Fq01(I,J) = Fq01(I,J) - S(I,J+1)*AMIN1( DF0(I,J+1),
     $          q(I,j+1)*DF(I,j+2),
     $          ABS(DFq01(I,J)*(SIGN(0.5,DFq01(I,J))+S0(I,J+1))))
        Fq10(I,J) = Fq10(I,J) + S(I+1,J+1)*AMIN1( DF0(I+1,J+1),
     $          q(I+1,J)*DF(I+1,J),
     $          ABS(DFq10(I,J)*(SIGN(0.5,DFq10(I,J))+S0(I+1,J+1))))
        Fq11(I,J) = Fq11(I,J) - S(I+1,J+1)*AMIN1( DF0(I+1,J+1),
     $          q(I+1,J+1)*DF(I+1,J+2),
     $          ABS(DFq11(I,J)*(SIGN(0.5,DFq11(I,J))+S0(I+1,J+1))))
  700   CONTINUE
C
C
C
C
        RETURN
        END
