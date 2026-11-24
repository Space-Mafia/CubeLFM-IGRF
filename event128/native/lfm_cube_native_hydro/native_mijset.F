        subroutine mijset(i,j)
#include "mdims.inc"
#ifdef DEBUG_MODE_ON
      write(*,*) "DEBUG:  in subdipol.F::mijset(...)"
#endif
        mi =i
        mi2 = mi/2
        mi2p1 = mi2+1
        mip1 = mi+1
        mip2 = mip1+1
        mip3 = mip2+1
        mip4 = mip3+1
        mip5 = mip4+1
        mim1 = mi-1
        mim2 = mim1-1
        mim3 = mim2-1
        mim4 = mim3-1
C
        mj =j
        mj2 = mj/2
        mj2p1 = mj2+1
        mj2p2 = mj2+2
        mjp1 = mj+1
        mjp2 = mjp1+1
        mjp3 = mjp2+1
        mjp4 = mjp3+1
        mjp5 = mjp4+1
        mjm1 = mj-1
        mjm2 = mjm1-1
        mjm3 = mjm2-1
        mjm4 = mjm3-1
C
        return
        end
