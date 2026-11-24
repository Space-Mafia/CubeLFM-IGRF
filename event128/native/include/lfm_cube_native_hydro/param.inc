C LFM-rescue Gate 1 compile-harness subset of the missing param.inc.
C Values are configurable test maxima, not recovered production dimensions.
      INTEGER NI, NJ, NK, NO, NO2, LI, LJ, LK
      INTEGER LLOW, LIHIGH, LJHIGH
      INTEGER NIP1, NJP1, NKP1
      INTEGER LIP1, LIP2, LIP3, LIP4, LIP5
      INTEGER LJP1, LJP2, LJP3, LJP4, LJP5
      INTEGER LIQQ, LJQQ
      PARAMETER (NI=64, NJ=64,
     $ NK=64, NO=8,
     $ NO2=NO/2, LI=64,
     $ LJ=64, LK=64)
      PARAMETER (LLOW=-NO/2+1, LIHIGH=LI+NO/2,
     $ LJHIGH=LJ+NO/2)
      PARAMETER (NIP1=NI+1, NJP1=NJ+1, NKP1=NK+1)
      PARAMETER (LIP1=LI+1, LIP2=LI+2, LIP3=LI+3,
     $ LIP4=LI+4, LIP5=LI+5)
      PARAMETER (LJP1=LJ+1, LJP2=LJ+2, LJP3=LJ+3,
     $ LJP4=LJ+4, LJP5=LJ+5)
C LIQQ/LJQQ are only capacity bounds in the recovered kernels.  The exact
C historical rest_of_param.inc spellings remain unresolved.
      PARAMETER (LIQQ=LI+NO, LJQQ=LJ+NO)
