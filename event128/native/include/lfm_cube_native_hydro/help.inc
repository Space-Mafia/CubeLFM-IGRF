C LFM-rescue Gate 1 compatibility storage for the recovered 2-D kernels.
C Names and initialization relations are proven by init-fortran.F.  The
C historical COMMON block name/order is unavailable, so this new block is
C confined to the compatibility harness until golden parity is established.
      INTEGER ORDER, ORD2
      REAL FCOEF(NO)
      REAL GAMMA, GAMM1, GMM1IN, GMQGM1, ALPHA, PI, PI4, PI8
      REAL PI4INV, PI8INV, CA, CAINV, ALCON, TWELVI
      REAL GMINV, GM1, GP1, GM1INV, WIDTH, GGM1, GGM1IN
      REAL BETA, PDMB
      COMMON /LFM_RESCUE_HELP_I/ ORDER, ORD2
      COMMON /LFM_RESCUE_HELP_R/ FCOEF, GAMMA, GAMM1, GMM1IN,
     $ GMQGM1, ALPHA, PI, PI4, PI8, PI4INV, PI8INV, CA, CAINV,
     $ ALCON, TWELVI, GMINV, GM1, GP1, GM1INV, WIDTH, GGM1,
     $ GGM1IN, BETA, PDMB
      SAVE /LFM_RESCUE_HELP_I/, /LFM_RESCUE_HELP_R/
