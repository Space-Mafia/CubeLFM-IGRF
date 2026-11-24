C LFM-rescue Gate 1 compatibility contract for the exact MIJSET routine.
C This is reconstructed from assignments in src/legacy/subdipol.F; it is not
C represented as a recovered copy of the missing historical mdims.inc.
      INTEGER MI, MI2, MI2P1
      INTEGER MIP1, MIP2, MIP3, MIP4, MIP5
      INTEGER MIM1, MIM2, MIM3, MIM4
      INTEGER MJ, MJ2, MJ2P1, MJ2P2
      INTEGER MJP1, MJP2, MJP3, MJP4, MJP5
      INTEGER MJM1, MJM2, MJM3, MJM4
      COMMON /LFM_RESCUE_MDIMS/ MI, MI2, MI2P1,
     $ MIP1, MIP2, MIP3, MIP4, MIP5,
     $ MIM1, MIM2, MIM3, MIM4,
     $ MJ, MJ2, MJ2P1, MJ2P2,
     $ MJP1, MJP2, MJP3, MJP4, MJP5,
     $ MJM1, MJM2, MJM3, MJM4
      SAVE /LFM_RESCUE_MDIMS/
