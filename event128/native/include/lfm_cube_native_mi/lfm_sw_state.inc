C Serial global-grid run-mode state.  This is new shell-owned storage, not a
C reconstruction of a missing legacy COMMON block.
      LOGICAL LFM_DIPOLE_BACKGROUND,LFM_IONOSPHERE_ENABLED
      COMMON /LFM_RESCUE_SW_MODE_L/ LFM_DIPOLE_BACKGROUND,
     $ LFM_IONOSPHERE_ENABLED
      SAVE /LFM_RESCUE_SW_MODE_L/
