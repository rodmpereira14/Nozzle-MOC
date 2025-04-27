""" ===============================================
                  INITIALIZATION
=============================================== """

import numpy as np

# ===================================================
# General parameters
# ===================================================

GAMMA = 1.4                       # Specific heat ratio
AE_AT = 3.5              # Nozzle exit to throat area ratio

N_CHARACTERISTICS = 45             # Number of characteristic lines
PRANDTL_MEYER_A = THETA_A = 0.19    # Prandtl-Meyer angle and flow angle at point a (just at the end of throat) (degrees)


R_THROAT = 0.1  # (m)              # Throat height 

P_C = 300000   # [Pa]              # Chamber pressure
T_C = 300      # [K]               # Chamber temperature
P_E = 20000    # [Pa]              # Static pressure ratio at the nozzle exit