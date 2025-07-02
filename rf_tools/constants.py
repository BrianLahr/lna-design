# rf_tools/constants.py
import numpy as np

# Physical constants
BOLTZMANN_CONST = 1.38065e-23
Z0 = 50  # Standard characteristic impedance

# Add these common RF constants:
C_LIGHT = 299792458  # Speed of light (m/s)
MU_0 = 4e-7 * np.pi  # Permeability of free space
EPSILON_0 = 1/(MU_0 * C_LIGHT**2)  # Permittivity of free space

# Standard reference temperatures
T_ROOM = 290  # Kelvin (17°C)
T_AMBIENT = 300  # Kelvin (27°C)

# Consider adding typical frequency ranges
RF_RANGES = {
    'HF': (3e6, 30e6),
    'VHF': (30e6, 300e6),
    'UHF': (300e6, 3e9),
    'SHF': (3e9, 30e9)
}