# rf_tools/conversions.py
import numpy as np
from typing import Union, Tuple
from cmath import log10, rect
from .constants import C_LIGHT, Z0

def celsius_to_kelvin(celsius: float) -> float:
    """Convert temperature from Celsius to Kelvin."""
    return celsius + 273.15

def polar_to_rect(amplitude: float, angle: float) -> complex:
    """Convert polar coordinates to rectangular form."""
    nprect = np.vectorize(rect)
    return nprect(amplitude, np.deg2rad(angle))

def rect_to_polar(x: complex) -> Tuple[float, float]:
    """Convert rectangular coordinates to polar form."""
    return abs(x), np.angle(x, deg=True)

# Make dB conversions numpy-friendly
def db10(x: Union[float, complex, np.ndarray]) -> Union[float, np.ndarray]:
    return 10 * np.log10(np.abs(x))

# Add warning for invalid inputs
def from_db10(x: float) -> float:
    if x > 100:  # Unrealistically high dB value
        import warnings
        warnings.warn(f"Unusually high dB value: {x}")
    return 10**(x/10)

def dbm10(x: Union[float, complex]) -> float:
    """Convert to power dBm (dB relative to 1mW)."""
    return db10(x) + 30

def db20(x: Union[float, complex]) -> float:
    """Convert to voltage dB (20*log10)."""
    x = abs(x)
    return 20 * log10(x)

def from_db20(x: float) -> float:
    """Convert from voltage dB to linear scale."""
    return round(10**(x/20), 3)

# Add these commonly needed conversions:
def wavelength(freq: float) -> float:
    """Calculate wavelength in meters for given frequency in Hz"""
    return C_LIGHT / freq

def vswr_to_gamma(vswr: float) -> float:
    """Convert VSWR to reflection coefficient"""
    return (vswr - 1) / (vswr + 1)

def normalize_impedance(z: complex, z0: float = Z0) -> complex:
    """Normalize impedance to reference Z0"""
    return z / z0

def denormalize_impedance(z_norm: complex, z0: float = Z0) -> complex:
    """Convert normalized impedance to ohms"""
    return z_norm * z0


