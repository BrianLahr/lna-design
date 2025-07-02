# rf_tools/network_analysis.py
import numpy as np
from typing import List, Tuple
from skrf import Network
from .constants import Z0
from .conversions import db10

class NetworkAnalyzer:
    """Class for analyzing RF networks using S-parameters."""
    
    @staticmethod
    def bandwidth(s21_list: List[complex], freq_step: float) -> float:
        """Calculate equivalent noise bandwidth from S21 parameters.
        
        Args:
            s21_list: List of complex S21 values across frequency
            freq_step: Frequency step size in Hz
            
        Returns:
            Noise bandwidth in Hz
        """
        s21_array = np.asarray(s21_list)  # Safer than array()
        if np.any(np.iscomplex(s21_array)):
            s21_array = np.abs(s21_array)
        return np.trapz(s21_array**2, dx=freq_step) / np.max(s21_array)**2
    
    # Add to both files
    def _validate_s_params(s_params):
        """Ensure S-parameters are physically realizable."""
        if not all(abs(s) < 10 for s in s_params):  # Arbitrary large limit
            raise ValueError("Unrealistic S-parameter magnitude")

    # In network_analysis.py
    @staticmethod
    def calculate_delta(s11: complex, s12: complex, 
                    s21: complex, s22: complex) -> complex:
        """Calculate S-parameter determinant."""
        return s11*s22 - s12*s21

    @staticmethod
    def stability_factor(s11: complex, s12: complex, 
                        s21: complex, s22: complex) -> float:
        """Calculate Rollet's K-factor for stability analysis.
        
        Args:
            s11: Input reflection coefficient
            s12: Reverse transmission coefficient
            s21: Forward transmission coefficient
            s22: Output reflection coefficient
            
        Returns:
            K-factor (device is unconditionally stable if K > 1)
            
        Formula:
            K = (1 - |S11|² - |S22|² + |Δ|²) / (2|S12S21|)
            where Δ = S11S22 - S12S21
        """
        delta = s11*s22 - s12*s21
        numerator = 1 - abs(s11)**2 - abs(s22)**2 + abs(delta)**2
        denominator = 2 * abs(s12*s21)
        return numerator / denominator

    @staticmethod
    def maximum_stable_gain(s_params: np.ndarray) -> float:
        """Calculate MSG for potentially unstable devices."""
        s11, s12, s21, s22 = s_params
        return db10(abs(s21/s12))
    
    @staticmethod
    def swr(reflection_coeff: complex) -> float:
        """Calculate Standing Wave Ratio from reflection coefficient."""
        gamma = abs(reflection_coeff)
        return (1 + gamma) / (1 - gamma)
    
    @staticmethod
    def transducer_gain(s11: complex, s12: complex, s21: complex, 
                       s22: complex, gamma_s: complex, gamma_l: complex) -> float:
        """Calculate transducer gain of a network."""
        qabs = lambda x: np.square(np.absolute(x))
        numerator = qabs(s21) * (1 - qabs(gamma_s)) * (1 - qabs(gamma_l))
        denominator = qabs((1 - s11*gamma_s)*(1 - s22*gamma_l) - s12*s21*gamma_l*gamma_s)
        return numerator / denominator
    
    @staticmethod
    def input_impedance(zs: complex) -> complex:
        """Calculate input impedance for negative resistance circuits.
        
        Args:
            zs: Source impedance (Ohms)
            
        Returns:
            Input impedance (Ohms)
            
        Raises:
            ValueError: If real part is non-negative
        """
        zin = -zs
        if zin.real >= 0:
            raise ValueError(f"Invalid impedance {zin}: Real part must be negative")
        return zin
    
    @staticmethod
    def output_impedance(zl: complex) -> complex:
        """Calculate output impedance from load impedance."""
        zout = -zl
        if zout.real < 0:
            return zout
        raise ValueError("Output impedance real part must be negative")
    
    @staticmethod
    def impedance_from_gamma(gamma: complex, z0: float = Z0) -> Tuple[complex, complex]:
        """Calculate impedance from reflection coefficient."""
        return (z0 * ((gamma + 1) / (1 - gamma)), 
                z0 * ((-gamma + 1) / (1 + gamma)))