# rf_tools/stability.py
import numpy as np
import math
from typing import Tuple
from .conversions import db10

class StabilityAnalyzer:
    """Class for analyzing stability of RF networks."""
    
    @staticmethod
    def maximum_available_gain(s11: complex, s12: complex, 
                             s21: complex, s22: complex) -> float:
        """
        Calculate Maximum Available Gain for unconditionally stable devices.
        
        Raises:
            ValueError: If the device is not unconditionally stable
        """
        delta = (s11 * s22) - (s12 * s21)
        amp_s11 = abs(s11)**2
        amp_s22 = abs(s22)**2
        amp_delta = abs(delta)**2
        
        k = (1 - amp_s11 - amp_s22 + amp_delta) / (2 * abs(s12 * s21))
        if k < 1:
            raise ValueError("Device is not unconditionally stable")
            
        b1 = 1 + amp_s11 - amp_s22 - amp_delta
        if b1 < 0:
            return db10(abs(s21 / s12)) + db10(abs(k + math.sqrt(k**2 - 1)))
        return db10(abs(s21 / s12)) + db10(abs(k - math.sqrt(k**2 - 1)))
    
    @staticmethod
    def rollet_condition(s11: complex, s12: complex, 
                        s21: complex, s22: complex) -> dict:
        """Analyze stability using Rollet's condition."""
        delta = (s11 * s22) - (s12 * s21)
        amp_s11 = abs(s11)**2
        amp_s22 = abs(s22)**2
        amp_delta = abs(delta)**2
        
        k = (1 - amp_s11 - amp_s22 + amp_delta) / (2 * abs(s12 * s21))
        
        result = {
            'delta': delta,
            'k_factor': k,
            'delta_magnitude': abs(delta),
            'is_stable': (amp_s11 < 1 and amp_s22 < 1 and 
                          abs(delta) < 1 and k > 1)
        }
        
        if result['is_stable']:
            result['message'] = 'Unconditionally stable'
        else:
            result['message'] = 'Potentially unstable - check stability circles'
        
        return result
    
    @staticmethod
    def mu_stability(s_params: np.ndarray, verbose: bool = True) -> Tuple[float, float]:
        """Calculate mu and mu' stability factors.
        
        Returns:
            Tuple of (mu_input, mu_output)
        """
        s11, s12, s21, s22 = s_params
        delta = s11*s22 - s12*s21
        
        mu1 = (1 - abs(s11)**2) / (abs(s22 - delta*s11.conj()) + abs(s12*s21))
        mu2 = (1 - abs(s22)**2) / (abs(s11 - delta*s22.conj()) + abs(s12*s21))
        
        if verbose:
            print(f"Input stability μ: {mu1:.3f}")
            print(f"Output stability μ: {mu2:.3f}")
            
        return min(mu1, mu2)
    
    @staticmethod
    def input_stability_circle(s_params: np.ndarray) -> Tuple[complex, float]:
        """Calculate input stability circle (center, radius)."""
        s11, s12, s21, s22 = s_params
        delta = s11*s22 - s12*s21
        center = (s11 - delta*s22.conj()).conj() / (abs(s11)**2 - abs(delta)**2)
        radius = abs(s12*s21) / abs(abs(s11)**2 - abs(delta)**2)
        return center, radius
    
    # In stability.py
    @classmethod
    def check_stability(cls, s_params: np.ndarray) -> dict:
        """Comprehensive stability analysis."""
        results = cls.rollet_condition(*s_params)
        results['mu'] = cls.mu_stability(*s_params, verbose=False)
        results['msg'] = cls.maximum_stable_gain(s_params) if not results['is_stable'] else None
        return results
    
    # Add to both files
    def _validate_s_params(s_params):
        """Ensure S-parameters are physically realizable."""
        if not all(abs(s) < 10 for s in s_params):  # Arbitrary large limit
            raise ValueError("Unrealistic S-parameter magnitude")