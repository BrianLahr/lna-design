# rf_tools/matching_networks.py
import numpy as np
import skrf as rf
from .constants import Z0
from typing import Tuple, Dict

class MatchingNetwork:
    """Class for designing and analyzing impedance matching networks."""
    
    def __init__(self, frequency, z0=Z0):
        """
        Initialize matching network designer.
        
        Args:
            frequency: Frequency object from skrf
            z0: Characteristic impedance (default: 50 ohms)
        """
        self.frequency = frequency
        self.z0 = z0
        self.media = rf.media.DefinedGammaZ0(frequency=frequency, z0=z0)
    
    @staticmethod
    def series_capacitor(freq: float, x: float) -> float:
        """Calculate series capacitance from reactance."""
        w = 2 * np.pi * freq
        return 1 / (w * x * Z0)
    
    @staticmethod
    def series_inductor(freq: float, x: float) -> float:
        """Calculate series inductance from reactance."""
        w = 2 * np.pi * freq
        return (x * Z0) / w
    
    @staticmethod
    def shunt_capacitor(freq: float, b: float) -> float:
        """Calculate shunt capacitance from susceptance."""
        w = 2 * np.pi * freq
        return b / (w * Z0)
    
    @staticmethod
    def shunt_inductor(freq: float, b: float) -> float:
        """Calculate shunt inductance from susceptance."""
        w = 2 * np.pi * freq
        return Z0 / (w * b)
    
    def l_match(self, zs: complex, zl: complex) -> Dict:
        """
        Design an L-section matching network.
        
        Args:
            zs: Source impedance
            zl: Load impedance
            
        Returns:
            Dictionary containing matching network parameters and components
        """
        rs, rl = zs.real, zl.real
        
        # Determine absorption type and component
        if zs.imag == 0:
            absorb = 'load'
            if zl.imag < 0:
                restype = 'inductor'
                res_value = self.shunt_inductor(self.frequency.f, -zl.imag)
                res_network = self.media.shunt_inductor(res_value * 1e-9)
            else:
                restype = 'capacitor'
                res_value = self.shunt_capacitor(self.frequency.f, zl.imag)
                res_network = self.media.shunt_capacitor(res_value * 1e-12)
        else:
            absorb = 'source'
            if zs.imag < 0:
                restype = 'inductor'
                res_value = self.shunt_inductor(self.frequency.f, -zs.imag)
                res_network = self.media.shunt_inductor(res_value * 1e-9)
            else:
                restype = 'capacitor'
                res_value = self.shunt_capacitor(self.frequency.f, zs.imag)
                res_network = self.media.shunt_capacitor(res_value * 1e-12)
        
        # Calculate matching components
        if rs > rl:
            m = rs / rl
            q = np.sqrt(m - 1)
            xs = q * rl
            xp = -xs * (1 + 1 / (q**2))
        else:
            m = rl / rs
            q = np.sqrt(m - 1)
            xp = rl / q
            xs = -xp / (1 + 1 / (q**2))
        
        # Determine component types
        xskind = 'Inductor' if xs > 0 else 'Capacitor'
        xpkind = 'Inductor' if xp > 0 else 'Capacitor'
        
        xs, xp = abs(xs), abs(xp)
        
        # Create matching network
        if rs > rl:
            if xskind == 'Inductor':
                jx1 = self.shunt_inductor(self.frequency.f, xs)
                jx2 = self.series_capacitor(self.frequency.f, xp)
                matching_network = self.media.shunt_inductor(jx1 * 1e-9) ** self.media.capacitor(jx2 * 1e-12)
            else:
                jx1 = self.shunt_capacitor(self.frequency.f, xs)
                jx2 = self.series_inductor(self.frequency.f, xp)
                matching_network = self.media.shunt_capacitor(jx1 * 1e-12) ** self.media.inductor(jx2 * 1e-9)
        else:
            if xskind == 'Inductor':
                jx1 = self.series_inductor(self.frequency.f, xs)
                jx2 = self.shunt_capacitor(self.frequency.f, xp)
                matching_network = self.media.inductor(jx1 * 1e-9) ** self.media.shunt_capacitor(jx2 * 1e-12)
            else:
                jx1 = self.series_capacitor(self.frequency.f, xs)
                jx2 = self.shunt_inductor(self.frequency.f, xp)
                matching_network = self.media.capacitor(jx1 * 1e-12) ** self.media.shunt_inductor(jx2 * 1e-9)
        
        # Combine with absorption network
        if absorb == 'source':
            matching_network = res_network ** matching_network
        else:
            matching_network = matching_network ** res_network
        
        return {
            'network': matching_network,
            'components': {
                'jx1': (jx1, xskind),
                'jx2': (jx2, xpkind),
                'res': (res_value, restype)
            },
            'parameters': {
                'Q': q,
                'matching_type': 'L-section',
                'absorption': absorb
            }
        }