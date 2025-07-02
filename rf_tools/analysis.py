# rf_tools/analysis.py
import numpy as np
import skrf as rf
import matplotlib.pyplot as plt
import pandas as pd
from typing import List, Dict, Tuple
from .constants import Z0
from .conversions import db10, db20, from_db10, polar_to_rect, rect_to_polar
from .network_analysis import NetworkAnalyzer
from .stability import StabilityAnalyzer
from .noise_analysis import NoiseAnalyzer
from .matching_networks import MatchingNetwork
from .plotting import plot_smith_circle, draw_matching_network

# Add to top of file
from .network_analysis import NetworkAnalyzer  # Fixed import
from .plotting import plot_smith_circle, calc_circle  # Added missing import

class AmplifierAnalysis:
    def __init__(self, network_file: str):
        self.network = rf.Network(network_file)
        self.frequency = self.network.frequency
        self.media = rf.media.DefinedGammaZ0(frequency=self.frequency, z0=Z0)
        self.noise_figure = None  # Initialize attribute
        
    def analyze_stability(self, freq_idx: int = None) -> Dict:
        """
        Analyze stability at a specific frequency or across all frequencies.
        
        Args:
            freq_idx: Frequency index (if None, analyzes all frequencies)
            
        Returns:
            Dictionary containing stability analysis results
        """
        if freq_idx is not None:
            s11, s12 = self.network.s[freq_idx].tolist()[0]
            s21, s22 = self.network.s[freq_idx].tolist()[1]
            return {
                'rollet': StabilityAnalyzer.rollet_condition(s11, s12, s21, s22),
                'mu': StabilityAnalyzer.mu_stability(s11, s12, s21, s22)
            }
        else:
            results = []
            for i in range(len(self.frequency.f)):
                s11, s12 = self.network.s[i].tolist()[0]
                s21, s22 = self.network.s[i].tolist()[1]
                results.append({
                    'frequency': self.frequency.f[i],
                    'mu': StabilityAnalyzer.mu_stability(s11, s12, s21, s22, verbose=False)
                })
            return results
    
    def plot_stability_over_frequency(self, reference_network: rf.Network = None):
        """
        Plot stability (mu) over frequency range.
        
        Args:
            reference_network: Optional reference network for comparison
        """
        freq_ghz = self.frequency.f / 1e9
        mu_values = [result['mu'] for result in self.analyze_stability()]
        
        plt.figure(figsize=(10, 8), dpi=80)
        plt.grid(visible=True, which='major', axis='both')
        plt.plot(freq_ghz, np.ones(len(freq_ghz)), label='Stability threshold')
        plt.plot(freq_ghz, mu_values, label='Current network')
        
        if reference_network:
            ref_analysis = AmplifierAnalysis(reference_network)
            ref_mu = [result['mu'] for result in ref_analysis.analyze_stability()]
            plt.plot(freq_ghz, ref_mu, label='Reference network', color='red')
        
        plt.title(f'{self.frequency.start/1e9} To {self.frequency.stop/1e9} GHz')
        plt.xlabel('Frequency (GHz)')
        plt.ylabel('Mu')
        plt.legend(loc=4)
        plt.xticks(np.arange(min(freq_ghz), max(freq_ghz)+1, 1))
        plt.show()
    
    def make_unconditionally_stable(self, freq_idx: int, rin: float = None, rout: float = None):
        """
        Add resistors to make the network unconditionally stable.
        
        Args:
            freq_idx: Frequency index for design
            rin: Input resistance (if None, calculates automatically)
            rout: Output resistance (if None, calculates automatically)
        """
        s11, s12 = self.network.s[freq_idx].tolist()[0]
        s21, s22 = self.network.s[freq_idx].tolist()[1]
        
        if rin is None or rout is None:
            stability = StabilityAnalyzer.rollet_condition(s11, s12, s21, s22)
            rin = (Z0 * 0.4).real  # Example calculation - adjust as needed
            rout = (Z0 * 0.22).real  # Example calculation - adjust as needed
        
        self.network = self.media.resistor(rin) ** self.network ** self.media.resistor(rout)
    
    def analyze_noise(self, freq_idx: int, rn: float, gamma_opt: complex, 
                    fmin: float, noise_figure: float, guess: complex):
        """
        Analyze noise performance at a specific frequency.
        
        Args:
            freq_idx: Frequency index
            rn: Normalized equivalent noise resistance
            gamma_opt: Optimum source reflection coefficient
            fmin: Minimum noise figure (dB)
            noise_figure: Desired noise figure (dB)
            guess: Initial guess point on Smith chart
        """
        return NoiseAnalyzer.noise_circle(
            rn=rn,
            gamma_opt=gamma_opt,
            fmin=fmin,
            noise_figure=noise_figure,
            guess=guess
        )
    
    def design_matching_network(self, freq_idx: int, zs: complex, zl: complex):
        """
        Design matching networks for given source and load impedances.
        
        Args:
            freq_idx: Frequency index
            zs: Source impedance
            zl: Load impedance
            
        Returns:
            Tuple of (input_matching, output_matching) networks
        """
        freq_hz = self.frequency.f[freq_idx]
        matcher = MatchingNetwork(self.frequency)
        
        input_match = matcher.l_match(Z0, np.conj(zs), freq_hz)
        output_match = matcher.l_match(np.conj(zl), Z0, freq_hz)
        
        return input_match, output_match
    
    def analyze_sensitivity(self, temperatures: List[float], snr: float) -> pd.DataFrame:
        """
        Analyze sensitivity across different temperatures.
        
        Args:
            temperatures: List of temperatures in Celsius
            snr: Required signal-to-noise ratio (dB)
            
        Returns:
            Pandas DataFrame with sensitivity analysis results
        """
        freq_step = self.frequency.f[1] - self.frequency.f[0]
        s21_list = self.network.s21.s
        
        results = {
            'Temp C': temperatures,
            'Bandwidth GHz': [],
            'KTB': [],
            'Thermal_noise_floor dBm': [],
            'NF dB': [],
            'Sensitivity dBm': []
        }
        
        for temp in temperatures:
            bandwidth = NetworkAnalyzer.bandwidth(s21_list, freq_step)
            ktb = NetworkAnalyzer.ktb(temp, bandwidth)
            noisefloor = NetworkAnalyzer.thermal_noise_floor(temp, bandwidth)
            sens = NetworkAnalyzer.sensitivity(temp, bandwidth, self.noise_figure, snr)
            
            results['Bandwidth GHz'].append(bandwidth[0]/1e9)
            results['KTB'].append(ktb[0])
            results['Thermal_noise_floor dBm'].append(noisefloor[0])
            results['Sensitivity dBm'].append(sens[0])
            results['NF dB'].append(self.noise_figure)
        
        return pd.DataFrame.from_dict(results)