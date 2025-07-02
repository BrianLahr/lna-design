# rf_tools/noise_analysis.py
import numpy as np
from typing import List, Tuple
import skrf as rf
from .constants import Z0
from .conversions import db10, from_db10
from .plotting import plot_smith_circle, calc_circle

class NoiseAnalyzer:
    """Class for analyzing noise performance of RF circuits."""
    
    @staticmethod
    def noise_circle(rn: float, gamma_opt: complex, fmin: float, 
                    noise_figure: float, guess: complex) -> Tuple[complex, float]:
        """
        Calculate and plot constant noise circles.
        
        Args:
            rn: Normalized equivalent noise resistance
            gamma_opt: Optimum source reflection coefficient
            fmin: Minimum noise figure (dB)
            noise_figure: Desired noise figure (dB)
            guess: Initial guess point on Smith chart
            
        Returns:
            Tuple of (circle center, circle radius)
        """
        N = ((from_db10(noise_figure) - from_db10(fmin)) / (4 * rn / Z0)) * (abs(1 + gamma_opt)**2)
        c_n = gamma_opt / (1 + N)
        r_n = (1 / (1 + N)) * np.sqrt(N**2 + N - N * (abs(gamma_opt)**2))
        
        # Plotting
        plot_smith_circle(
            center=c_n,
            radius=r_n,
            name=f"Noise {noise_figure}dB",
            color="purple",
            linewidth=2
        )
        plot_smith_circle(
            center=guess,
            radius=0.01,
            name="Guess",
            color="black",
            linewidth=1
        )
        plot_smith_circle(
            center=gamma_opt,
            radius=0.01,
            name="Gamma_opt",
            color="brown",
            linewidth=1
        )
        
        print(f"The optimum source reflection coefficient is {gamma_opt}")
        return c_n, r_n
    
    @staticmethod
    def gain_noise_intersection(gamma_s: complex, noise_figure: float, 
                               s11: complex, gamma_opt: complex, rn: float, 
                               fmin: float) -> List[complex]:
        """
        Find intersections between gain and noise circles.
        
        Args:
            gamma_s: Source reflection coefficient
            noise_figure: Desired noise figure (dB)
            s11: S11 parameter
            gamma_opt: Optimum source reflection coefficient
            rn: Normalized equivalent noise resistance
            fmin: Minimum noise figure (dB)
            
        Returns:
            List of intersection points
        """
        sqabs = lambda x: np.square(np.absolute(x))
        
        # Calculate noise circle
        N = ((from_db10(noise_figure) - from_db10(fmin)) / (4 * rn / Z0)) * (abs(1 + gamma_opt)**2)
        c_n = gamma_opt / (1 + N)
        r_n = (1 / (1 + N)) * np.sqrt(N**2 + N - N * (abs(gamma_opt)**2))
        
        # Calculate gain circle
        gs = ((1 - sqabs(gamma_s)) / sqabs(1 - s11 * gamma_s)) * (1 - sqabs(s11))
        Cs = (gs * np.conjugate(s11)) / (1 - (1 - gs) * sqabs(s11))
        Rs = (np.sqrt(1 - gs) * (1 - sqabs(s11))) / (1 - (1 - gs) * sqabs(s11))
        
        # Plotting
        plot_smith_circle(c_n, r_n, f"Noise {noise_figure}dB", "purple", 2)
        plot_smith_circle(Cs, Rs, "Input Gain", "red", 2)
        plot_smith_circle(gamma_opt, 0.01, "Gamma_opt", "brown", 1)
        
        # Find intersections
        noise_circle = rf.Network(s=calc_circle(c_n, r_n))
        gain_circle = rf.Network(s=calc_circle(Cs, Rs))
        
        intersections = []
        for i, point1 in enumerate(gain_circle.s):
            for point2 in noise_circle.s:
                if np.allclose(point1, point2, atol=0.05):
                    intersections.append(point1[0][0])
        
        # Plot intersections
        for i, point in enumerate(intersections):
            plot_smith_circle(point, 0.02, f"Intersection {i}", linewidth=1)
        
        return intersections