#!/usr/bin/env python3
"""
RF Amplifier Design Tool - Main Entry Point
"""
from rf_tools.analysis import AmplifierAnalysis
from rf_tools.conversions import polar_to_rect
import matplotlib.pyplot as plt

def main():
    """Main execution workflow"""
    try:
        # Initialize with S-parameter file
        analyzer = AmplifierAnalysis('data/f360772a.s2p')
        analyzer.noise_figure = 0.4  # Set default noise figure
        
        # Configure analysis at 10 GHz
        freq_idx = 10
        noise_params = {
            'rn': 0.09 * 50,
            'gamma_opt': polar_to_rect(0.66, 102),
            'fmin': 0.37,
            'noise_figure': 0.4,
            'guess': 0.3309 + 0.5188j
        }
        
        # Run analyses
        analyzer.make_unconditionally_stable(freq_idx)
        analyzer.analyze_noise(freq_idx, **noise_params)
        
        # Generate outputs
        analyzer.plot_stability_over_frequency()
        plt.show()
        
    except FileNotFoundError:
        print("Error: S-parameter file not found")
    except Exception as e:
        print(f"Analysis failed: {str(e)}")

if __name__ == "__main__":
    main()