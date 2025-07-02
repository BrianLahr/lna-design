# Explicit version control
__version__ = "1.0.0"

# Export conversion utilities too
from .conversions import (
    db10, db20, dbm10,
    from_db10, from_db20,
    polar_to_rect, rect_to_polar
)

# Add package-level documentation
__doc__ = """
RF Tools - A Python package for microwave network analysis

Features:
- Amplifier stability analysis
- Noise figure calculations
- Impedance matching
- S-parameter manipulation
"""

# Group related exports
__all__ = [
    # Core classes
    'AmplifierAnalysis',
    'NoiseAnalyzer',
    'StabilityAnalyzer',
    
    # Utilities
    'network_analysis',
    'matching_networks'
]