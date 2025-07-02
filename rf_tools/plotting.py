# rf_tools/plotting.py
import schemdraw
import schemdraw.elements as elm
import skrf as rf
import matplotlib.pyplot as plt
from typing import Optional, Tuple, List
from .constants import Z0
import numpy as np

def calc_circle(c, r):
  theta = np.linspace(0, 2*np.pi, 1000)
  return c + r*np.exp(1.0j*theta)

def plot_smith_circle(center: complex, radius: float, name: str = None, 
                     color: str = 'blue', linewidth: int = 1, 
                     figsize: Tuple[int, int] = (10, 8), dpi: int = 80):
    """
    Plot a circle on the Smith chart.
    
    Args:
        center: Center of the circle (complex number)
        radius: Radius of the circle
        name: Name for the circle legend
        color: Color of the circle
        linewidth: Line width
        figsize: Figure size
        dpi: Figure resolution
    """
    plt.figure(figsize=figsize, dpi=dpi)
    n = rf.Network(name=name, s=calc_circle(center, radius))
    n.plot_s_smith(lw=linewidth, draw_labels=True, color=color)


def draw_matching_network(zs: complex, zl: complex, rin: float, rout: float,
                         input_pkg: Tuple, output_pkg: Tuple) -> schemdraw.Drawing:
    """
    Draw a schematic of the matching network.
    
    Args:
        zs: Source impedance
        zl: Load impedance
        rin: Input resistance
        rout: Output resistance
        input_pkg: Input package parameters
        output_pkg: Output package parameters
        
    Returns:
        schemdraw.Drawing object
    """
    _, res_in, jx1_in, jx2_in, xskind_in, xpkind_in, restype_in, absorb_in = input_pkg
    _, res_out, jx1_out, jx2_out, xskind_out, xpkind_out, restype_out, absorb_out = output_pkg
    
    rs, rl = zs.real, zl.real
    d = schemdraw.Drawing()
    
    # Source side
    d += elm.Ground()
    d += (V1 := elm.SourceV())
    d += elm.Resistor().right().label(f'{Z0} Ohm')
    
    # Input matching
    if absorb_in == 'source':
        if restype_in == 'inductor':
            d.push()
            d += elm.Inductor().down().label(f'{res_in*1e3:.3f} nH')
            d += elm.Ground()
            d.pop()
        else:
            d.push()
            d += elm.Capacitor().down().label(f'{res_in:.3f} pF')
            d += elm.Ground()
            d.pop()
    
    # Main matching components
    if rs > rl:
        if xskind_in == 'Capacitor':
            d += elm.Line().right()
            d.push()
            d += elm.Capacitor().down().label(f'{jx1_in:.3f} pF')
            d += elm.Ground()
            d.pop()
            d += elm.Inductor().right().label(f'{jx2_in:.3f} nH')
        else:
            d += elm.Line().right()
            d.push()
            d += elm.Inductor().down().label(f'{jx1_in:.3f} nH')
            d += elm.Ground()
            d.pop()
            d += elm.Capacitor().right().label(f'{jx2_in:.3f} pF')
    else:
        if xskind_in == 'Capacitor':
            d += elm.Capacitor().label(f'{jx1_in:.3f} pF')
            d.push()
            d += elm.Inductor().down().label(f'{jx2_in:.3f} nH')
            d += elm.Ground()
            d.pop()
        else:
            d += elm.Inductor().label(f'{jx1_in:.3f} nH')
            d.push()
            d += elm.Capacitor().down().label(f'{jx2_in:.3f} pF')
            d += elm.Ground()
            d.pop()
    
    # Transistor representation
    d += elm.Resistor().right().label(f'{rin} Ohm')
    d += (Q1 := elm.Bjt())
    d += elm.Ground().at((Q1, 'emitter'))
    d += elm.Resistor().right().at((Q1, 'collector')).label(f'{rout:.3f} Ohm')
    
    # Output matching (similar structure as input)
    # ... (output matching components would go here)
    
    d += elm.Line()
    d += elm.Resistor().down().label(f'{Z0} Ohm')
    d += elm.Ground()
    
    return d