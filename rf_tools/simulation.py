# rf_tools/simulation.py
import skrf as rf
import numpy as np
from typing import List, Tuple

class CircuitSimulator:
    """Class for simulating RF circuits."""
    
    def __init__(self, frequency, z0=50):
        """
        Initialize circuit simulator.
        
        Args:
            frequency: Frequency object from skrf
            z0: Characteristic impedance
        """
        self.frequency = frequency
        self.z0 = z0
        self.media = rf.media.DefinedGammaZ0(frequency=frequency, z0=z0)
        self.components = {}
        self.connections = []
    
    def add_component(self, name: str, comp_type: str, value: float):
        """
        Add a component to the circuit.
        
        Args:
            name: Component name
            comp_type: Type ('capacitor', 'inductor', 'resistor', 'port', 'ground')
            value: Component value (Farads, Henries, Ohms for RLC)
        """
        if comp_type == 'capacitor':
            self.components[name] = self.media.capacitor(value, name=name)
        elif comp_type == 'inductor':
            self.components[name] = self.media.inductor(value, name=name)
        elif comp_type == 'resistor':
            self.components[name] = self.media.resistor(value, name=name)
        elif comp_type == 'port':
            self.components[name] = rf.Circuit.Port(
                frequency=self.frequency, name=name, z0=self.z0)
        elif comp_type == 'ground':
            self.components[name] = rf.Circuit.Ground(
                frequency=self.frequency, name=name, z0=self.z0)
    
    def connect(self, connections: List[Tuple]):
        """
        Connect components in the circuit.
        
        Args:
            connections: List of connection tuples (component_name, port)
        """
        self.connections.append(connections)
    
    def simulate(self, power: List[float] = None, phase: List[float] = None):
        """
        Simulate the circuit and return network.
        
        Args:
            power: Optional input powers for each port
            phase: Optional input phases for each port
            
        Returns:
            Network representing the circuit
        """
        circuit = rf.Circuit(self.connections)
        network = circuit.network
        
        if power is not None and phase is not None:
            self.voltages = circuit.voltages_external(power, phase)
            self.currents = circuit.currents_external(power, phase)
        
        return network
    
    def plot_circuit(self):
        """Plot the circuit diagram."""
        circuit = rf.Circuit(self.connections)
        circuit.plot_graph(
            network_labels=True,
            edge_labels=True,
            port_labels=True,
            inter_labels=True,
            network_fontsize=20,
            edges_fontsize=20,
            port_fontsize=20
        )