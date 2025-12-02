"""
Results Visualization Module
Creates publication-quality plots for quantum chemistry results
"""

import numpy as np
from typing import List, Dict, Optional
import os


class ResultsVisualizer:
    """Create visualizations for quantum chemistry results"""
    
    def __init__(self):
        try:
            import matplotlib.pyplot as plt
            from matplotlib import rcParams
            self.plt = plt
            self.rcParams = rcParams
            self.available = True
            
            # Configure matplotlib for publication-quality plots
            rcParams['figure.dpi'] = 150
            rcParams['font.size'] = 11
            rcParams['axes.grid'] = True
            rcParams['grid.alpha'] = 0.3
            rcParams['font.family'] = 'sans-serif'
        except ImportError:
            self.available = False
            print("⚠️  matplotlib not available")
    
    def plot_pes_curve(self, 
                      bond_lengths: List[float],
                      energies: List[float],
                      title: str = "Potential Energy Surface",
                      filename: Optional[str] = None):
        """Plot potential energy surface"""
        if not self.available:
            return
        
        plt = self.plt
        fig, ax = plt.subplots(figsize=(10, 6))
        
        ax.plot(bond_lengths, energies, 'o-', 
               linewidth=2.5, markersize=8, color='#0066CC', 
               label='VQE Energy')
        
        # Mark minimum
        min_idx = np.argmin(energies)
        ax.plot(bond_lengths[min_idx], energies[min_idx], 
               'r*', markersize=20, label=f'Minimum: {bond_lengths[min_idx]:.2f} Å')
        
        ax.set_xlabel('Bond Length (Å)', fontsize=13, fontweight='bold')
        ax.set_ylabel('Energy (Hartree)', fontsize=13, fontweight='bold')
        ax.set_title(title, fontsize=15, fontweight='bold')
        ax.legend(fontsize=11, framealpha=0.9)
        ax.grid(alpha=0.3)
        
        plt.tight_layout()
        
        if filename:
            os.makedirs(os.path.dirname(filename), exist_ok=True)
            plt.savefig(filename, dpi=300, bbox_inches='tight')
        else:
            plt.show()
        
        plt.close()
    
    def plot_convergence(self,
                        iterations: List[int],
                        energies: List[float],
                        title: str = "VQE Convergence",
                        filename: Optional[str] = None):
        """Plot VQE convergence"""
        if not self.available:
            return
        
        plt = self.plt
        fig, ax = plt.subplots(figsize=(10, 6))
        
        ax.plot(iterations, energies, 'o-', 
               linewidth=2, markersize=5, color='#0066CC')
        
        final_energy = energies[-1]
        ax.axhline(y=final_energy, color='#FF4444', linestyle='--', 
                  linewidth=2, alpha=0.7,
                  label=f'Final: {final_energy:.6f} Ha')
        
        ax.set_xlabel('Iteration', fontsize=13, fontweight='bold')
        ax.set_ylabel('Energy (Hartree)', fontsize=13, fontweight='bold')
        ax.set_title(title, fontsize=15, fontweight='bold')
        ax.legend(fontsize=11, framealpha=0.9)
        ax.grid(alpha=0.3)
        
        plt.tight_layout()
        
        if filename:
            os.makedirs(os.path.dirname(filename), exist_ok=True)
            plt.savefig(filename, dpi=300, bbox_inches='tight')
        else:
            plt.show()
        
        plt.close()
    
    def plot_binding_comparison(self,
                               results: Dict,
                               filename: Optional[str] = None):
        """Plot binding energy comparison"""
        if not self.available:
            return
        
        plt = self.plt
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 7))
        
        names = list(results.keys())
        delta_e_elec = [results[n].delta_e_electronic_kj_mol for n in names]
        delta_g_solv = [results[n].delta_g_solvation_kj_mol for n in names]
        t_delta_s = [results[n].entropy_contribution_kj_mol for n in names]
        
        x = np.arange(len(names))
        width = 0.6
        
        # Stacked components
        ax1.bar(x, delta_e_elec, width, label='ΔE_elec', color='#0066CC')
        ax1.bar(x, delta_g_solv, width, bottom=delta_e_elec, 
               label='ΔG_solv', color='#4DA8FF')
        ax1.bar(x, t_delta_s, width, 
               bottom=np.array(delta_e_elec)+np.array(delta_g_solv), 
               label='-TΔS', color='#CCE5FF')
        
        ax1.set_xlabel('Drug Candidate', fontsize=13, fontweight='bold')
        ax1.set_ylabel('Energy (kJ/mol)', fontsize=13, fontweight='bold')
        ax1.set_title('ΔG Components Breakdown', fontsize=15, fontweight='bold')
        ax1.set_xticks(x)
        ax1.set_xticklabels([n.split()[0] for n in names], rotation=20, ha='right')
        ax1.legend(fontsize=11, framealpha=0.9)
        ax1.axhline(y=0, color='black', linestyle='-', linewidth=1)
        ax1.grid(alpha=0.3)
        
        # Total ΔG
        delta_g_total = [results[n].delta_g_binding_kj_mol for n in names]
        colors = ['#00AA00' if dg < -35 else '#FFA500' if dg < -25 else '#FF4444' 
                 for dg in delta_g_total]
        
        bars = ax2.bar(x, delta_g_total, width, color=colors, 
                      edgecolor='black', linewidth=1.5, alpha=0.85)
        
        ax2.set_xlabel('Drug Candidate', fontsize=13, fontweight='bold')
        ax2.set_ylabel('ΔG_binding (kJ/mol)', fontsize=13, fontweight='bold')
        ax2.set_title('Total Binding Free Energy', fontsize=15, fontweight='bold')
        ax2.set_xticks(x)
        ax2.set_xticklabels([n.split()[0] for n in names], rotation=20, ha='right')
        ax2.grid(alpha=0.3)
        
        # Value labels
        for bar, val in zip(bars, delta_g_total):
            height = bar.get_height()
            ax2.text(bar.get_x() + bar.get_width()/2., height,
                    f'{val:.1f}', ha='center', 
                    va='bottom' if val < 0 else 'top', 
                    fontweight='bold', fontsize=10)
        
        plt.tight_layout()
        
        if filename:
            os.makedirs(os.path.dirname(filename), exist_ok=True)
            plt.savefig(filename, dpi=300, bbox_inches='tight')
        else:
            plt.show()
        
        plt.close()
    
    def plot_molecule_3d(self, molecule, filename: Optional[str] = None):
        """Plot 3D molecular structure"""
        if not self.available:
            return
        
        try:
            from mpl_toolkits.mplot3d import Axes3D
            plt = self.plt
            
            fig = plt.figure(figsize=(10, 10))
            ax = fig.add_subplot(111, projection='3d')
            
            # Color scheme for atoms
            atom_colors = {
                'H': '#FFFFFF', 'C': '#808080', 'N': '#0000FF',
                'O': '#FF0000', 'Li': '#800080', 'F': '#90EE90'
            }
            
            atom_sizes = {
                'H': 100, 'C': 200, 'N': 200,
                'O': 200, 'Li': 250, 'F': 180
            }
            
            for atom in molecule.atoms:
                x, y, z = atom.position
                color = atom_colors.get(atom.symbol, '#CCCCCC')
                size = atom_sizes.get(atom.symbol, 150)
                ax.scatter(x, y, z, c=color, s=size, 
                          edgecolors='black', linewidth=2,
                          label=atom.symbol if atom.symbol not in [a.symbol for a in molecule.atoms[:molecule.atoms.index(atom)]] else "")
            
            ax.set_xlabel('X (Å)', fontsize=11, fontweight='bold')
            ax.set_ylabel('Y (Å)', fontsize=11, fontweight='bold')
            ax.set_zlabel('Z (Å)', fontsize=11, fontweight='bold')
            ax.set_title(f'3D Structure: {molecule.name}', fontsize=13, fontweight='bold')
            ax.legend(fontsize=10)
            
            plt.tight_layout()
            
            if filename:
                os.makedirs(os.path.dirname(filename), exist_ok=True)
                plt.savefig(filename, dpi=300, bbox_inches='tight')
            else:
                plt.show()
            
            plt.close()
        except ImportError:
            print("⚠️  3D plotting requires mpl_toolkits")
