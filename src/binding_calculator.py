"""
Drug Binding Energy Calculator
Computes binding free energies with thermodynamic corrections
"""

import numpy as np
from typing import Dict, Optional
from dataclasses import dataclass


@dataclass
class BindingEnergyResult:
    """Results from binding energy calculation"""
    delta_e_electronic: float  # Electronic energy change (Hartree)
    delta_g_solvation: float    # Solvation free energy (Hartree)
    entropy_contribution: float  # -TΔS term (Hartree)
    delta_g_binding: float      # Total binding free energy (Hartree)
    
    # Convenience properties in common units
    delta_e_electronic_kj_mol: float
    delta_g_solvation_kj_mol: float
    entropy_contribution_kj_mol: float
    delta_g_binding_kj_mol: float
    delta_g_binding_kcal_mol: float
    
    # Estimated binding affinity
    estimated_kd_nm: float
    estimated_ic50_nm: float


class BindingEnergyCalculator:
    """Calculate protein-ligand binding energies"""
    
    HARTREE_TO_KJ_MOL = 2625.49962
    HARTREE_TO_KCAL_MOL = 627.509474
    R_GAS = 8.314  # J/(mol·K) = 0.008314 kJ/(mol·K)
    
    def __init__(self, temperature: float = 298.15):
        """
        Args:
            temperature: Temperature in Kelvin (default 298.15 K = 25°C)
        """
        self.temperature = temperature
        self.RT = self.R_GAS * temperature / 1000  # in kJ/mol
    
    def calculate_binding_energy(self,
                                 complex_energy: Dict,
                                 protein_energy: Dict,
                                 ligand_energy: Dict,
                                 include_solvation: bool = True,
                                 include_entropy: bool = True,
                                 custom_corrections: Optional[Dict] = None) -> BindingEnergyResult:
        """
        Calculate binding free energy: ΔG = E(complex) - E(protein) - E(ligand) + corrections
        
        Args:
            complex_energy: Energy dict with 'energy' key (in Hartree)
            protein_energy: Energy dict with 'energy' key (in Hartree)
            ligand_energy: Energy dict with 'energy' key (in Hartree)
            include_solvation: Include solvation free energy correction
            include_entropy: Include entropy contribution
            custom_corrections: Dict with optional 'solvation' and 'entropy' values (in Hartree)
        
        Returns:
            BindingEnergyResult with all energy components
        """
        # Electronic energy difference
        delta_e_electronic = (complex_energy['energy'] - 
                            protein_energy['energy'] - 
                            ligand_energy['energy'])
        
        # Solvation correction (empirical)
        if include_solvation:
            if custom_corrections and 'solvation' in custom_corrections:
                delta_g_solvation = custom_corrections['solvation']
            else:
                # Empirical estimate: favorable solvation ~-20 kJ/mol
                delta_g_solvation = -20.0 / self.HARTREE_TO_KJ_MOL
        else:
            delta_g_solvation = 0.0
        
        # Entropy contribution (empirical)
        if include_entropy:
            if custom_corrections and 'entropy' in custom_corrections:
                entropy_contribution = custom_corrections['entropy']
            else:
                # Typical entropy penalty for binding: +50 kJ/mol (unfavorable)
                # This accounts for loss of translational/rotational freedom
                entropy_contribution = 50.0 / self.HARTREE_TO_KJ_MOL
        else:
            entropy_contribution = 0.0
        
        # Total binding free energy
        delta_g_binding = delta_e_electronic + delta_g_solvation + entropy_contribution
        
        # Convert to common units
        delta_e_electronic_kj = delta_e_electronic * self.HARTREE_TO_KJ_MOL
        delta_g_solvation_kj = delta_g_solvation * self.HARTREE_TO_KJ_MOL
        entropy_contribution_kj = entropy_contribution * self.HARTREE_TO_KJ_MOL
        delta_g_binding_kj = delta_g_binding * self.HARTREE_TO_KJ_MOL
        delta_g_binding_kcal = delta_g_binding * self.HARTREE_TO_KCAL_MOL
        
        # Estimate binding affinity (Kd and IC50)
        estimated_kd = self._estimate_kd(delta_g_binding_kj)
        estimated_ic50 = self._estimate_ic50(delta_g_binding_kj)
        
        return BindingEnergyResult(
            delta_e_electronic=delta_e_electronic,
            delta_g_solvation=delta_g_solvation,
            entropy_contribution=entropy_contribution,
            delta_g_binding=delta_g_binding,
            delta_e_electronic_kj_mol=delta_e_electronic_kj,
            delta_g_solvation_kj_mol=delta_g_solvation_kj,
            entropy_contribution_kj_mol=entropy_contribution_kj,
            delta_g_binding_kj_mol=delta_g_binding_kj,
            delta_g_binding_kcal_mol=delta_g_binding_kcal,
            estimated_kd_nm=estimated_kd,
            estimated_ic50_nm=estimated_ic50
        )
    
    def _estimate_kd(self, delta_g_kj_mol: float) -> float:
        """
        Estimate dissociation constant (Kd) from ΔG
        ΔG = RT ln(Kd)  =>  Kd = exp(ΔG/RT)
        
        Returns:
            Kd in nanomolar (nM)
        """
        # Convert to M
        kd_molar = np.exp(delta_g_kj_mol / self.RT)
        # Convert to nM
        kd_nm = kd_molar * 1e9
        return kd_nm
    
    def _estimate_ic50(self, delta_g_kj_mol: float) -> float:
        """
        Estimate IC50 from ΔG (rough approximation: IC50 ≈ 2 * Kd)
        
        Returns:
            IC50 in nanomolar (nM)
        """
        kd_nm = self._estimate_kd(delta_g_kj_mol)
        return 2.0 * kd_nm
    
    def interpret_binding_affinity(self, delta_g_kj_mol: float) -> Dict[str, str]:
        """
        Interpret binding affinity strength
        
        Returns:
            Dict with interpretation and category
        """
        if delta_g_kj_mol < -50:
            return {
                'category': 'Very Strong',
                'description': 'Exceptional binding affinity, potential lead compound',
                'emoji': '⭐⭐⭐'
            }
        elif delta_g_kj_mol < -35:
            return {
                'category': 'Strong',
                'description': 'Good binding affinity, promising candidate',
                'emoji': '⭐⭐'
            }
        elif delta_g_kj_mol < -20:
            return {
                'category': 'Moderate',
                'description': 'Moderate binding, may need optimization',
                'emoji': '⭐'
            }
        elif delta_g_kj_mol < 0:
            return {
                'category': 'Weak',
                'description': 'Weak binding, significant optimization needed',
                'emoji': '⚠️'
            }
        else:
            return {
                'category': 'No Binding',
                'description': 'Unfavorable binding, likely not a binder',
                'emoji': '❌'
            }
    
    def compare_ligands(self, results: Dict[str, BindingEnergyResult]) -> Dict:
        """
        Compare multiple ligands and rank by binding affinity
        
        Args:
            results: Dict mapping ligand name to BindingEnergyResult
        
        Returns:
            Dict with ranking and analysis
        """
        # Sort by ΔG (most negative = strongest binding)
        sorted_ligands = sorted(
            results.items(),
            key=lambda x: x[1].delta_g_binding_kj_mol
        )
        
        best_ligand = sorted_ligands[0]
        worst_ligand = sorted_ligands[-1]
        
        rankings = []
        for rank, (name, result) in enumerate(sorted_ligands, 1):
            interpretation = self.interpret_binding_affinity(result.delta_g_binding_kj_mol)
            rankings.append({
                'rank': rank,
                'name': name,
                'delta_g_kj_mol': result.delta_g_binding_kj_mol,
                'estimated_kd_nm': result.estimated_kd_nm,
                'category': interpretation['category'],
                'emoji': interpretation['emoji']
            })
        
        return {
            'rankings': rankings,
            'best_ligand': {
                'name': best_ligand[0],
                'delta_g': best_ligand[1].delta_g_binding_kj_mol,
                'kd_nm': best_ligand[1].estimated_kd_nm
            },
            'worst_ligand': {
                'name': worst_ligand[0],
                'delta_g': worst_ligand[1].delta_g_binding_kj_mol,
                'kd_nm': worst_ligand[1].estimated_kd_nm
            },
            'delta_g_range': {
                'min': best_ligand[1].delta_g_binding_kj_mol,
                'max': worst_ligand[1].delta_g_binding_kj_mol,
                'span': worst_ligand[1].delta_g_binding_kj_mol - best_ligand[1].delta_g_binding_kj_mol
            }
        }


class ThermodynamicCycle:
    """Calculate binding energies using thermodynamic cycles"""
    
    def __init__(self, calculator: BindingEnergyCalculator):
        self.calculator = calculator
    
    def compute_relative_binding(self,
                                 ligand_a_result: BindingEnergyResult,
                                 ligand_b_result: BindingEnergyResult) -> Dict:
        """
        Calculate relative binding free energy between two ligands
        ΔΔG = ΔG(B) - ΔG(A)
        
        Positive ΔΔG means ligand A binds better than ligand B
        """
        ddg = ligand_b_result.delta_g_binding_kj_mol - ligand_a_result.delta_g_binding_kj_mol
        
        # Calculate fold-change in binding affinity
        fold_change = np.exp(-ddg / self.calculator.RT)
        
        return {
            'delta_delta_g_kj_mol': ddg,
            'fold_change': fold_change,
            'interpretation': (
                f"Ligand A binds {fold_change:.1f}x {'stronger' if ddg > 0 else 'weaker'} than Ligand B"
            )
        }
