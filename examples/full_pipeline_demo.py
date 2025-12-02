"""
Complete Quantum Drug Discovery Pipeline Example
Demonstrates end-to-end workflow with multiple molecules
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

import numpy as np
from quantum_chemistry import MoleculeBuilder, QuantumSimulator, EnergyConverter
from vqe_engine import VQEEngine, QuantumCircuitAnalyzer
from binding_calculator import BindingEnergyCalculator
from visualizer import ResultsVisualizer

try:
    import matplotlib.pyplot as plt
    MATPLOTLIB_AVAILABLE = True
except ImportError:
    print("⚠️  matplotlib not installed. Install with: pip install matplotlib")
    MATPLOTLIB_AVAILABLE = False


def print_section(title):
    """Print formatted section header"""
    print("\n" + "=" * 80)
    print(f"  {title}")
    print("=" * 80 + "\n")


def main():
    """Run complete quantum drug discovery pipeline"""
    
    print("\n🧬 QUANTUM DRUG DISCOVERY PIPELINE")
    print("   Comprehensive Molecular Simulation & Binding Energy Analysis\n")
    
    # Initialize components
    print("Initializing quantum simulation framework...")
    simulator = QuantumSimulator(use_mock=True)
    vqe_engine = VQEEngine(use_mock=True)
    calculator = BindingEnergyCalculator(temperature=298.15)
    
    if MATPLOTLIB_AVAILABLE:
        visualizer = ResultsVisualizer()
    
    print("✓ Framework initialized\n")
    
    # =================================================================
    # PART 1: Validate VQE on Small Molecules
    # =================================================================
    print_section("PART 1: VQE Validation on Test Molecules")
    
    test_molecules = [
        ('H2', MoleculeBuilder.create_h2()),
        ('LiH', MoleculeBuilder.create_lih()),
        ('H2O', MoleculeBuilder.create_h2o())
    ]
    
    validation_results = {}
    
    for name, molecule in test_molecules:
        print(f"🔬 Analyzing {molecule.name}...")
        print(f"   Atoms: {molecule.num_atoms()}, Electrons: {molecule.num_electrons()}")
        
        # Compute Hamiltonian
        ham_data = simulator.compute_hamiltonian(molecule, basis="sto-3g")
        print(f"   System: {ham_data.get_system_size()}")
        print(f"   HF Energy: {ham_data.hf_energy:.6f} Ha")
        
        # Run VQE
        vqe_result = vqe_engine.run_vqe(
            hamiltonian_data=ham_data.__dict__,
            ansatz="efficient_su2",
            optimizer="slsqp",
            reps=2,
            max_iterations=50
        )
        
        print(f"   VQE Energy: {vqe_result.energy:.6f} Ha")
        print(f"   Iterations: {vqe_result.num_iterations}")
        print(f"   Time: {vqe_result.execution_time:.2f}s")
        
        # Validate against exact
        validation = vqe_engine.validate_against_exact(ham_data.__dict__, vqe_result)
        print(f"   Exact Energy: {validation['exact_energy']:.6f} Ha")
        print(f"   Error: {validation['error']:.6f} Ha ({validation['error_percent']:.2f}%)")
        print(f"   Chemical Accuracy: {'✓ YES' if validation['within_chemical_accuracy'] else '✗ NO'}")
        
        validation_results[name] = {
            'molecule': molecule,
            'hamiltonian': ham_data,
            'vqe_result': vqe_result,
            'validation': validation
        }
        print()
    
    # =================================================================
    # PART 2: Potential Energy Surface Scan
    # =================================================================
    print_section("PART 2: Potential Energy Surface (H2)")
    
    bond_lengths = np.linspace(0.5, 2.5, 12)
    pes_energies = []
    
    print(f"Scanning {len(bond_lengths)} bond lengths from {bond_lengths[0]:.2f} to {bond_lengths[-1]:.2f} Å...")
    
    for r in bond_lengths:
        h2 = MoleculeBuilder.create_h2(bond_length=r)
        ham = simulator.compute_hamiltonian(h2, basis="sto-3g")
        vqe_res = vqe_engine.run_vqe(
            hamiltonian_data=ham.__dict__,
            ansatz="efficient_su2",
            reps=2,
            max_iterations=30
        )
        pes_energies.append(vqe_res.energy)
        print(f"   r = {r:.2f} Å: E = {vqe_res.energy:.6f} Ha")
    
    # Find equilibrium
    min_idx = np.argmin(pes_energies)
    eq_length = bond_lengths[min_idx]
    eq_energy = pes_energies[min_idx]
    print(f"\n✓ Equilibrium: r_eq = {eq_length:.2f} Å, E_eq = {eq_energy:.6f} Ha")
    
    # =================================================================
    # PART 3: Drug Binding Energy Calculations
    # =================================================================
    print_section("PART 3: Drug-Target Binding Energy Calculations")
    
    # Simulate drug candidates with mock energies
    drug_systems = {
        'Drug-A (Ibuprofen-like)': {
            'complex': -125.8, 'protein': -95.2, 'ligand': -28.5
        },
        'Drug-B (Aspirin-like)': {
            'complex': -118.3, 'protein': -92.0, 'ligand': -25.1
        },
        'Drug-C (Novel Compound)': {
            'complex': -132.1, 'protein': -97.8, 'ligand': -31.4
        },
        'Drug-D (Reference)': {
            'complex': -115.0, 'protein': -90.5, 'ligand': -23.8
        }
    }
    
    binding_results = {}
    
    for drug_name, energies in drug_systems.items():
        print(f"\n💊 {drug_name}")
        print(f"   E(complex): {energies['complex']:.2f} Ha")
        print(f"   E(protein): {energies['protein']:.2f} Ha")
        print(f"   E(ligand):  {energies['ligand']:.2f} Ha")
        
        result = calculator.calculate_binding_energy(
            complex_energy={'energy': energies['complex']},
            protein_energy={'energy': energies['protein']},
            ligand_energy={'energy': energies['ligand']},
            include_solvation=True,
            include_entropy=True
        )
        
        print(f"\n   Energy Components:")
        print(f"   • ΔE_electronic: {result.delta_e_electronic_kj_mol:+.2f} kJ/mol")
        print(f"   • ΔG_solvation:  {result.delta_g_solvation_kj_mol:+.2f} kJ/mol")
        print(f"   • -TΔS:          {result.entropy_contribution_kj_mol:+.2f} kJ/mol")
        print(f"   ─────────────────────────────────")
        print(f"   • ΔG_binding:    {result.delta_g_binding_kj_mol:+.2f} kJ/mol")
        print(f"                    ({result.delta_g_binding_kcal_mol:+.2f} kcal/mol)")
        
        print(f"\n   Binding Affinity:")
        print(f"   • Est. Kd:  {result.estimated_kd_nm:.2e} nM")
        print(f"   • Est. IC50: {result.estimated_ic50_nm:.2e} nM")
        
        interpretation = calculator.interpret_binding_affinity(result.delta_g_binding_kj_mol)
        print(f"   • Strength: {interpretation['emoji']} {interpretation['category']}")
        print(f"   • Notes: {interpretation['description']}")
        
        binding_results[drug_name] = result
    
    # =================================================================
    # PART 4: Comparative Analysis
    # =================================================================
    print_section("PART 4: Comparative Analysis & Ranking")
    
    comparison = calculator.compare_ligands(binding_results)
    
    print("🏆 Drug Candidate Rankings:\n")
    for item in comparison['rankings']:
        print(f"   #{item['rank']} {item['name']}")
        print(f"      ΔG = {item['delta_g_kj_mol']:+.2f} kJ/mol | "
              f"Kd ≈ {item['estimated_kd_nm']:.2e} nM | "
              f"{item['emoji']} {item['category']}")
        print()
    
    print(f"Best Candidate: {comparison['best_ligand']['name']}")
    print(f"   ΔG = {comparison['best_ligand']['delta_g']:.2f} kJ/mol")
    print(f"   Kd ≈ {comparison['best_ligand']['kd_nm']:.2e} nM")
    
    # =================================================================
    # PART 5: Generate Visualizations
    # =================================================================
    if MATPLOTLIB_AVAILABLE:
        print_section("PART 5: Generating Visualizations")
        
        # PES curve
        visualizer.plot_pes_curve(
            bond_lengths, pes_energies,
            title="H₂ Potential Energy Surface",
            filename="examples/outputs/pes_curve.png"
        )
        print("✓ Generated: pes_curve.png")
        
        # VQE convergence
        h2_vqe = validation_results['H2']['vqe_result']
        iterations = [h['iteration'] for h in h2_vqe.convergence_history]
        energies = [h['energy'] for h in h2_vqe.convergence_history]
        visualizer.plot_convergence(
            iterations, energies,
            title="VQE Convergence (H₂)",
            filename="examples/outputs/vqe_convergence.png"
        )
        print("✓ Generated: vqe_convergence.png")
        
        # Binding comparison
        visualizer.plot_binding_comparison(
            binding_results,
            filename="examples/outputs/binding_comparison.png"
        )
        print("✓ Generated: binding_comparison.png")
        
        print("\n📁 All visualizations saved to examples/outputs/")
    
    # =================================================================
    # Summary
    # =================================================================
    print_section("SUMMARY")
    
    print("✅ Pipeline Execution Complete!\n")
    print(f"📊 Results:")
    print(f"   • Validated {len(test_molecules)} molecules with VQE")
    print(f"   • Scanned {len(bond_lengths)} points on H₂ PES")
    print(f"   • Evaluated {len(drug_systems)} drug candidates")
    print(f"   • Best candidate: {comparison['best_ligand']['name']}")
    print(f"   • ΔG range: {comparison['delta_g_range']['min']:.1f} to {comparison['delta_g_range']['max']:.1f} kJ/mol")
    
    if MATPLOTLIB_AVAILABLE:
        print(f"\n📈 Generated 3 visualization plots")
    
    print("\n🎯 Next Steps:")
    print("   1. Refine candidate structures with molecular dynamics")
    print("   2. Validate with larger basis sets (6-31G*, cc-pVDZ)")
    print("   3. Run on real quantum hardware for selected candidates")
    print("   4. Perform experimental validation of top candidates")
    
    print("\n" + "=" * 80)
    print("  Thank you for using Quantum Drug Discovery Pipeline!")
    print("=" * 80 + "\n")


if __name__ == "__main__":
    main()
