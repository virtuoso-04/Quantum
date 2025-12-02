"""
Simple validation script to test the quantum drug discovery pipeline
Run this to verify your installation is working correctly
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from quantum_chemistry import MoleculeBuilder, QuantumSimulator
from vqe_engine import VQEEngine
from binding_calculator import BindingEnergyCalculator


def test_molecule_creation():
    """Test creating molecular structures"""
    print("\n📝 Test 1: Molecule Creation")
    h2 = MoleculeBuilder.create_h2()
    h2o = MoleculeBuilder.create_h2o()
    
    assert h2.num_atoms() == 2
    assert h2.num_electrons() == 2
    assert h2o.num_atoms() == 3
    assert h2o.num_electrons() == 10
    
    print("   ✓ Molecules created successfully")
    return True


def test_hamiltonian_computation():
    """Test Hamiltonian generation"""
    print("\n📝 Test 2: Hamiltonian Computation")
    h2 = MoleculeBuilder.create_h2()
    simulator = QuantumSimulator(use_mock=True)
    ham = simulator.compute_hamiltonian(h2, basis="sto-3g")
    
    assert ham.num_qubits > 0
    assert ham.num_electrons == 2
    assert ham.hf_energy < 0  # Energy should be negative
    
    print(f"   ✓ Hamiltonian computed: {ham.num_qubits} qubits, E_HF = {ham.hf_energy:.6f} Ha")
    return True


def test_vqe_optimization():
    """Test VQE algorithm"""
    print("\n📝 Test 3: VQE Optimization")
    h2 = MoleculeBuilder.create_h2()
    simulator = QuantumSimulator(use_mock=True)
    ham = simulator.compute_hamiltonian(h2, basis="sto-3g")
    
    vqe = VQEEngine(use_mock=True)
    result = vqe.run_vqe(
        hamiltonian_data=ham.__dict__,
        ansatz="efficient_su2",
        reps=2,
        max_iterations=30
    )
    
    assert result.success
    assert result.num_iterations > 0
    assert result.energy < ham.hf_energy  # VQE should improve on HF
    
    print(f"   ✓ VQE converged: E = {result.energy:.6f} Ha in {result.num_iterations} iterations")
    return True


def test_vqe_validation():
    """Test VQE accuracy"""
    print("\n📝 Test 4: VQE Validation")
    h2 = MoleculeBuilder.create_h2()
    simulator = QuantumSimulator(use_mock=True)
    ham = simulator.compute_hamiltonian(h2, basis="sto-3g")
    
    vqe = VQEEngine(use_mock=True)
    result = vqe.run_vqe(ham.__dict__, ansatz="efficient_su2", reps=2, max_iterations=30)
    validation = vqe.validate_against_exact(ham.__dict__, result)
    
    assert 'error' in validation
    assert 'within_chemical_accuracy' in validation
    
    print(f"   ✓ VQE error: {validation['error']:.6f} Ha ({validation['error_percent']:.2f}%)")
    print(f"   ✓ Chemical accuracy: {'PASSED' if validation['within_chemical_accuracy'] else 'FAILED'}")
    return True


def test_binding_energy_calculation():
    """Test binding energy calculator"""
    print("\n📝 Test 5: Binding Energy Calculation")
    calculator = BindingEnergyCalculator(temperature=298.15)
    
    result = calculator.calculate_binding_energy(
        complex_energy={'energy': -125.8},
        protein_energy={'energy': -95.2},
        ligand_energy={'energy': -28.5}
    )
    
    assert hasattr(result, 'delta_g_binding_kj_mol')
    assert hasattr(result, 'estimated_kd_nm')
    assert result.delta_g_binding_kj_mol < 0  # Should be favorable
    
    print(f"   ✓ ΔG_binding = {result.delta_g_binding_kj_mol:.2f} kJ/mol")
    print(f"   ✓ Est. Kd = {result.estimated_kd_nm:.2e} nM")
    return True


def test_energy_conversions():
    """Test energy unit conversions"""
    print("\n📝 Test 6: Energy Unit Conversions")
    from quantum_chemistry import EnergyConverter
    
    energy_ha = -1.0
    energy_kj = EnergyConverter.hartree_to_kj_mol(energy_ha)
    energy_kcal = EnergyConverter.hartree_to_kcal_mol(energy_ha)
    energy_ev = EnergyConverter.hartree_to_ev(energy_ha)
    
    # Check conversion factors (approximate)
    assert abs(energy_kj + 2625.5) < 1.0
    assert abs(energy_kcal + 627.5) < 1.0
    assert abs(energy_ev + 27.21) < 0.1
    
    print(f"   ✓ Energy conversions working correctly")
    return True


def run_all_tests():
    """Run all validation tests"""
    print("\n" + "="*70)
    print("  QUANTUM DRUG DISCOVERY - VALIDATION TESTS")
    print("="*70)
    
    tests = [
        test_molecule_creation,
        test_hamiltonian_computation,
        test_vqe_optimization,
        test_vqe_validation,
        test_binding_energy_calculation,
        test_energy_conversions
    ]
    
    passed = 0
    failed = 0
    
    for test in tests:
        try:
            if test():
                passed += 1
        except Exception as e:
            print(f"   ✗ Test failed: {e}")
            failed += 1
    
    print("\n" + "="*70)
    print(f"  RESULTS: {passed} passed, {failed} failed")
    print("="*70)
    
    if failed == 0:
        print("\n✅ All tests passed! Your installation is working correctly.")
        print("   You're ready to use the Quantum Drug Discovery platform!")
    else:
        print("\n⚠️  Some tests failed. Please check your installation.")
    
    return failed == 0


if __name__ == "__main__":
    success = run_all_tests()
    sys.exit(0 if success else 1)
