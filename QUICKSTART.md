# 🚀 Quick Start Guide

Welcome to the Quantum Drug Discovery Platform! Here's how to get started quickly.

## ⚡ 5-Minute Setup

### 1. Install Dependencies
```bash
pip install -r requirements.txt
```

### 2. Validate Installation
```bash
python examples/validate_installation.py
```

### 3. Run Complete Demo
```bash
python examples/full_pipeline_demo.py
```

### 4. Explore Interactive Tutorial
```bash
jupyter notebook tutorials/quantum_drug_discovery_tutorial.ipynb
```

## 📚 What to Read First

1. **[ENHANCEMENT_SUMMARY.md](ENHANCEMENT_SUMMARY.md)** - See what's new and improved
2. **[README.md](README.md)** - Complete documentation and API reference
3. **[CONTRIBUTING.md](CONTRIBUTING.md)** - How to contribute

## 🎯 Common Tasks

### Create a Molecule
```python
from quantum_chemistry import MoleculeBuilder

h2 = MoleculeBuilder.create_h2(bond_length=0.74)
h2o = MoleculeBuilder.create_h2o()
```

### Run VQE Calculation
```python
from quantum_chemistry import QuantumSimulator
from vqe_engine import VQEEngine

simulator = QuantumSimulator()
hamiltonian = simulator.compute_hamiltonian(h2o, basis="sto-3g")

vqe = VQEEngine()
result = vqe.run_vqe(hamiltonian.__dict__, ansatz="efficient_su2")
print(f"Energy: {result.energy:.6f} Ha")
```

### Calculate Binding Energy
```python
from binding_calculator import BindingEnergyCalculator

calculator = BindingEnergyCalculator(temperature=298.15)
result = calculator.calculate_binding_energy(
    complex_energy={'energy': -125.8},
    protein_energy={'energy': -95.2},
    ligand_energy={'energy': -28.5}
)
print(f"ΔG: {result.delta_g_binding_kj_mol:.2f} kJ/mol")
print(f"Kd: {result.estimated_kd_nm:.2e} nM")
```

## 🔍 Explore Examples

- **`examples/full_pipeline_demo.py`** - Complete workflow (VQE validation, PES scan, binding calculations)
- **`examples/validate_installation.py`** - Test suite (6 validation tests)
- **`tutorials/quantum_drug_discovery_tutorial.ipynb`** - Interactive learning (25+ cells)

## 🆘 Need Help?

- Read the full [README.md](README.md) for detailed documentation
- Check [ENHANCEMENT_SUMMARY.md](ENHANCEMENT_SUMMARY.md) for feature overview
- Open an issue on GitHub
- Email: your.email@example.com

## 🎓 Learning Path

1. **Beginner**: Start with Jupyter tutorial
2. **Intermediate**: Run full pipeline demo
3. **Advanced**: Explore source code in `src/`
4. **Expert**: Contribute new features!

---

**Happy Quantum Computing! ⚛️**
