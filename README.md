# Quantum Drug Discovery Pipeline

Real quantum chemistry for protein-ligand binding calculations using VQE.

## Overview

This project implements a genuine quantum chemistry pipeline for drug discovery:
- **PDB Parsing**: Extract active site fragments from protein structures
- **Hamiltonian Generation**: PySCF for electronic integrals → Qiskit Nature for qubit mapping
- **VQE Execution**: Variational quantum eigensolver for ground state energy
- **Binding Calculation**: ΔG = E(complex) - E(protein) - E(ligand) + corrections
- **Real-time Dashboard**: WebSocket streaming of VQE convergence

## Installation

```bash
# Create virtual environment
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate

# Install quantum chemistry packages
pip install numpy scipy
pip install pyscf
pip install qiskit qiskit-algorithms qiskit-nature[pyscf]

# Install web server
pip install fastapi uvicorn websockets

# Optional: Biopython for PDB parsing
pip install biopython
```

## Quick Start

### 1. Run the API Server

```bash
cd qdrug/src
python api_server.py
```

Server runs at `http://localhost:8000`

### 2. Open Dashboard

Open `qdrug/web/index.html` in your browser, or visit:
```
http://localhost:8000
```

### 3. Test Validation Molecules

```bash
# From qdrug/src
python vqe_runner.py
```

This validates VQE accuracy on H₂ and LiH by comparing against exact diagonalization.

## Architecture

```
qdrug/
├── src/
│   ├── utils.py           # Energy conversions (Hartree ↔ kJ/mol)
│   ├── parse_pdb.py       # Fragment extraction from PDB
│   ├── hamiltonian.py     # PySCF → Qiskit Nature pipeline
│   ├── vqe_runner.py      # VQE with convergence tracking
│   ├── binding_calc.py    # ΔG calculation with corrections
│   └── api_server.py      # FastAPI + WebSocket server
├── web/
│   └── index.html         # Interactive dashboard
├── notebook/
│   └── demo.ipynb         # Jupyter examples
└── data/
    └── *.pdb              # Protein structures
```

## Scientific Method

### Energy Calculation

$$E_{\text{complex}}, E_{\text{protein}}, E_{\text{ligand}} \leftarrow \text{VQE}(\hat{H}_{\text{qubit}})$$

$$\Delta E_{\text{elec}} = E_{\text{complex}} - (E_{\text{protein}} + E_{\text{ligand}})$$

$$\Delta G_{\text{bind}} = \Delta E_{\text{elec}} + \Delta G_{\text{solv}} + T\Delta S$$

### Workflow

1. **Fragment Extraction**: Extract 8-15 heavy atoms around ligand from PDB
2. **Hamiltonian Building**: 
   - Geometry → PySCF (RHF, electronic integrals)
   - Second quantization → Fermion operators
   - Jordan-Wigner or Bravyi-Kitaev mapping → Qubit operators
3. **VQE Optimization**:
   - Ansatz: EfficientSU2 or TwoLocal
   - Optimizer: SLSQP, COBYLA, or L-BFGS-B
   - Convergence tracking via callback
4. **Binding Energy**:
   - Convert Hartree → kJ/mol (1 Ha = 2625.49962 kJ/mol)
   - Add empirical solvation (~-20 kJ/mol)
   - Add entropy penalty (~-50 kJ/mol)

## API Endpoints

### REST API

```bash
# Health check
GET /

# List available molecules
GET /api/molecules

# Calculate binding energy
POST /api/calculate
{
  "molecule": "h2",
  "basis": "sto-3g",
  "mapper": "jordan_wigner",
  "ansatz": "efficient_su2",
  "optimizer": "slsqp",
  "reps": 2,
  "max_iterations": 100
}

# Validation (VQE vs exact)
POST /api/validate
```

### WebSocket

```javascript
const ws = new WebSocket('ws://localhost:8000/ws/calculate');

ws.send(JSON.stringify({
  molecule: 'h2',
  ansatz: 'efficient_su2',
  reps: 2,
  max_iterations: 50
}));

ws.onmessage = (event) => {
  const msg = JSON.parse(event.data);
  // msg.type: 'status', 'hamiltonian', 'convergence', 'complete'
};
```

## Validation

The pipeline includes validation against exact diagonalization:

```python
from parse_pdb import create_h2_molecule, create_lih_molecule
from hamiltonian import HamiltonianBuilder
from vqe_runner import VQERunner

builder = HamiltonianBuilder(basis='sto-3g')
runner = VQERunner()

h2_ham = builder.build_from_geometry(create_h2_molecule())
comparison = runner.compare_with_exact(h2_ham)

print(f"VQE Error: {comparison['error']:.6f} Hartree")
print(f"Chemical Accuracy: {comparison['within_chemical_accuracy']}")
```

Expected results:
- **H₂**: VQE error < 0.001 Ha (~2.6 kJ/mol)
- **LiH**: VQE error < 0.002 Ha (~5.2 kJ/mol)

## Limitations & Transparency

This is a **research demonstration**, not production drug discovery:

### Approximations
1. **Fragment-based**: Only 8-15 atoms around binding site (not full protein)
2. **Mock corrections**: Solvation (~-20 kJ/mol) and entropy (~-50 kJ/mol) are empirical estimates
3. **Small basis**: sto-3g for speed; real calculations need larger basis sets
4. **No dynamics**: Single geometry snapshot, ignoring conformational flexibility

### When Libraries Missing
If PySCF/Qiskit unavailable, the code generates **mock Hamiltonians** that simulate realistic VQE behavior. This allows exploration of the dashboard and workflow without full quantum chemistry installation.

## Dashboard Features

- **Real-time VQE convergence**: Live chart updates during optimization
- **Parameter control**: Sliders for ansatz reps, max iterations
- **Molecule selection**: H₂, LiH, H₂O, mock fragment
- **WebSocket status**: Connection indicator
- **Result display**: Final energy in Hartree and kJ/mol
- **Message log**: Timestamped calculation progress

## Example Usage

```python
# Calculate binding energy for water molecule
from parse_pdb import create_water_molecule
from hamiltonian import HamiltonianBuilder
from vqe_runner import VQERunner
from binding_calc import BindingCalculator

# Build Hamiltonian
geometry = create_water_molecule()
builder = HamiltonianBuilder(basis='sto-3g')
h2o_ham = builder.build_from_geometry(geometry)

# Run VQE
runner = VQERunner()
result = runner.run_vqe(h2o_ham, reps=2, max_iterations=100)

print(f"Ground state energy: {result['energy']:.6f} Hartree")
print(f"Converged in {result['optimizer_evals']} iterations")
```

## Converting to Real Drug Discovery

To adapt for actual drug binding:

1. **Replace fragment extraction**: Use real PDB parser with Biopython
2. **Larger basis sets**: Use 6-31G* or cc-pVDZ instead of sto-3g
3. **Advanced solvation**: Integrate GBSA or Poisson-Boltzmann solvers
4. **Conformational sampling**: Run VQE on multiple protein/ligand geometries
5. **Real device**: Submit to IBM Quantum hardware via Qiskit Runtime

## References

- **PySCF**: Q. Sun et al., *WIREs Comput Mol Sci* (2018)
- **Qiskit Nature**: [qiskit.org/ecosystem/nature](https://qiskit.org/ecosystem/nature)
- **VQE**: A. Peruzzo et al., *Nat Commun* **5**, 4213 (2014)
- **Drug binding**: Gilson & Zhou, *Chem Rev* **107**, 1557 (2007)

## License

MIT License - See LICENSE file

## Contributing

Contributions welcome! Focus areas:
- Enhanced solvation models
- Integration with AutoDock or similar docking tools
- Quantum error mitigation for real devices
- Expanded molecular fragment library

---

**Note**: This project demonstrates quantum chemistry principles. For actual drug discovery, consult domain experts and use validated computational chemistry suites.
