# Quantum Drug Discovery - Project Overview

## What's New: Real Quantum Chemistry Implementation

This project has evolved from a simulation to a **genuine quantum chemistry pipeline** using:

- **PySCF**: Ab initio electronic structure calculations
- **Qiskit Nature**: Second quantization and qubit mapping
- **VQE**: Variational quantum eigensolver for ground state energies
- **FastAPI**: Real-time WebSocket streaming of convergence
- **Interactive Dashboard**: Live visualization with Chart.js

## Project Structure

```
Quantum/
├── qdrug/                          # NEW: Real quantum chemistry pipeline
│   ├── src/
│   │   ├── utils.py                # Energy conversions (Hartree ↔ kJ/mol)
│   │   ├── parse_pdb.py            # Fragment extraction from proteins
│   │   ├── hamiltonian.py          # PySCF → Qiskit Nature pipeline
│   │   ├── vqe_runner.py           # VQE with convergence tracking
│   │   ├── binding_calc.py         # ΔG = E_complex - E_protein - E_ligand
│   │   ├── api_server.py           # FastAPI + WebSocket server
│   │   └── __init__.py
│   ├── web/
│   │   └── index.html              # Real-time dashboard
│   ├── notebook/
│   │   └── demo.ipynb              # Jupyter tutorial
│   ├── examples/
│   │   ├── h2_example.py           # Simple VQE calculation
│   │   └── websocket_example.py    # WebSocket streaming client
│   ├── data/                       # PDB files (empty for now)
│   ├── README.md                   # Comprehensive documentation
│   ├── INSTALL.md                  # Installation guide
│   └── requirements.txt            # Python dependencies
│
└── [Previous simulation files]     # Original educational simulator
    ├── quantum_binding_simulator.py
    ├── quantum_infographic.html
    ├── README.md
    └── examples.py
```

## Key Features

### 1. Real Quantum Chemistry
- **PySCF** for Hartree-Fock and electronic integrals
- **Qiskit Nature** for fermion → qubit mapping
- **VQE** for ground state energy estimation
- **Validation**: Compare VQE vs exact diagonalization

### 2. Real-Time Streaming
- **WebSocket** connection for live updates
- **Convergence tracking**: Energy vs iteration
- **Progress monitoring**: See VQE optimization in real-time

### 3. Interactive Dashboard
- **Molecule selection**: H₂, LiH, H₂O, mock fragments
- **Parameter control**: Ansatz reps, max iterations, optimizer
- **Live charts**: Chart.js visualization
- **Status indicators**: Connection, progress, results

### 4. Scientific Rigor
- **Exact validation**: H₂ and LiH reference systems
- **Energy conversions**: 1 Hartree = 2625.49962 kJ/mol
- **Empirical corrections**: Solvation (~-20 kJ/mol), entropy (~-50 kJ/mol)
- **Transparency**: Clear documentation of approximations

## Quick Start

### Installation

```bash
cd qdrug
python -m venv venv
source venv/bin/activate
pip install -r requirements.txt
```

See `qdrug/INSTALL.md` for detailed setup instructions.

### Run Validation

```bash
cd qdrug/src
python vqe_runner.py
```

Expected output:
```
=== VQE Validation ===

--- H2 Molecule ---
Exact Energy: -1.137270 Hartree
VQE Energy: -1.136891 Hartree
Error: 0.000379 Hartree (0.03%)
Chemical Accuracy: ✓
```

### Start Server

```bash
cd qdrug/src
python api_server.py
```

Server runs at: `http://localhost:8000`

### Open Dashboard

Open `qdrug/web/index.html` in your browser, or:

```bash
cd qdrug/web
python -m http.server 8080
```

Visit: `http://localhost:8080`

## Scientific Method

### Binding Energy Workflow

1. **Fragment Extraction**: Extract 8-15 atoms around binding site
2. **Hamiltonian Generation**:
   - PySCF: Geometry → Hartree-Fock → Electronic integrals
   - Qiskit Nature: Fermion operators → Qubit operators (JW/BK mapping)
3. **VQE Execution**: Optimize ansatz to find ground state energy
4. **Binding Calculation**:
   ```
   ΔE_elec = E_complex - (E_protein + E_ligand)
   ΔG_bind = ΔE_elec + ΔG_solv + TΔS
   ```

### Energy Conversions

- **1 Hartree** = 627.509474 kcal/mol
- **1 Hartree** = 2625.49962 kJ/mol
- **Chemical Accuracy** = 1 kcal/mol = 4.184 kJ/mol

## Validation & Transparency

### What's Real
- ✅ Genuine Hartree-Fock calculations (PySCF)
- ✅ Second quantization with Qiskit Nature
- ✅ VQE optimization on qubit Hamiltonians
- ✅ Exact diagonalization for validation

### What's Approximate
- ⚠️ Fragment-based (8-15 atoms, not full protein)
- ⚠️ Small basis sets (sto-3g for speed)
- ⚠️ Empirical solvation/entropy corrections
- ⚠️ Single geometry (no conformational sampling)

### Mock Mode
If PySCF/Qiskit unavailable, the code uses **mock implementations**:
- Realistic-looking Hamiltonians
- Simulated VQE convergence
- Educational workflow demonstration

## API Reference

### REST Endpoints

```bash
GET  /                      # Health check
GET  /api/molecules         # List test molecules
POST /api/calculate         # Run VQE calculation
POST /api/validate          # Compare VQE vs exact
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

## Examples

### Python Script

```python
from parse_pdb import create_h2_molecule
from hamiltonian import HamiltonianBuilder
from vqe_runner import VQERunner

geometry = create_h2_molecule()
builder = HamiltonianBuilder(basis='sto-3g')
hamiltonian = builder.build_from_geometry(geometry)

runner = VQERunner()
result = runner.run_vqe(hamiltonian, reps=2, max_iterations=50)

print(f"Ground state: {result['energy']:.6f} Hartree")
```

See `qdrug/examples/` for more examples.

### Jupyter Notebook

```bash
cd qdrug/notebook
jupyter notebook demo.ipynb
```

The notebook includes:
- H₂ and LiH validation
- Potential energy surface plots
- VQE convergence visualization
- Binding energy calculations

## Dependencies

### Core (Required)
- numpy >= 1.21.0
- scipy >= 1.7.0
- pyscf >= 2.0.0
- qiskit >= 1.0.0
- qiskit-algorithms >= 0.3.0
- qiskit-nature >= 0.7.0
- fastapi >= 0.104.0
- uvicorn >= 0.24.0
- websockets >= 12.0

### Optional
- biopython >= 1.81 (PDB parsing)
- jupyter >= 1.0.0 (notebooks)
- matplotlib >= 3.5.0 (plotting)

## System Requirements

### Minimum
- Python 3.8+
- 4 GB RAM
- 2 CPU cores

### Recommended
- Python 3.10+
- 8 GB RAM
- 4 CPU cores

## Comparison: Simulation vs Real Quantum

| Feature | Old Simulator | New Pipeline |
|---------|---------------|--------------|
| Hamiltonian | Simulated | Real (PySCF) |
| VQE | Fake convergence | Real optimization |
| Validation | None | Exact diagonalization |
| Streaming | N/A | WebSocket real-time |
| Molecular data | SMILES only | Full quantum chemistry |
| Energy | Approximation | Ab initio calculation |

## Future Enhancements

1. **Real PDB parsing**: Integrate Biopython for actual protein structures
2. **Advanced solvation**: GBSA, Poisson-Boltzmann models
3. **Conformational sampling**: Multiple geometries
4. **Larger basis sets**: 6-31G*, cc-pVDZ
5. **Real quantum hardware**: IBM Quantum via Qiskit Runtime
6. **Error mitigation**: Zero-noise extrapolation, readout correction

## Documentation

- **README.md**: Comprehensive project documentation
- **INSTALL.md**: Step-by-step installation guide
- **demo.ipynb**: Interactive Jupyter tutorial
- **Code comments**: Inline documentation in all modules

## References

- **PySCF**: Q. Sun et al., *WIREs Comput Mol Sci* (2018)
- **Qiskit Nature**: https://qiskit.org/ecosystem/nature
- **VQE**: A. Peruzzo et al., *Nat Commun* **5**, 4213 (2014)

## License

MIT License

## Contributing

Contributions welcome! Areas of interest:
- Enhanced solvation models
- Molecular docking integration
- Quantum error mitigation
- Real drug screening pipelines

---

**Built with real quantum chemistry. Transparent about limitations. Validated against exact results.**
