# Contributing to Quantum Catalyst for Alzheimer's Drug Discovery

Thank you for your interest in contributing to this research project! We welcome contributions from quantum computing researchers, computational chemists, drug discovery scientists, and software developers.

## 🎯 Project Mission

This project aims to demonstrate that **Variational Quantum Eigensolver (VQE) algorithms can achieve chemical accuracy** (<1 kcal/mol) for binding energy calculations in Alzheimer's drug discovery—a precision level classical methods fundamentally cannot reach.

## 🔬 Research Focus Areas

We particularly welcome contributions in:

### 1. **Quantum Algorithms**
- Novel VQE ansätze optimized for Amyloid-beta (Aβ) protein fragments
- Quantum error mitigation techniques (ZNE, CDR, PEC)
- Hybrid quantum-classical optimization strategies
- Qubit-efficient encodings for molecular systems

### 2. **Computational Chemistry**
- Advanced fragmentation methods (FMO, ONIOM, QM/MM)
- Improved solvation models (GBSA, PBSA, explicit water)
- Thermodynamic property calculations (entropy, free energy)
- Benchmarking against experimental Aβ inhibitor data

### 3. **Alzheimer's Research**
- Aβ peptide structure database (PDB conformations)
- Known inhibitor datasets with experimental Kd/IC50 values
- Aggregation mechanism modeling
- Structure-activity relationship (SAR) analysis

### 4. **Software Engineering**
- Performance optimization for large-scale calculations
- Integration with quantum hardware (IBM Quantum, AWS Braket, IonQ)
- Visualization tools for molecular interactions
- CI/CD pipelines and automated testing

## 🚀 Getting Started

### Prerequisites

```bash
# Clone the repository
git clone https://github.com/virtuoso-04/Quantum.git
cd Quantum

# Create virtual environment
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate

# Install dependencies
pip install -r requirements.txt

# Optional: Install quantum chemistry libraries
pip install pyscf qiskit qiskit-nature qiskit-algorithms

# Run the demonstration
python examples/alzheimers_vqe_demo.py
```

### Development Setup

```bash
# Install development dependencies
pip install pytest black flake8 mypy sphinx

# Run tests
pytest tests/

# Format code
black src/ examples/ tests/

# Lint code
flake8 src/ examples/ tests/

# Type checking
mypy src/
```

## 📝 Contribution Guidelines

### Code Contributions

#### 1. Fork and Clone
```bash
# Fork on GitHub, then clone your fork
git clone https://github.com/YOUR_USERNAME/Quantum.git
cd Quantum
git remote add upstream https://github.com/virtuoso-04/Quantum.git
```

#### 2. Create a Feature Branch
```bash
git checkout -b feature/your-feature-name
# Or for bug fixes:
git checkout -b fix/issue-description
```

#### 3. Make Your Changes

**Code Standards:**
- Follow PEP 8 Python style guidelines
- Use type hints for function signatures
- Include docstrings (Google or NumPy style)
- Keep functions focused and modular (<50 lines)
- Add unit tests for new functionality

**Example Code Structure:**
```python
def calculate_binding_energy(
    complex_energy: float,
    protein_energy: float,
    ligand_energy: float,
    temperature: float = 310.15
) -> Dict[str, float]:
    """
    Calculate binding free energy for Aβ inhibitor.
    
    Args:
        complex_energy: VQE energy of Aβ-inhibitor complex (Hartree)
        protein_energy: VQE energy of Aβ fragment alone (Hartree)
        ligand_energy: VQE energy of inhibitor alone (Hartree)
        temperature: Temperature in Kelvin (default: 310.15 K, body temp)
    
    Returns:
        Dictionary containing:
            - delta_g_binding: Free energy of binding (kJ/mol)
            - kd_estimated: Dissociation constant (nM)
            - classification: "strong" | "moderate" | "weak"
    
    References:
        Gilson & Zhou (2007) Chem. Rev. 107:1557-1576
    """
    # Implementation with clear variable names
    delta_e = (complex_energy - protein_energy - ligand_energy) * 2625.5  # Ha to kJ/mol
    # ... rest of implementation
```

#### 4. Add Tests

Create tests in `tests/` directory:

```python
# tests/test_vqe_engine.py
import pytest
from src.vqe_engine import VQEEngine

def test_vqe_convergence_on_h2():
    """Test VQE achieves chemical accuracy on H2 molecule."""
    vqe = VQEEngine(use_mock=True)
    result = vqe.run_vqe(
        hamiltonian_data={'num_qubits': 2, 'energy_exact': -1.137},
        max_iterations=100
    )
    assert result.success
    assert abs(result.energy - (-1.137)) < 0.0016  # Chemical accuracy
```

#### 5. Document Your Changes

- Update relevant docstrings
- Add examples to `examples/` if introducing new features
- Update `README.md` if adding major functionality
- Add references to `LITERATURE_REVIEW.md` if citing new papers

#### 6. Commit and Push

```bash
# Stage your changes
git add .

# Commit with descriptive message
git commit -m "Add fragment-based QM/MM integration for Aβ peptides

- Implement ONIOM-style multilayer fragmentation
- Add tests for 50-atom Aβ fragments
- Benchmark against MP2 reference data
- Achieves <1 kcal/mol accuracy with 15 qubits

Refs: Dapprich et al. (1999) J. Mol. Struct. 461:1-21"

# Push to your fork
git push origin feature/your-feature-name
```

#### 7. Submit Pull Request

- Go to GitHub and create a Pull Request
- Fill out the PR template with:
  - **Description**: What does this PR do?
  - **Motivation**: Why is this change needed?
  - **Testing**: How was it tested?
  - **References**: Cite relevant papers
  - **Checklist**: Tests pass, code formatted, docs updated

### Research Contributions

#### Publishing Results
If your contribution leads to publishable research:

1. **Coordinate with maintainers** before submission
2. **Cite this repository** using `CITATION.cff`
3. **Share preprints** (arXiv) with the community
4. **Add results** to `RESEARCH.md` with reproducible scripts

#### Benchmarking Standards
When adding benchmarks:

- Use standardized test molecules (H₂, LiH, H₂O, NH₃)
- Report metrics: energy error (mHa), convergence iterations, runtime
- Compare against exact solutions or high-level methods (CCSD(T))
- Include hardware details (CPU/GPU, quantum device specs)
- Provide random seeds for reproducibility

## 🧪 Research Protocols

### Experimental Validation

When integrating experimental data:

1. **Source**: Cite original papers (DOI, PubMed ID)
2. **Quality**: Use peer-reviewed data from reputable journals
3. **Format**: Store in `data/experimental/` with metadata JSON
4. **License**: Verify data can be redistributed under MIT license

Example metadata:
```json
{
  "compound": "Aβ42 inhibitor X",
  "kd_nm": 2.3,
  "method": "Isothermal Titration Calorimetry (ITC)",
  "temperature_k": 310.15,
  "reference": "Smith et al. (2023) J. Med. Chem. 66:1234-1245",
  "doi": "10.1021/acs.jmedchem.3b12345"
}
```

### Computational Standards

#### Energy Calculations
- **Basis sets**: Start with STO-3G, validate with 6-31G*, cc-pVDZ
- **Convergence**: Max force <0.001 Ha/Bohr, energy <1e-6 Ha
- **Accuracy target**: <1 kcal/mol (0.0016 Ha) vs reference

#### VQE Optimization
- **Ansatz**: Document depth, gate count, parameter count
- **Optimizer**: Report algorithm, tolerance, max iterations
- **Initial guess**: Specify random seed or initial parameters
- **Shots**: Report number of circuit evaluations

## 📚 Documentation Standards

### Code Documentation

Use Google-style docstrings:

```python
def optimize_vqe_parameters(
    ansatz: str,
    hamiltonian: np.ndarray,
    initial_params: Optional[np.ndarray] = None
) -> VQEResult:
    """
    Optimize VQE parameters using classical optimizer.
    
    This function implements the variational quantum eigensolver (VQE)
    algorithm to find the ground state energy of a molecular Hamiltonian.
    It uses a hybrid quantum-classical approach suitable for NISQ devices.
    
    Args:
        ansatz: Parameterized quantum circuit type. Supported values:
            - "efficient_su2": Hardware-efficient ansatz
            - "uccsd": Unitary coupled-cluster with singles and doubles
        hamiltonian: Molecular Hamiltonian in qubit representation (Pauli basis)
        initial_params: Starting parameters for optimization. If None, uses
            random initialization from uniform distribution [-π, π].
    
    Returns:
        VQEResult object containing:
            - energy: Optimized ground state energy (Hartree)
            - optimal_parameters: Final circuit parameters
            - num_iterations: Number of optimizer steps
            - convergence_history: Energy at each iteration
    
    Raises:
        ValueError: If ansatz type is not recognized
        RuntimeError: If optimization fails to converge
    
    Example:
        >>> from src.vqe_engine import VQEEngine
        >>> vqe = VQEEngine(use_mock=False)
        >>> hamiltonian = get_h2_hamiltonian()
        >>> result = vqe.optimize_vqe_parameters("efficient_su2", hamiltonian)
        >>> print(f"Ground state: {result.energy:.6f} Ha")
    
    References:
        Peruzzo et al. (2014) Nature Communications 5:4213
        Kandala et al. (2017) Nature 549:242-246
    
    Notes:
        - For production use, consider quantum error mitigation (ZNE, CDR)
        - Convergence strongly depends on ansatz choice and initial parameters
        - Typical runtime: 30-100 iterations for small molecules (<10 qubits)
    """
    # Implementation
```

### README Updates

When adding major features, update README sections:
- **Quick Start**: Add usage examples
- **Project Structure**: Document new modules
- **Scientific Background**: Explain new methods
- **References**: Cite new papers

## 🤝 Code of Conduct

This project follows the [Contributor Covenant Code of Conduct](CODE_OF_CONDUCT.md). 

**Key Principles:**
- Respectful and inclusive communication
- Constructive feedback on research and code
- Credit original authors and cite sources properly
- Share knowledge and mentor newcomers
- Prioritize scientific rigor and reproducibility

## 🏆 Recognition

Contributors will be:
- Listed in `README.md` contributors section
- Acknowledged in research papers (as appropriate)
- Invited to co-author publications for substantial contributions
- Credited in `CITATION.cff` for algorithmic contributions

## 📊 Issue Reporting

### Bug Reports

Use the bug report template:

```markdown
**Describe the bug**
VQE fails to converge for H2O molecule with efficient_su2 ansatz

**To Reproduce**
1. Run `python examples/alzheimers_vqe_demo.py`
2. Observe convergence failure at iteration 87
3. Error message: "Optimizer reached max iterations"

**Expected behavior**
VQE should converge within 100 iterations for H2O

**Environment**
- OS: macOS 14.2
- Python: 3.11.5
- NumPy: 1.24.3
- Qiskit: 0.45.1 (if applicable)

**Additional context**
Works fine for H2 and LiH, only fails on H2O
```

### Feature Requests

Propose new features with:
- **Research motivation**: Why is this needed?
- **Scientific background**: Cite relevant papers
- **Implementation ideas**: Suggest approach
- **Expected impact**: How does this advance the research?

## 🔬 Research Collaboration

### Proposing Research Projects

Open an issue with the `research-proposal` label:

```markdown
**Research Question**
Can VQE achieve sub-chemical accuracy (<0.5 kcal/mol) for Aβ fragment 
binding energies using error-mitigated quantum hardware?

**Hypothesis**
Zero-noise extrapolation (ZNE) combined with probabilistic error 
cancellation (PEC) can reduce hardware noise below chemical accuracy threshold.

**Methodology**
1. Implement ZNE and PEC in src/error_mitigation.py
2. Benchmark on IBM Quantum devices (ibm_brisbane, ibm_kyoto)
3. Compare noisy vs mitigated VQE energies for Aβ(16-21) fragment
4. Validate against CCSD(T)/cc-pVDZ reference

**Expected Outcomes**
- Error reduction from ~10 kcal/mol (raw) to <0.5 kcal/mol (mitigated)
- Publication in J. Chem. Theory Comput. or Quantum Sci. Technol.

**Timeline**
3-6 months (hardware access dependent)

**Collaborators Needed**
- Quantum error mitigation expert
- Access to IBM Quantum premium tier
```

### Data Sharing

Share computational results in `data/results/`:

```
data/
├── experimental/          # Experimental Kd, IC50 values
├── reference/             # CCSD(T), MP2 benchmark energies
└── results/
    ├── vqe_h2o_scan.json  # Reproducible VQE results
    └── metadata.yaml       # Computation details
```

## 📖 Additional Resources

### Learning Materials
- **Quantum Computing**: Qiskit Textbook (qiskit.org/learn)
- **Quantum Chemistry**: Szabo & Ostlund "Modern Quantum Chemistry"
- **VQE Tutorials**: IBM Quantum Challenges
- **Drug Discovery**: Computational approaches to Alzheimer's research

### Communication Channels
- **GitHub Issues**: Bug reports, feature requests
- **Discussions**: Research questions, brainstorming
- **Pull Requests**: Code review and feedback

### Related Projects
- PySCF: Python quantum chemistry (github.com/pyscf/pyscf)
- Qiskit Nature: Quantum chemistry module (github.com/Qiskit/qiskit-nature)
- OpenFermion: Quantum simulation of chemistry (github.com/quantumlib/OpenFermion)

## 📜 License

By contributing, you agree that your contributions will be licensed under the MIT License. You retain copyright to your contributions but grant this project a perpetual, worldwide, non-exclusive license to use, modify, and distribute your work.

## ❓ Questions?

- Open a GitHub Discussion for general questions
- Email maintainers for private research collaboration inquiries
- Join our community meetings (schedule in Discussions)

---

**Thank you for contributing to quantum-accelerated Alzheimer's drug discovery!** 🧬⚛️

*Together, we're working to reduce the 15-year drug development timeline and bring hope to millions affected by Alzheimer's disease.*
