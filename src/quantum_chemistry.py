"""
Quantum Chemistry Module - Core molecular simulation capabilities
Supports both real quantum computations and educational mock simulations
"""

import numpy as np
from typing import Dict, List, Tuple, Optional
from dataclasses import dataclass
from enum import Enum


class BasisSet(Enum):
    """Available basis sets for molecular calculations"""
    STO_3G = "sto-3g"
    MINIMAL = "minimal"
    DOUBLE_ZETA = "6-31g"
    TRIPLE_ZETA = "6-311g"


class MappingType(Enum):
    """Fermion to qubit mapping strategies"""
    JORDAN_WIGNER = "jordan_wigner"
    PARITY = "parity"
    BRAVYI_KITAEV = "bravyi_kitaev"


@dataclass
class Atom:
    """Represents an atom in 3D space"""
    symbol: str
    position: Tuple[float, float, float]  # x, y, z in Angstroms
    
    def distance_to(self, other: 'Atom') -> float:
        """Calculate distance to another atom"""
        dx = self.position[0] - other.position[0]
        dy = self.position[1] - other.position[1]
        dz = self.position[2] - other.position[2]
        return np.sqrt(dx**2 + dy**2 + dz**2)


@dataclass
class Molecule:
    """Molecular structure representation"""
    atoms: List[Atom]
    charge: int = 0
    multiplicity: int = 1
    name: str = "Unknown"
    
    def num_atoms(self) -> int:
        return len(self.atoms)
    
    def num_electrons(self) -> int:
        """Calculate total number of electrons"""
        atomic_numbers = {
            'H': 1, 'He': 2, 'Li': 3, 'Be': 4, 'B': 5, 'C': 6,
            'N': 7, 'O': 8, 'F': 9, 'Ne': 10, 'Na': 11, 'Mg': 12
        }
        total = sum(atomic_numbers.get(atom.symbol, 0) for atom in self.atoms)
        return total - self.charge
    
    def to_xyz_string(self) -> str:
        """Convert to XYZ format"""
        lines = [f"{self.num_atoms()}", f"{self.name}"]
        for atom in self.atoms:
            x, y, z = atom.position
            lines.append(f"{atom.symbol:4s} {x:12.6f} {y:12.6f} {z:12.6f}")
        return "\n".join(lines)
    
    def get_geometry_string(self) -> str:
        """Get geometry in format suitable for quantum chemistry codes"""
        lines = []
        for atom in self.atoms:
            x, y, z = atom.position
            lines.append(f"{atom.symbol} {x:.6f} {y:.6f} {z:.6f}")
        return "; ".join(lines)


class MoleculeBuilder:
    """Factory for creating common molecular structures"""
    
    @staticmethod
    def create_h2(bond_length: float = 0.74) -> Molecule:
        """Create H2 molecule with specified bond length (Angstroms)"""
        atoms = [
            Atom('H', (0.0, 0.0, 0.0)),
            Atom('H', (0.0, 0.0, bond_length))
        ]
        return Molecule(atoms=atoms, charge=0, multiplicity=1, name="H2")
    
    @staticmethod
    def create_lih(bond_length: float = 1.596) -> Molecule:
        """Create LiH molecule"""
        atoms = [
            Atom('Li', (0.0, 0.0, 0.0)),
            Atom('H', (0.0, 0.0, bond_length))
        ]
        return Molecule(atoms=atoms, charge=0, multiplicity=1, name="LiH")
    
    @staticmethod
    def create_h2o(angle_deg: float = 104.5, bond_length: float = 0.96) -> Molecule:
        """Create water molecule"""
        angle_rad = np.deg2rad(angle_deg / 2)
        atoms = [
            Atom('O', (0.0, 0.0, 0.0)),
            Atom('H', (bond_length * np.sin(angle_rad), bond_length * np.cos(angle_rad), 0.0)),
            Atom('H', (-bond_length * np.sin(angle_rad), bond_length * np.cos(angle_rad), 0.0))
        ]
        return Molecule(atoms=atoms, charge=0, multiplicity=1, name="H2O")
    
    @staticmethod
    def create_nh3() -> Molecule:
        """Create ammonia molecule"""
        # Pyramidal structure
        bond_length = 1.012  # Angstroms
        angle = 106.67  # degrees
        h = bond_length * np.cos(np.deg2rad(angle / 2))
        r = bond_length * np.sin(np.deg2rad(angle / 2))
        
        atoms = [
            Atom('N', (0.0, 0.0, 0.0)),
            Atom('H', (0.0, r, h)),
            Atom('H', (r * np.sqrt(3) / 2, -r / 2, h)),
            Atom('H', (-r * np.sqrt(3) / 2, -r / 2, h))
        ]
        return Molecule(atoms=atoms, charge=0, multiplicity=1, name="NH3")
    
    @staticmethod
    def create_ch4() -> Molecule:
        """Create methane molecule (tetrahedral)"""
        bond_length = 1.089  # Angstroms
        # Tetrahedral coordinates
        t = bond_length / np.sqrt(3)
        atoms = [
            Atom('C', (0.0, 0.0, 0.0)),
            Atom('H', (t, t, t)),
            Atom('H', (-t, -t, t)),
            Atom('H', (-t, t, -t)),
            Atom('H', (t, -t, -t))
        ]
        return Molecule(atoms=atoms, charge=0, multiplicity=1, name="CH4")
    
    @staticmethod
    def create_benzene() -> Molecule:
        """Create benzene molecule (planar hexagon)"""
        r = 1.40  # C-C bond length in Angstroms
        atoms = [Atom('C', (0.0, 0.0, 0.0))]  # Center C
        
        # Create hexagonal ring
        for i in range(6):
            angle = i * np.pi / 3
            x = r * np.cos(angle)
            y = r * np.sin(angle)
            atoms.append(Atom('C', (x, y, 0.0)))
            # Add hydrogen
            atoms.append(Atom('H', (1.4 * x, 1.4 * y, 0.0)))
        
        return Molecule(atoms=atoms[1:], charge=0, multiplicity=1, name="C6H6")


@dataclass
class HamiltonianData:
    """Container for molecular Hamiltonian information"""
    num_qubits: int
    num_electrons: int
    nuclear_repulsion: float
    hf_energy: float
    one_body_integrals: Optional[np.ndarray] = None
    two_body_integrals: Optional[np.ndarray] = None
    basis: str = "sto-3g"
    mapping: str = "jordan_wigner"
    
    def get_system_size(self) -> str:
        """Get human-readable system size"""
        return f"{self.num_qubits} qubits, {self.num_electrons} electrons"


class QuantumSimulator:
    """Simulates quantum chemistry calculations"""
    
    def __init__(self, use_mock: bool = True):
        self.use_mock = use_mock
        self.pyscf_available = False
        self.qiskit_available = False
        
        if not use_mock:
            try:
                import pyscf
                self.pyscf_available = True
            except ImportError:
                print("⚠️  PySCF not available, using mock simulations")
            
            try:
                import qiskit
                self.qiskit_available = True
            except ImportError:
                print("⚠️  Qiskit not available, using mock simulations")
    
    def compute_hamiltonian(self, molecule: Molecule, basis: str = "sto-3g") -> HamiltonianData:
        """Compute molecular Hamiltonian"""
        if self.pyscf_available and not self.use_mock:
            return self._compute_real_hamiltonian(molecule, basis)
        else:
            return self._compute_mock_hamiltonian(molecule, basis)
    
    def _compute_mock_hamiltonian(self, molecule: Molecule, basis: str) -> HamiltonianData:
        """Generate realistic mock Hamiltonian for educational purposes"""
        num_electrons = molecule.num_electrons()
        num_orbitals = self._estimate_orbitals(molecule, basis)
        num_qubits = 2 * num_orbitals  # Spin orbitals
        
        # Realistic energy estimates based on molecule type
        energy_per_electron = -0.5  # Rough atomic unit estimate
        nuclear_repulsion = self._estimate_nuclear_repulsion(molecule)
        hf_energy = energy_per_electron * num_electrons + nuclear_repulsion
        
        return HamiltonianData(
            num_qubits=num_qubits,
            num_electrons=num_electrons,
            nuclear_repulsion=nuclear_repulsion,
            hf_energy=hf_energy,
            basis=basis,
            mapping="jordan_wigner"
        )
    
    def _compute_real_hamiltonian(self, molecule: Molecule, basis: str) -> HamiltonianData:
        """Compute actual Hamiltonian using PySCF"""
        try:
            from pyscf import gto, scf
            
            # Build PySCF molecule
            mol = gto.Mole()
            mol.atom = molecule.get_geometry_string()
            mol.basis = basis
            mol.charge = molecule.charge
            mol.spin = molecule.multiplicity - 1
            mol.build()
            
            # Run Hartree-Fock
            mf = scf.RHF(mol)
            hf_energy = mf.kernel()
            
            # Get integrals
            h1e = mol.intor('int1e_kin') + mol.intor('int1e_nuc')
            h2e = mol.intor('int2e')
            
            num_orbitals = mol.nao
            num_qubits = 2 * num_orbitals
            
            return HamiltonianData(
                num_qubits=num_qubits,
                num_electrons=molecule.num_electrons(),
                nuclear_repulsion=mol.energy_nuc(),
                hf_energy=hf_energy,
                one_body_integrals=h1e,
                two_body_integrals=h2e,
                basis=basis
            )
        except Exception as e:
            print(f"⚠️  PySCF calculation failed: {e}")
            return self._compute_mock_hamiltonian(molecule, basis)
    
    def _estimate_orbitals(self, molecule: Molecule, basis: str) -> int:
        """Estimate number of orbitals for basis set"""
        basis_multipliers = {
            "sto-3g": 1,
            "minimal": 1,
            "6-31g": 2,
            "6-311g": 3
        }
        multiplier = basis_multipliers.get(basis, 1)
        
        # Count valence orbitals per atom
        valence_orbitals = {'H': 1, 'He': 1, 'Li': 2, 'C': 4, 'N': 5, 'O': 6}
        total_orbitals = sum(
            valence_orbitals.get(atom.symbol, 4) * multiplier
            for atom in molecule.atoms
        )
        return total_orbitals
    
    def _estimate_nuclear_repulsion(self, molecule: Molecule) -> float:
        """Estimate nuclear repulsion energy"""
        atomic_numbers = {'H': 1, 'Li': 3, 'C': 6, 'N': 7, 'O': 8}
        
        total = 0.0
        for i, atom1 in enumerate(molecule.atoms):
            for atom2 in molecule.atoms[i+1:]:
                z1 = atomic_numbers.get(atom1.symbol, 1)
                z2 = atomic_numbers.get(atom2.symbol, 1)
                distance = atom1.distance_to(atom2)
                if distance > 0.01:  # Avoid division by zero
                    total += (z1 * z2) / distance
        
        return total


class EnergyConverter:
    """Conversion utilities for energy units"""
    
    HARTREE_TO_EV = 27.211386245988
    HARTREE_TO_KJ_MOL = 2625.49962
    HARTREE_TO_KCAL_MOL = 627.509474
    
    @classmethod
    def hartree_to_kj_mol(cls, energy: float) -> float:
        return energy * cls.HARTREE_TO_KJ_MOL
    
    @classmethod
    def hartree_to_kcal_mol(cls, energy: float) -> float:
        return energy * cls.HARTREE_TO_KCAL_MOL
    
    @classmethod
    def hartree_to_ev(cls, energy: float) -> float:
        return energy * cls.HARTREE_TO_EV
    
    @classmethod
    def kj_mol_to_hartree(cls, energy: float) -> float:
        return energy / cls.HARTREE_TO_KJ_MOL
