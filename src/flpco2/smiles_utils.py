"""
SMILES and molecular identifier generation utilities.

Generates SMILES, InChI, and InChIKey from XYZ files using OpenBabel and RDKit.
Includes round-trip validation to ensure structural fidelity.
"""

import subprocess
import tempfile
from pathlib import Path
from typing import Dict, Optional, Tuple
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors
from rdkit.Chem.MolStandardize import rdMolStandardize


def xyz_to_smiles_openbabel(xyz_path: Path) -> Optional[str]:
    """
    Convert XYZ file to SMILES using OpenBabel.
    
    Args:
        xyz_path: Path to XYZ file
        
    Returns:
        SMILES string or None if conversion fails
    """
    try:
        result = subprocess.run(
            ["obabel", str(xyz_path), "-osmi"],
            capture_output=True,
            text=True,
            timeout=10
        )
        
        if result.returncode == 0 and result.stdout.strip():
            # Extract SMILES (first field before any whitespace/newline)
            smiles = result.stdout.strip().split()[0]
            return smiles
        return None
    except (subprocess.TimeoutExpired, FileNotFoundError, IndexError):
        return None


def smiles_to_inchi(smiles: str) -> Tuple[Optional[str], Optional[str]]:
    """
    Convert SMILES to InChI and InChIKey using RDKit.
    
    Args:
        smiles: SMILES string
        
    Returns:
        Tuple of (InChI, InChIKey) or (None, None) if conversion fails
    """
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None, None
        
        inchi = Chem.MolToInchi(mol)
        inchikey = Chem.MolToInchiKey(mol)
        
        return inchi, inchikey
    except:
        return None, None


def sanitize_smiles_rdkit(smiles: str) -> Optional[str]:
    """
    Sanitize and canonicalize SMILES using RDKit.
    
    Args:
        smiles: Input SMILES string
        
    Returns:
        Canonical SMILES or None if sanitization fails
    """
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
        
        # Sanitize
        Chem.SanitizeMol(mol)
        
        # Return canonical SMILES
        return Chem.MolToSmiles(mol, canonical=True)
    except:
        return None


def read_xyz_coordinates(xyz_path: Path) -> Tuple[list, np.ndarray]:
    """
    Read atomic symbols and coordinates from XYZ file.
    
    Args:
        xyz_path: Path to XYZ file
        
    Returns:
        Tuple of (atomic_symbols list, coordinates array Nx3)
    """
    with open(xyz_path, 'r') as f:
        lines = f.readlines()
    
    n_atoms = int(lines[0].strip())
    
    symbols = []
    coords = []
    
    for line in lines[2:2+n_atoms]:
        parts = line.strip().split()
        if len(parts) >= 4:
            symbols.append(parts[0])
            coords.append([float(parts[1]), float(parts[2]), float(parts[3])])
    
    return symbols, np.array(coords)


def smiles_to_xyz_coordinates(smiles: str) -> Tuple[Optional[list], Optional[np.ndarray]]:
    """
    Generate 3D coordinates from SMILES.
    
    Args:
        smiles: SMILES string
        
    Returns:
        Tuple of (atomic_symbols, coordinates) or (None, None) if fails
    """
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None, None
        
        # Add hydrogens
        mol = Chem.AddHs(mol)
        
        # Generate 3D coordinates
        if AllChem.EmbedMolecule(mol, randomSeed=42) != 0:
            return None, None
        
        # Optimize geometry
        AllChem.MMFFOptimizeMolecule(mol, maxIters=200)
        
        # Extract coordinates
        conf = mol.GetConformer()
        symbols = [atom.GetSymbol() for atom in mol.GetAtoms()]
        coords = np.array([list(conf.GetAtomPosition(i)) for i in range(mol.GetNumAtoms())])
        
        return symbols, coords
    except:
        return None, None


def compute_rmsd(coords1: np.ndarray, coords2: np.ndarray) -> float:
    """
    Compute RMSD between two sets of coordinates after optimal alignment.
    
    Uses Kabsch algorithm for alignment.
    
    Args:
        coords1: First coordinate set (Nx3)
        coords2: Second coordinate set (Nx3)
        
    Returns:
        RMSD value in Angstroms
    """
    if coords1.shape != coords2.shape:
        return float('inf')
    
    # Center both coordinate sets
    coords1_centered = coords1 - coords1.mean(axis=0)
    coords2_centered = coords2 - coords2.mean(axis=0)
    
    # Kabsch algorithm
    H = coords1_centered.T @ coords2_centered
    U, S, Vt = np.linalg.svd(H)
    R = Vt.T @ U.T
    
    # Ensure right-handed coordinate system
    if np.linalg.det(R) < 0:
        Vt[-1, :] *= -1
        R = Vt.T @ U.T
    
    # Apply rotation
    coords2_aligned = coords2_centered @ R
    
    # Compute RMSD
    diff = coords1_centered - coords2_aligned
    rmsd = np.sqrt((diff ** 2).sum() / len(coords1))
    
    return rmsd


def validate_smiles_roundtrip(
    xyz_path: Path,
    smiles: str,
    rmsd_threshold: float = 1.0
) -> Dict[str, any]:
    """
    Validate SMILES by round-trip conversion and RMSD check.
    
    Args:
        xyz_path: Original XYZ file
        smiles: SMILES to validate
        rmsd_threshold: Maximum acceptable RMSD in Angstroms
        
    Returns:
        Dictionary with validation results
    """
    result = {
        "passed": False,
        "rmsd": None,
        "error": None,
    }
    
    try:
        # Read original coordinates
        symbols_orig, coords_orig = read_xyz_coordinates(xyz_path)
        
        # Generate coordinates from SMILES
        symbols_gen, coords_gen = smiles_to_xyz_coordinates(smiles)
        
        if symbols_gen is None or coords_gen is None:
            result["error"] = "Failed to generate 3D coordinates from SMILES"
            return result
        
        # Check atom count match
        if len(symbols_orig) != len(symbols_gen):
            result["error"] = f"Atom count mismatch: {len(symbols_orig)} vs {len(symbols_gen)}"
            return result
        
        # Check element match (order may differ)
        if sorted(symbols_orig) != sorted(symbols_gen):
            result["error"] = "Element composition mismatch"
            return result
        
        # Compute RMSD
        rmsd = compute_rmsd(coords_orig, coords_gen)
        result["rmsd"] = float(rmsd)
        result["passed"] = rmsd <= rmsd_threshold
        
        if not result["passed"]:
            result["error"] = f"RMSD {rmsd:.2f} Å exceeds threshold {rmsd_threshold} Å"
        
    except Exception as e:
        result["error"] = f"Validation error: {str(e)}"
    
    return result


def infer_smiles_from_xyz(
    xyz_path: Path,
    validate: bool = True,
    rmsd_threshold: float = 1.0
) -> Dict[str, any]:
    """
    Infer molecular identifiers from XYZ file using multiple strategies.
    
    Args:
        xyz_path: Path to XYZ file
        validate: Whether to perform round-trip validation
        rmsd_threshold: RMSD threshold for validation (Angstroms)
        
    Returns:
        Dictionary with SMILES, InChI, InChIKey, and validation results
    """
    result = {
        "openbabel_smiles": None,
        "rdkit_smiles": None,
        "canonical_smiles": None,
        "inchi": None,
        "inchikey": None,
        "validation": {
            "attempted": False,
            "passed": False,
            "rmsd": None,
            "error": None,
        },
    }
    
    # Get SMILES from OpenBabel
    ob_smiles = xyz_to_smiles_openbabel(xyz_path)
    result["openbabel_smiles"] = ob_smiles
    
    if ob_smiles is None:
        result["validation"]["error"] = "OpenBabel conversion failed"
        return result
    
    # Sanitize with RDKit
    rdkit_smiles = sanitize_smiles_rdkit(ob_smiles)
    result["rdkit_smiles"] = rdkit_smiles
    result["canonical_smiles"] = rdkit_smiles  # RDKit returns canonical
    
    if rdkit_smiles is None:
        result["validation"]["error"] = "RDKit sanitization failed"
        return result
    
    # Generate InChI and InChIKey
    inchi, inchikey = smiles_to_inchi(rdkit_smiles)
    result["inchi"] = inchi
    result["inchikey"] = inchikey
    
    # Validate if requested
    if validate:
        result["validation"]["attempted"] = True
        validation_result = validate_smiles_roundtrip(xyz_path, rdkit_smiles, rmsd_threshold)
        result["validation"].update(validation_result)
    
    return result


if __name__ == "__main__":
    # Test the SMILES inference
    import sys
    from rich.console import Console
    
    console = Console()
    
    if len(sys.argv) > 1:
        xyz_file = Path(sys.argv[1])
        if xyz_file.exists():
            console.print(f"[bold]Processing:[/bold] {xyz_file}")
            result = infer_smiles_from_xyz(xyz_file, validate=True)
            
            console.print(f"\n[cyan]OpenBabel SMILES:[/cyan] {result['openbabel_smiles']}")
            console.print(f"[cyan]RDKit SMILES:[/cyan] {result['rdkit_smiles']}")
            console.print(f"[cyan]InChI:[/cyan] {result['inchi']}")
            console.print(f"[cyan]InChIKey:[/cyan] {result['inchikey']}")
            
            if result["validation"]["attempted"]:
                status = "✓ PASSED" if result["validation"]["passed"] else "✗ FAILED"
                color = "green" if result["validation"]["passed"] else "red"
                console.print(f"\n[{color}]Validation:[/{color}] {status}")
                if result["validation"]["rmsd"] is not None:
                    console.print(f"RMSD: {result['validation']['rmsd']:.3f} Å")
                if result["validation"]["error"]:
                    console.print(f"Error: {result['validation']['error']}")
        else:
            console.print(f"[red]File not found:[/red] {xyz_file}")
    else:
        console.print("Usage: python smiles_utils.py <xyz_file>")

