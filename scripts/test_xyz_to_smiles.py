from pathlib import Path
import subprocess

from rdkit import Chem


def xyz_to_smiles(xyz_path: Path) -> str | None:
    """
    Convert a single XYZ file to a canonical SMILES string.
    Uses OpenBabel (obabel) + RDKit.
    """
    # 1) Call OpenBabel: XYZ -> SMILES
    cmd = ["obabel", str(xyz_path), "-osmi"]
    result = subprocess.run(cmd, capture_output=True, text=True)

    if result.returncode != 0:
        print(f"[OpenBabel error] {result.stderr}")
        return None

    out = result.stdout.strip()
    if not out:
        print("[OpenBabel] No output")
        return None

    # OpenBabel prints lines like: "CCO 1CO2"
    first_line = out.splitlines()[0]
    parts = first_line.split()
    if not parts:
        print("[OpenBabel] Could not parse SMILES line:", first_line)
        return None

    obabel_smiles = parts[0]
    print("Raw SMILES from OpenBabel:", obabel_smiles)

    # 2) Clean / canonicalize with RDKit
    mol = Chem.MolFromSmiles(obabel_smiles)
    if mol is None:
        print("[RDKit] Could not parse SMILES")
        return None

    canon = Chem.MolToSmiles(mol)
    return canon


def main():
    # Path to one XYZ file in your repo
    repo_root = Path(__file__).resolve().parents[1]
    xyz_path = repo_root / "data" / "raw" / "xyz" / "1" / "1CO2.xyz"

    print("Testing XYZ -> SMILES for:", xyz_path)

    if not xyz_path.exists():
        print("XYZ file does not exist. Check the path!")
        return

    smiles = xyz_to_smiles(xyz_path)
    if smiles:
        print("✅ Final canonical SMILES:", smiles)
    else:
        print("❌ Failed to convert XYZ to SMILES.")


if __name__ == "__main__":
    main()
