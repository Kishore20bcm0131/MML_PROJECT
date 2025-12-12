from pathlib import Path
import csv
import subprocess

from rdkit import Chem


def xyz_to_smiles(xyz_path: Path) -> tuple[str | None, str | None]:
    """
    Convert XYZ -> SMILES using OpenBabel + RDKit.
    Returns (smiles, error_message).
    """
    cmd = ["obabel", str(xyz_path), "-osmi"]
    result = subprocess.run(cmd, capture_output=True, text=True)

    if result.returncode != 0:
        return None, f"OpenBabel error: {result.stderr.strip()}"

    out = result.stdout.strip()
    if not out:
        return None, "OpenBabel produced no output"

    first_line = out.splitlines()[0]
    parts = first_line.split()
    if not parts:
        return None, f"Could not parse OpenBabel line: {first_line}"

    obabel_smiles = parts[0]

    # RDKit cleanup
    mol = Chem.MolFromSmiles(obabel_smiles)
    if mol is None:
        return None, f"RDKit could not parse SMILES: {obabel_smiles}"

    canon = Chem.MolToSmiles(mol)
    return canon, None


def main():
    repo_root = Path(__file__).resolve().parents[1]
    xyz_root = repo_root / "data" / "raw" / "xyz"
    out_csv = repo_root / "data" / "processed" / "xyz_smiles_map.csv"

    print("Scanning XYZ files in:", xyz_root)
    if not xyz_root.exists():
        print("Folder does not exist! Check data/raw/xyz path.")
        return

    xyz_paths = sorted(xyz_root.rglob("*.xyz"))
    print(f"Found {len(xyz_paths)} XYZ files")

    rows = []
    for i, xyz_path in enumerate(xyz_paths, start=1):
        print(f"[{i}/{len(xyz_paths)}] {xyz_path}")
        smiles, error = xyz_to_smiles(xyz_path)
        rows.append(
            {
                "xyz_path": str(xyz_path.relative_to(repo_root)),
                "smiles": smiles,
                "error": error,
            }
        )

    out_csv.parent.mkdir(parents=True, exist_ok=True)

    fieldnames = ["xyz_path", "smiles", "error"]
    with open(out_csv, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    print("Done. Wrote:", out_csv)


if __name__ == "__main__":
    main()
