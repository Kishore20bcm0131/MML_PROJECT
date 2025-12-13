from pathlib import Path
import csv
import subprocess


def is_flp_xyz(xyz_path: Path) -> bool:
    """
    Return True if this XYZ is the bare FLP, not a CO2 adduct.

    Rule:
    - Skip files whose name contains 'CO2'
    - Keep files where the filename stem == folder name, e.g.
      data/raw/xyz/1/1.xyz   -> folder '1', stem '1'  -> FLP
      data/raw/xyz/1/1CO2.xyz -> folder '1', stem '1CO2' -> NOT FLP
    """
    name = xyz_path.name  # e.g. '1.xyz', '1CO2.xyz'
    if "CO2" in name.upper():
        return False

    folder = xyz_path.parent.name      # e.g. '1'
    stem = xyz_path.stem               # e.g. '1' or '1CO2'
    return stem == folder              # FLP if '1' == '1'


def obabel_xyz_to_smiles_raw(xyz_path: Path) -> tuple[str | None, str | None]:
    """
    XYZ -> SMILES using OpenBabel only.
    Returns (smiles, error_message).
    We DO NOT run RDKit here; we keep the raw OpenBabel SMILES.
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
    return obabel_smiles, None


def main():
    # repo root = parent of scripts/
    repo_root = Path(__file__).resolve().parents[1]
    xyz_root = repo_root / "data" / "raw" / "xyz"
    out_csv = repo_root / "data" / "processed" / "flp_xyz_smiles_obabelonly.csv"

    print("Scanning XYZ files in:", xyz_root)
    xyz_paths = sorted(xyz_root.rglob("*.xyz"))
    print(f"Found {len(xyz_paths)} total XYZ files")

    rows = []
    idx = 0
    for xyz_path in xyz_paths:
        if not is_flp_xyz(xyz_path):
            continue  # skip CO2 and other non-FLPs

        idx += 1
        print(f"[{idx}] {xyz_path}")

        smiles, error = obabel_xyz_to_smiles_raw(xyz_path)

        rows.append(
            {
                "xyz_path": str(xyz_path.relative_to(repo_root)),
                "filename": xyz_path.name,
                "smiles": smiles,
                "error": error,
            }
        )

    out_csv.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = ["xyz_path", "filename", "smiles", "error"]
    with open(out_csv, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    print("Done. Wrote:", out_csv)


if __name__ == "__main__":
    main()
