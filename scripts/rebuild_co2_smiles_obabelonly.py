from pathlib import Path
import csv
import subprocess


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
    repo_root = Path(__file__).resolve().parents[1]
    xyz_root = repo_root / "data" / "raw" / "xyz"
    out_csv = repo_root / "data" / "processed" / "co2_xyz_smiles_obabelonly.csv"

    print("Scanning XYZ files in:", xyz_root)
    xyz_paths = sorted(xyz_root.rglob("*.xyz"))
    print(f"Found {len(xyz_paths)} total XYZ files")

    rows = []
    idx = 0
    for xyz_path in xyz_paths:
        filename = xyz_path.name

        # Only CO2 files; adjust condition if your naming is different
        if "CO2" not in filename.upper():
            continue

        idx += 1
        print(f"[{idx}] {xyz_path}")

        smiles, error = obabel_xyz_to_smiles_raw(xyz_path)

        rows.append(
            {
                "xyz_path": str(xyz_path.relative_to(repo_root)),
                "filename": filename,
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
