"""
Registry builder for FLPCO2DB.

Generates per-entry YAML files and a central registry from raw data.
"""

import hashlib
import yaml
from pathlib import Path
from datetime import datetime
from typing import Dict, List, Optional
import pandas as pd
import numpy as np
from rich.console import Console
from rich.progress import track

from . import RAW_DATA_DIR, PROCESSED_DATA_DIR, ENTRIES_DIR, REGISTRY_PATH
from .staging import build_xyz_inventory, compute_file_sha256
from .smiles_utils import infer_smiles_from_xyz, read_xyz_coordinates

console = Console()


def convert_numpy_types(obj):
    """
    Recursively convert numpy types to native Python types.

    This is necessary for YAML serialization because yaml.safe_load()
    cannot deserialize numpy-specific types.

    Args:
        obj: Object that may contain numpy types

    Returns:
        Object with all numpy types converted to Python types
    """
    if isinstance(obj, dict):
        return {key: convert_numpy_types(value) for key, value in obj.items()}
    elif isinstance(obj, list):
        return [convert_numpy_types(item) for item in obj]
    elif isinstance(obj, np.integer):
        return int(obj)
    elif isinstance(obj, np.floating):
        return float(obj)
    elif isinstance(obj, np.ndarray):
        return obj.tolist()
    elif isinstance(obj, np.bool_):
        return bool(obj)
    else:
        return obj


def load_co2_energies(graphs_csv_dir: Path) -> pd.DataFrame:
    """
    Load CO2 energy data from graphs CSV.
    
    Args:
        graphs_csv_dir: Path to graphs_csv directory
        
    Returns:
        DataFrame with CO2 energy data
    """
    co2_csv = graphs_csv_dir / "CO2.csv"
    if not co2_csv.exists():
        raise FileNotFoundError(f"CO2.csv not found in {graphs_csv_dir}")
    
    df = pd.read_csv(co2_csv)
    
    # Rename OLD ID to flp_id for consistency
    if "OLD ID" in df.columns:
        df = df.rename(columns={"OLD ID": "flp_id"})
    
    return df


def compute_la_lb_distance(xyz_path: Path) -> Optional[float]:
    """
    Compute Lewis acid - Lewis base distance from XYZ file.
    
    Note: This is a placeholder. Actual implementation would require
    identifying LA and LB atoms, which needs chemical knowledge.
    For now, we'll return None and students can implement this later.
    
    Args:
        xyz_path: Path to XYZ file
        
    Returns:
        Distance in Angstroms or None
    """
    # TODO: Implement LA-LB identification and distance calculation
    # This requires domain knowledge about which atoms are LA/LB
    # See mml_studio_07/utils.py for potential patterns
    return None


def build_entry(
    flp_id: int,
    xyz_inventory: pd.DataFrame,
    co2_energies: pd.DataFrame,
    graphs_csv_dir: Path,
    html_pages_dir: Path,
    validate_smiles: bool = True
) -> Dict[str, any]:
    """
    Build a single registry entry for an FLP.
    
    Args:
        flp_id: FLP identifier
        xyz_inventory: DataFrame with XYZ file inventory
        co2_energies: DataFrame with CO2 energy data
        graphs_csv_dir: Path to graphs CSV directory
        html_pages_dir: Path to HTML pages directory
        validate_smiles: Whether to validate SMILES with round-trip check
        
    Returns:
        Dictionary with entry data
    """
    entry = {
        "flp_id": flp_id,
        "flp_code": None,
        "provenance": {
            "xyz_paths": {},
            "graphs_row_source": str(graphs_csv_dir / "CO2.csv"),
            "html_page": None,
        },
        "structure": {
            "la_lb_distance_A": None,
            "smiles": {},
        },
        "energies_recovered": {
            "units": "kcal/mol",
            "gas": {},
            "solution": {},
        },
        "compute_plan": {
            "status": "not_started",
            "notes": "Calculations to be performed by students using mml_studio_07 patterns or custom workflows",
        },
        "qc_flags": {
            "has_xyz_flp": False,
            "has_xyz_co2": False,
            "has_recovered_energy": False,
            "join_ok": False,
            "smiles_validation_passed": None,
        },
        "notes": None,
    }
    
    # Get XYZ files for this FLP
    flp_xyz = xyz_inventory[xyz_inventory["flp_id"] == flp_id]
    
    # Get bare FLP path
    flp_bare = flp_xyz[flp_xyz["file_type"] == "flp_bare"]
    if not flp_bare.empty:
        entry["provenance"]["xyz_paths"]["flp"] = flp_bare.iloc[0]["file_path"]
        entry["qc_flags"]["has_xyz_flp"] = True
        
        # Infer SMILES for bare FLP
        flp_xyz_path = Path(entry["provenance"]["xyz_paths"]["flp"])
        if flp_xyz_path.exists():
            smiles_result = infer_smiles_from_xyz(flp_xyz_path, validate=validate_smiles)
            entry["structure"]["smiles"]["flp_bare"] = {
                "rdkit": smiles_result["rdkit_smiles"],
                "openbabel": smiles_result["openbabel_smiles"],
                "inchi": smiles_result["inchi"],
                "inchikey": smiles_result["inchikey"],
                "validation": smiles_result["validation"],
            }
            
            if smiles_result["validation"]["attempted"]:
                entry["qc_flags"]["smiles_validation_passed"] = smiles_result["validation"]["passed"]
            
            # Compute LA-LB distance
            entry["structure"]["la_lb_distance_A"] = compute_la_lb_distance(flp_xyz_path)
    
    # Get CO2 adduct path
    co2_adduct = flp_xyz[flp_xyz["file_type"] == "co2_adduct"]
    if not co2_adduct.empty:
        entry["provenance"]["xyz_paths"]["co2"] = co2_adduct.iloc[0]["file_path"]
        entry["qc_flags"]["has_xyz_co2"] = True
        
        # Infer SMILES for CO2 adduct
        co2_xyz_path = Path(entry["provenance"]["xyz_paths"]["co2"])
        if co2_xyz_path.exists():
            smiles_result = infer_smiles_from_xyz(co2_xyz_path, validate=validate_smiles)
            
            entry["structure"]["smiles"]["co2_adducts"] = [{
                "tag": "bidentate_default",
                "rdkit": smiles_result["rdkit_smiles"],
                "openbabel": smiles_result["openbabel_smiles"],
                "inchi": smiles_result["inchi"],
                "inchikey": smiles_result["inchikey"],
                "xyz_source": entry["provenance"]["xyz_paths"]["co2"],
                "validation": smiles_result["validation"],
            }]
    
    # Get energy data
    energy_row = co2_energies[co2_energies["flp_id"] == flp_id]
    if not energy_row.empty:
        entry["qc_flags"]["has_recovered_energy"] = True
        row = energy_row.iloc[0]
        
        # Extract FLP code if available
        if "FLP" in row and pd.notna(row["FLP"]):
            entry["flp_code"] = str(row["FLP"])
        
        # Extract energies (column names may vary, adjust as needed)
        # Gas phase
        for key in ["E", "H", "G"]:
            col_name = f"dG_{key}_gas"  # Adjust based on actual column names
            if col_name in row and pd.notna(row[col_name]):
                entry["energies_recovered"]["gas"][key] = float(row[col_name])
        
        # Solution phase
        for key in ["E", "H", "G"]:
            col_name = f"dG_{key}_sol"  # Adjust based on actual column names
            if col_name in row and pd.notna(row[col_name]):
                entry["energies_recovered"]["solution"][key] = float(row[col_name])
        
        # If energies not found with those column names, try alternatives
        if not entry["energies_recovered"]["gas"] and not entry["energies_recovered"]["solution"]:
            # Try simpler column names
            for key in row.index:
                if "gas" in key.lower() and any(e in key for e in ["E", "H", "G"]):
                    value = row[key]
                    if pd.notna(value):
                        entry["energies_recovered"]["gas"][key.split("_")[-1]] = float(value)
                elif "sol" in key.lower() and any(e in key for e in ["E", "H", "G"]):
                    value = row[key]
                    if pd.notna(value):
                        entry["energies_recovered"]["solution"][key.split("_")[-1]] = float(value)
    
    # Check HTML page
    html_file = html_pages_dir / f"{flp_id}.html"
    if html_file.exists():
        entry["provenance"]["html_page"] = str(html_file.relative_to(RAW_DATA_DIR.parent.parent))
    
    # Set join_ok flag
    entry["qc_flags"]["join_ok"] = (
        entry["qc_flags"]["has_xyz_flp"] and
        (entry["qc_flags"]["has_xyz_co2"] or entry["qc_flags"]["has_recovered_energy"])
    )
    
    return entry


def build_registry(
    validate_smiles: bool = True,
    max_entries: Optional[int] = None
) -> Dict[str, any]:
    """
    Build the complete CO2 registry.
    
    Args:
        validate_smiles: Whether to validate SMILES with round-trip checks
        max_entries: Maximum number of entries to process (for testing)
        
    Returns:
        Dictionary with registry data
    """
    console.print("[bold]Building FLPCO2DB Registry[/bold]\n")
    
    # Load raw data
    console.print("Loading raw data...")
    graphs_csv_dir = RAW_DATA_DIR / "graphs_csv"
    html_pages_dir = RAW_DATA_DIR / "html_pages"
    
    xyz_inventory = build_xyz_inventory(RAW_DATA_DIR)
    co2_energies = load_co2_energies(graphs_csv_dir)
    
    # Get list of FLP IDs
    flp_ids = sorted(xyz_inventory["flp_id"].unique())
    if max_entries:
        flp_ids = flp_ids[:max_entries]
        console.print(f"[yellow]Processing only first {max_entries} entries for testing[/yellow]")
    
    console.print(f"Processing {len(flp_ids)} FLP entries...")
    
    # Build entries
    entries = []
    failed_entries = []
    
    for flp_id in track(flp_ids, description="Building entries"):
        try:
            entry = build_entry(
                flp_id,
                xyz_inventory,
                co2_energies,
                graphs_csv_dir,
                html_pages_dir,
                validate_smiles=validate_smiles
            )
            entries.append(entry)
            
            # Save per-entry YAML
            entry_path = ENTRIES_DIR / f"{flp_id}.yaml"
            entry_path.parent.mkdir(parents=True, exist_ok=True)
            with open(entry_path, 'w') as f:
                # Convert numpy types to native Python types before saving
                entry_clean = convert_numpy_types(entry)
                yaml.dump(entry_clean, f, default_flow_style=False, sort_keys=False)
        
        except Exception as e:
            console.print(f"[red]Error processing FLP {flp_id}:[/red] {e}")
            failed_entries.append((flp_id, str(e)))
    
    # Build summary statistics
    stats = {
        "flps_total": len(entries),
        "with_xyz_co2": sum(1 for e in entries if e["qc_flags"]["has_xyz_co2"]),
        "with_energy_co2": sum(1 for e in entries if e["qc_flags"]["has_recovered_energy"]),
        "overlap": sum(1 for e in entries if e["qc_flags"]["has_xyz_co2"] and e["qc_flags"]["has_recovered_energy"]),
        "smiles_validated": sum(1 for e in entries if e["qc_flags"].get("smiles_validation_passed") is True),
        "smiles_failed": sum(1 for e in entries if e["qc_flags"].get("smiles_validation_passed") is False),
    }
    
    # Build central registry
    registry = {
        "generated_at": datetime.now().isoformat(),
        "dataset_version": "0.1.0",
        "source_zip": "FLPDB-d549732.zip",
        "source_zip_sha256": None,  # Compute if needed
        "method_profile": "raw_data_only",
        "notes": "Initial registry built from raw FLPDB data. No calculations performed yet.",
        "counts": stats,
        "entries": [
            {
                "flp_id": e["flp_id"],
                "flp_code": e["flp_code"],
                "entry_file": f"entries/{e['flp_id']}.yaml",
                "qc_flags": e["qc_flags"],
            }
            for e in entries
        ],
        "failed_entries": failed_entries,
    }
    
    # Save central registry
    REGISTRY_PATH.parent.mkdir(parents=True, exist_ok=True)
    with open(REGISTRY_PATH, 'w') as f:
        # Convert numpy types to native Python types before saving
        registry_clean = convert_numpy_types(registry)
        yaml.dump(registry_clean, f, default_flow_style=False, sort_keys=False)
    
    console.print(f"\n[bold green]✓ Registry built successfully![/bold green]")
    console.print(f"  Entries: {stats['flps_total']}")
    console.print(f"  With CO₂ XYZ: {stats['with_xyz_co2']}")
    console.print(f"  With CO₂ energy: {stats['with_energy_co2']}")
    console.print(f"  Overlap: {stats['overlap']}")
    if validate_smiles:
        console.print(f"  SMILES validated: {stats['smiles_validated']}")
        console.print(f"  SMILES failed: {stats['smiles_failed']}")
    
    if failed_entries:
        console.print(f"\n[yellow]Warning: {len(failed_entries)} entries failed[/yellow]")
    
    console.print(f"\nRegistry saved to: {REGISTRY_PATH}")
    console.print(f"Entry files saved to: {ENTRIES_DIR}/")
    
    return registry


if __name__ == "__main__":
    # Build the registry
    import argparse
    
    parser = argparse.ArgumentParser(description="Build FLPCO2DB registry")
    parser.add_argument(
        "--no-validate",
        action="store_true",
        help="Skip SMILES validation (faster)"
    )
    parser.add_argument(
        "--max-entries",
        type=int,
        help="Maximum number of entries to process (for testing)"
    )
    
    args = parser.parse_args()
    
    registry = build_registry(
        validate_smiles=not args.no_validate,
        max_entries=args.max_entries
    )

