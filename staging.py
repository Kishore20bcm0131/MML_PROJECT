"""
Data staging utilities for FLPCO2DB

Functions to organize and validate raw data from the FLPDB snapshot.
"""

import os
import hashlib
from pathlib import Path
from typing import Dict, List, Tuple
import pandas as pd
from rich.console import Console
from rich.table import Table

console = Console()


def compute_file_sha256(filepath: Path) -> str:
    """Compute SHA256 hash of a file."""
    sha256_hash = hashlib.sha256()
    with open(filepath, "rb") as f:
        for byte_block in iter(lambda: f.read(4096), b""):
            sha256_hash.update(byte_block)
    return sha256_hash.hexdigest()


def verify_raw_data_structure(raw_data_dir: Path) -> Dict[str, any]:
    """
    Verify the raw data directory structure and count files.
    
    Args:
        raw_data_dir: Path to data/raw directory
        
    Returns:
        Dictionary with verification statistics
    """
    stats = {
        "flp_directories": 0,
        "total_xyz_files": 0,
        "co2_xyz_files": 0,
        "csv_files": [],
        "html_files": 0,
        "errors": [],
    }
    
    # Check XYZ directory
    xyz_dir = raw_data_dir / "xyz"
    if xyz_dir.exists():
        flp_dirs = [d for d in xyz_dir.iterdir() if d.is_dir()]
        stats["flp_directories"] = len(flp_dirs)
        
        for flp_dir in flp_dirs:
            xyz_files = list(flp_dir.glob("*.xyz"))
            stats["total_xyz_files"] += len(xyz_files)
            
            # Check for CO2 xyz
            co2_files = [f for f in xyz_files if "CO2" in f.stem]
            stats["co2_xyz_files"] += len(co2_files)
    else:
        stats["errors"].append(f"xyz directory not found: {xyz_dir}")
    
    # Check CSV files
    graphs_csv_dir = raw_data_dir / "graphs_csv"
    if graphs_csv_dir.exists():
        stats["csv_files"] = [f.name for f in graphs_csv_dir.glob("*.csv")]
    else:
        stats["errors"].append(f"graphs_csv directory not found: {graphs_csv_dir}")
    
    # Check HTML files
    html_dir = raw_data_dir / "html_pages"
    if html_dir.exists():
        stats["html_files"] = len(list(html_dir.glob("*.html")))
    else:
        stats["errors"].append(f"html_pages directory not found: {html_dir}")
    
    return stats


def build_xyz_inventory(raw_data_dir: Path) -> pd.DataFrame:
    """
    Build an inventory of all XYZ files organized by FLP ID.
    
    Args:
        raw_data_dir: Path to data/raw directory
        
    Returns:
        DataFrame with columns: flp_id, file_type, file_path
    """
    xyz_dir = raw_data_dir / "xyz"
    inventory = []
    
    for flp_dir in sorted(xyz_dir.iterdir()):
        if not flp_dir.is_dir():
            continue
            
        flp_id = flp_dir.name
        
        for xyz_file in flp_dir.glob("*.xyz"):
            # Determine file type based on name
            filename = xyz_file.stem
            
            if filename == flp_id:
                file_type = "flp_bare"
            elif "CO2" in filename:
                file_type = "co2_adduct"
            elif "H2O" in filename:
                file_type = "h2o_adduct"
            elif "H2" in filename and "H2O" not in filename and "H2CO" not in filename:
                file_type = "h2_adduct"
            elif "CH3OH" in filename:
                file_type = "ch3oh_adduct"
            elif "H2CO" in filename:
                file_type = "h2co_adduct"
            elif "HCOOH" in filename:
                file_type = "hcooh_adduct"
            else:
                file_type = "other"
            
            inventory.append({
                "flp_id": int(flp_id),
                "file_type": file_type,
                "file_path": str(xyz_file.relative_to(raw_data_dir.parent.parent)),
            })
    
    return pd.DataFrame(inventory)


def print_verification_report(stats: Dict[str, any]):
    """Print a formatted verification report."""
    console.print("\n[bold]FLPCO2DB Data Verification Report[/bold]\n")
    
    table = Table(show_header=True, header_style="bold magenta")
    table.add_column("Metric", style="cyan")
    table.add_column("Count", style="green", justify="right")
    table.add_column("Expected", style="yellow", justify="right")
    table.add_column("Status", justify="center")
    
    # FLP directories
    status = "✓" if stats["flp_directories"] >= 133 else "✗"
    table.add_row(
        "FLP Directories",
        str(stats["flp_directories"]),
        "133",
        f"[green]{status}[/green]" if status == "✓" else f"[red]{status}[/red]"
    )
    
    # Total XYZ files
    status = "✓" if stats["total_xyz_files"] >= 927 else "✗"
    table.add_row(
        "Total XYZ Files",
        str(stats["total_xyz_files"]),
        "927",
        f"[green]{status}[/green]" if status == "✓" else f"[red]{status}[/red]"
    )
    
    # CO2 XYZ files
    status = "✓" if stats["co2_xyz_files"] >= 132 else "✗"
    table.add_row(
        "CO₂ XYZ Files",
        str(stats["co2_xyz_files"]),
        "132",
        f"[green]{status}[/green]" if status == "✓" else f"[red]{status}[/red]"
    )
    
    # HTML files
    status = "✓" if stats["html_files"] >= 100 else "?"
    table.add_row(
        "HTML Pages",
        str(stats["html_files"]),
        "~133",
        f"[yellow]{status}[/yellow]"
    )
    
    console.print(table)
    
    # CSV files
    if stats["csv_files"]:
        console.print(f"\n[bold]CSV Files Found:[/bold] {', '.join(stats['csv_files'])}")
    
    # Errors
    if stats["errors"]:
        console.print("\n[bold red]Errors:[/bold red]")
        for error in stats["errors"]:
            console.print(f"  • {error}")
    else:
        console.print("\n[bold green]✓ All expected directories found![/bold green]")


if __name__ == "__main__":
    # Test the staging functions
    from flpco2 import RAW_DATA_DIR
    
    console.print("[bold]Running data verification...[/bold]")
    stats = verify_raw_data_structure(RAW_DATA_DIR)
    print_verification_report(stats)
    
    console.print("\n[bold]Building XYZ inventory...[/bold]")
    inventory = build_xyz_inventory(RAW_DATA_DIR)
    console.print(f"Inventory contains {len(inventory)} XYZ files")
    console.print(f"File types: {inventory['file_type'].value_counts().to_dict()}")

