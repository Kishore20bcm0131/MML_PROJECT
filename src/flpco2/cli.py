"""
Command-line interface for FLPCO2DB.

Provides commands for data management, registry building, and inspection.
"""

import typer
import yaml
from pathlib import Path
from rich.console import Console
from rich.table import Table
from rich import print_json
import json

from . import REGISTRY_PATH, ENTRIES_DIR, RAW_DATA_DIR, PROCESSED_DATA_DIR
from .staging import verify_raw_data_structure, print_verification_report, build_xyz_inventory
from .registry_builder import build_registry

app = typer.Typer(
    name="flpco2",
    help="FLPCO2DB - Frustrated Lewis Pairs CO₂ Database Management Tool",
    add_completion=False,
)

console = Console()


@app.command()
def stage(verify: bool = typer.Option(True, help="Verify staged data")):
    """
    Stage and verify raw data from FLPDB snapshot.
    
    This command checks that all expected files are present:
    - XYZ geometries (133 FLPs, 927 files total, 132 with CO₂)
    - CSV files with energies and metadata
    - HTML pages
    """
    console.print("[bold]FLPCO2DB Data Staging[/bold]\n")
    
    if not RAW_DATA_DIR.exists():
        console.print(f"[red]Error:[/red] Raw data directory not found: {RAW_DATA_DIR}")
        console.print("\nExpected structure:")
        console.print("  data/raw/xyz/")
        console.print("  data/raw/graphs_csv/")
        console.print("  data/raw/html_pages/")
        raise typer.Exit(1)
    
    if verify:
        console.print("Verifying data structure...")
        stats = verify_raw_data_structure(RAW_DATA_DIR)
        print_verification_report(stats)
        
        if stats["errors"]:
            console.print("\n[red]Data verification failed. Please check the errors above.[/red]")
            raise typer.Exit(1)
        else:
            console.print("\n[bold green]✓ Data staging verified successfully![/bold green]")


@app.command()
def build_reg(
    validate: bool = typer.Option(True, help="Validate SMILES with round-trip checks"),
    max_entries: int = typer.Option(None, help="Maximum entries to process (for testing)"),
    output_format: str = typer.Option("yaml", help="Output format: yaml or json")
):
    """
    Build the CO₂ registry from raw data.
    
    This generates:
    - Per-entry YAML files in data/processed/entries/
    - Central registry at data/processed/co2_registry.yaml
    
    Each entry includes:
    - Provenance (XYZ paths, source files)
    - SMILES, InChI, InChIKey for structures
    - Recovered energies from original data
    - QC flags
    """
    try:
        registry = build_registry(
            validate_smiles=validate,
            max_entries=max_entries
        )
        
        # Optionally save as JSON as well
        if output_format == "json":
            json_path = REGISTRY_PATH.with_suffix(".json")
            with open(json_path, 'w') as f:
                json.dump(registry, f, indent=2)
            console.print(f"Registry also saved as JSON: {json_path}")
        
        console.print("\n[green]✓ Registry built successfully![/green]")
        
    except Exception as e:
        console.print(f"\n[red]Error building registry:[/red] {e}")
        raise typer.Exit(1)


@app.command()
def validate():
    """
    Validate the registry against expected schema and counts.
    
    Checks:
    - Registry file exists and is valid YAML
    - Expected overlap counts (132 entries with XYZ + energy)
    - Entry files exist for all listed entries
    - QC flags are consistent
    """
    console.print("[bold]Validating FLPCO2DB Registry[/bold]\n")
    
    if not REGISTRY_PATH.exists():
        console.print(f"[red]Error:[/red] Registry not found: {REGISTRY_PATH}")
        console.print("\nRun 'flpco2 build-reg' first to create the registry.")
        raise typer.Exit(1)
    
    # Load registry
    try:
        with open(REGISTRY_PATH, 'r') as f:
            registry = yaml.safe_load(f)
    except Exception as e:
        console.print(f"[red]Error loading registry:[/red] {e}")
        raise typer.Exit(1)
    
    # Validate counts
    counts = registry.get("counts", {})
    
    table = Table(title="Registry Validation", show_header=True, header_style="bold magenta")
    table.add_column("Check", style="cyan")
    table.add_column("Value", justify="right")
    table.add_column("Expected", justify="right")
    table.add_column("Status", justify="center")
    
    # Check overlap
    overlap = counts.get("overlap", 0)
    status = "✓" if overlap >= 132 else "✗"
    color = "green" if overlap >= 132 else "red"
    table.add_row(
        "Entries with XYZ + Energy",
        str(overlap),
        "132",
        f"[{color}]{status}[/{color}]"
    )
    
    # Check total entries
    total = counts.get("flps_total", 0)
    status = "✓" if total >= 133 else "✗"
    color = "green" if total >= 133 else "yellow"
    table.add_row(
        "Total FLP Entries",
        str(total),
        "133",
        f"[{color}]{status}[/{color}]"
    )
    
    # Check entry files
    entries = registry.get("entries", [])
    missing_files = []
    for entry in entries:
        entry_file = PROCESSED_DATA_DIR / entry["entry_file"]
        if not entry_file.exists():
            missing_files.append(entry_file)
    
    status = "✓" if not missing_files else "✗"
    color = "green" if not missing_files else "red"
    table.add_row(
        "Entry Files Present",
        f"{len(entries) - len(missing_files)}/{len(entries)}",
        f"{len(entries)}",
        f"[{color}]{status}[/{color}]"
    )
    
    console.print(table)
    
    if missing_files:
        console.print(f"\n[yellow]Warning:[/yellow] {len(missing_files)} entry files missing")
        console.print("[red]Validation failed.[/red]")
        raise typer.Exit(1)
    else:
        console.print("\n[bold green]✓ Registry validation passed![/bold green]")


@app.command()
def inspect(
    flp_id: int = typer.Argument(..., help="FLP ID to inspect"),
    format: str = typer.Option("yaml", help="Output format: yaml or json")
):
    """
    Inspect a specific FLP entry.
    
    Displays the complete entry information including:
    - Provenance and file paths
    - SMILES and molecular identifiers
    - Recovered energies
    - QC flags
    """
    entry_file = ENTRIES_DIR / f"{flp_id}.yaml"
    
    if not entry_file.exists():
        console.print(f"[red]Error:[/red] Entry not found: {entry_file}")
        console.print(f"\nFLP ID {flp_id} may not exist in the registry.")
        raise typer.Exit(1)
    
    # Load entry
    with open(entry_file, 'r') as f:
        entry = yaml.safe_load(f)
    
    console.print(f"\n[bold]FLP Entry: {flp_id}[/bold]")
    console.print(f"Code: {entry.get('flp_code', 'N/A')}\n")
    
    if format == "json":
        print_json(data=entry)
    else:
        # Pretty print YAML
        console.print(yaml.dump(entry, default_flow_style=False, sort_keys=False))


@app.command()
def stats():
    """
    Display summary statistics about the dataset.
    
    Shows:
    - Total entries
    - Coverage (XYZ files, energies, overlap)
    - SMILES validation results
    - Missing data summary
    """
    if not REGISTRY_PATH.exists():
        console.print(f"[red]Error:[/red] Registry not found: {REGISTRY_PATH}")
        console.print("\nRun 'flpco2 build-reg' first to create the registry.")
        raise typer.Exit(1)
    
    # Load registry
    with open(REGISTRY_PATH, 'r') as f:
        registry = yaml.safe_load(f)
    
    counts = registry.get("counts", {})
    
    console.print("\n[bold]FLPCO2DB Statistics[/bold]\n")
    
    table = Table(show_header=True, header_style="bold magenta")
    table.add_column("Metric", style="cyan")
    table.add_column("Count", justify="right", style="green")
    
    table.add_row("Total FLP Entries", str(counts.get("flps_total", 0)))
    table.add_row("With CO₂ XYZ Files", str(counts.get("with_xyz_co2", 0)))
    table.add_row("With CO₂ Energies", str(counts.get("with_energy_co2", 0)))
    table.add_row("Overlap (XYZ + Energy)", str(counts.get("overlap", 0)))
    
    if "smiles_validated" in counts:
        table.add_row("SMILES Validated", str(counts.get("smiles_validated", 0)))
        table.add_row("SMILES Failed", str(counts.get("smiles_failed", 0)))
    
    console.print(table)
    
    # Show failed entries if any
    failed = registry.get("failed_entries", [])
    if failed:
        console.print(f"\n[yellow]Warning:[/yellow] {len(failed)} entries failed during registry build")


@app.command()
def export(
    output_path: Path = typer.Option(..., "--output", "-o", help="Output file path"),
    format: str = typer.Option("csv", help="Export format: csv, json")
):
    """
    Export registry data to various formats for analysis.
    
    Exports a flattened view of the registry suitable for:
    - CSV: Data analysis, spreadsheets
    - JSON: Programmatic access, web APIs
    """
    if not REGISTRY_PATH.exists():
        console.print(f"[red]Error:[/red] Registry not found: {REGISTRY_PATH}")
        raise typer.Exit(1)
    
    # Load all entries
    with open(REGISTRY_PATH, 'r') as f:
        registry = yaml.safe_load(f)
    
    entries_data = []
    for entry_ref in registry.get("entries", []):
        entry_file = PROCESSED_DATA_DIR / entry_ref["entry_file"]
        if entry_file.exists():
            with open(entry_file, 'r') as f:
                entry = yaml.safe_load(f)
                entries_data.append(entry)
    
    # Flatten for export
    flattened = []
    for entry in entries_data:
        flat = {
            "flp_id": entry["flp_id"],
            "flp_code": entry.get("flp_code"),
            "has_xyz_flp": entry["qc_flags"]["has_xyz_flp"],
            "has_xyz_co2": entry["qc_flags"]["has_xyz_co2"],
            "has_energy": entry["qc_flags"]["has_recovered_energy"],
            "smiles_flp": entry["structure"]["smiles"].get("flp_bare", {}).get("rdkit"),
            "inchikey_flp": entry["structure"]["smiles"].get("flp_bare", {}).get("inchikey"),
        }
        flattened.append(flat)
    
    # Export
    import pandas as pd
    df = pd.DataFrame(flattened)
    
    if format == "csv":
        df.to_csv(output_path, index=False)
    elif format == "json":
        df.to_json(output_path, orient="records", indent=2)
    else:
        console.print(f"[red]Unknown format:[/red] {format}")
        raise typer.Exit(1)
    
    console.print(f"[green]✓ Exported {len(df)} entries to {output_path}[/green]")


@app.command()
def version():
    """Show version information."""
    from . import __version__
    console.print(f"FLPCO2DB version {__version__}")


def main():
    """Main entry point for the CLI."""
    app()


if __name__ == "__main__":
    main()

