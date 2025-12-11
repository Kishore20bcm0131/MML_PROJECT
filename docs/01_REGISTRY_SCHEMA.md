# Registry Schema Documentation

This document describes the structure and format of the FLPCO2DB registry.

## Overview

The registry consists of two parts:
1. **Central registry** (`data/processed/co2_registry.yaml`): Summary and index
2. **Per-entry files** (`data/processed/entries/<id>.yaml`): Detailed data for each FLP

## Central Registry Schema

Location: `data/processed/co2_registry.yaml`

```yaml
generated_at: "2025-11-06T12:00:00"  # ISO 8601 timestamp
dataset_version: "0.1.0"              # Semantic version
source_zip: "FLPDB-d549732.zip"       # Source file name
source_zip_sha256: null               # Optional: SHA256 hash
method_profile: "raw_data_only"       # Computation status

notes: "Initial registry built from raw FLPDB data..."

counts:
  flps_total: 133                     # Total FLP entries
  with_xyz_co2: 132                   # Entries with CO₂ XYZ
  with_energy_co2: 137                # Entries with CO₂ energies
  overlap: 132                        # Entries with both
  smiles_validated: 120               # SMILES round-trip passed
  smiles_failed: 12                   # SMILES round-trip failed

entries:
  - flp_id: 108
    flp_code: "0BN08011"
    entry_file: "entries/108.yaml"
    qc_flags:
      has_xyz_flp: true
      has_xyz_co2: true
      has_recovered_energy: true
      join_ok: true
      smiles_validation_passed: true
  # ... more entries ...

failed_entries: []  # List of (flp_id, error_message) tuples
```

## Per-Entry Schema

Location: `data/processed/entries/<flp_id>.yaml`

### Complete Example

```yaml
flp_id: 108
flp_code: "0BN08011"

provenance:
  xyz_paths:
    flp: "data/raw/xyz/108/108.xyz"
    co2: "data/raw/xyz/108/108CO2.xyz"
  graphs_row_source: "data/raw/graphs_csv/CO2.csv"
  html_page: "data/raw/html_pages/108.html"

structure:
  la_lb_distance_A: null  # To be computed by students
  smiles:
    flp_bare:
      rdkit: "CC1=CC=C(C=C1)P(C2=CC=CC=C2)C3=CC=CC=C3..."
      openbabel: "CC1=CC=C(C=C1)P(C2=CC=CC=C2)C3=CC=CC=C3..."
      inchi: "InChI=1S/C19H17P/c1-16-12-14-17(15-13-16)..."
      inchikey: "XXXXXXXXXXXXXXXXXXX-UHFFFAOYSA-N"
      validation:
        attempted: true
        passed: true
        rmsd: 0.234
        error: null
    co2_adducts:
      - tag: "bidentate_default"
        rdkit: "CC1=CC=C(C=C1)P(C2=CC=CC=C2)C3=CC=CC=C3..."
        openbabel: "..."
        inchi: "..."
        inchikey: "..."
        xyz_source: "data/raw/xyz/108/108CO2.xyz"
        validation:
          attempted: true
          passed: false
          rmsd: 1.234
          error: "RMSD 1.23 Å exceeds threshold 1.0 Å"

energies_recovered:
  units: "kcal/mol"
  gas:
    E: -12.34
    H: -11.23
    G: -10.12
  solution:
    E: -8.91
    H: -7.80
    G: 6.70

compute_plan:
  status: "not_started"
  notes: "Calculations to be performed by students..."

qc_flags:
  has_xyz_flp: true
  has_xyz_co2: true
  has_recovered_energy: true
  join_ok: true
  smiles_validation_passed: true

notes: null
```

## Field Descriptions

### Top-Level Fields

| Field | Type | Description |
|-------|------|-------------|
| `flp_id` | integer | Unique FLP identifier from FLPDB |
| `flp_code` | string | Alphanumeric code from original dataset (e.g., "0BN08011") |
| `provenance` | object | Source file paths and metadata |
| `structure` | object | Molecular structure information |
| `energies_recovered` | object | Original energy data from FLPDB |
| `compute_plan` | object | Status of calculations to be performed |
| `qc_flags` | object | Quality control and data availability flags |
| `notes` | string | Additional notes or warnings |

### Provenance

| Field | Type | Description |
|-------|------|-------------|
| `xyz_paths.flp` | string | Path to bare FLP XYZ file |
| `xyz_paths.co2` | string | Path to FLP-CO₂ adduct XYZ file |
| `graphs_row_source` | string | Path to CSV file with energy data |
| `html_page` | string | Path to HTML description page |

**Note**: All paths are relative to the project root.

### Structure

| Field | Type | Description |
|-------|------|-------------|
| `la_lb_distance_A` | float or null | Distance between Lewis acid and Lewis base in Angstroms |
| `smiles.flp_bare` | object | SMILES and identifiers for bare FLP |
| `smiles.co2_adducts` | list | List of SMILES for CO₂ adducts (may have multiple binding modes) |

**SMILES Object**:
```yaml
rdkit: "..."           # Canonical SMILES from RDKit
openbabel: "..."       # SMILES from OpenBabel
inchi: "..."           # InChI identifier
inchikey: "..."        # InChIKey (hash)
validation:            # Round-trip validation results
  attempted: bool
  passed: bool
  rmsd: float          # RMSD in Angstroms
  error: string        # Error message if failed
```

### Energies Recovered

| Field | Type | Description |
|-------|------|-------------|
| `units` | string | Energy units (typically "kcal/mol") |
| `gas` | object | Gas-phase energies {E, H, G} |
| `solution` | object | Solution-phase energies {E, H, G} |

**Energy Types**:
- `E`: Electronic energy
- `H`: Enthalpy
- `G`: Gibbs free energy

**Note**: These are binding energies from the original FLPDB data. Students will recompute these.

### QC Flags

| Field | Type | Description |
|-------|------|-------------|
| `has_xyz_flp` | boolean | Bare FLP XYZ file exists |
| `has_xyz_co2` | boolean | CO₂ adduct XYZ file exists |
| `has_recovered_energy` | boolean | Original energy data available |
| `join_ok` | boolean | Entry has sufficient data for analysis |
| `smiles_validation_passed` | boolean or null | SMILES round-trip validation passed |

## Unit Conventions

- **Distances**: Angstroms (Å)
- **Energies**: kcal/mol (as reported in original FLPDB)
- **Temperature**: Kelvin (K)
- **Coordinates**: XYZ files in Angstroms

## Data Provenance

All data originates from:
- **Source**: FLPDB database snapshot (FLPDB-d549732.zip)
- **Original publication**: See `reference/papers/`
- **XYZ geometries**: Optimized structures from original study
- **Energies**: DFT calculations reported in original paper

## Validation Rules

### SMILES Round-Trip Validation

1. XYZ file → OpenBabel → SMILES
2. SMILES → RDKit → 3D coordinates
3. Compare original vs regenerated coordinates (RMSD)
4. **Pass threshold**: RMSD ≤ 1.0 Å

**Interpretation**:
- `rmsd < 0.5 Å`: Excellent match
- `0.5 ≤ rmsd ≤ 1.0 Å`: Good match
- `rmsd > 1.0 Å`: Failed (possible bonding ambiguity)

### Data Completeness

An entry is considered **complete** if:
- `join_ok = true`
- At least one of: `has_xyz_flp` or `has_recovered_energy`

Expected complete entries: **132**

## Accessing the Registry

### Python

```python
import yaml
from pathlib import Path

# Load entry
with open("data/processed/entries/108.yaml") as f:
    entry = yaml.safe_load(f)

# Access fields
flp_id = entry["flp_id"]
smiles = entry["structure"]["smiles"]["flp_bare"]["rdkit"]
has_co2 = entry["qc_flags"]["has_xyz_co2"]

# Get XYZ file path
xyz_path = Path(entry["provenance"]["xyz_paths"]["flp"])
```

### Command Line

```bash
# View full entry
flpco2 inspect 108

# Export to CSV
flpco2 export -o data.csv --format csv

# Get statistics
flpco2 stats
```

## Schema Version History

- **v0.1.0** (2025-11-06): Initial schema with basic provenance and SMILES
- Future versions will add:
  - Computed descriptors
  - Recomputed energies with different methods
  - Additional binding modes
  - Stereochemistry annotations

