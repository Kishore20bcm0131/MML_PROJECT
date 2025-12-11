# FLPCO2DB - Getting Started

Welcome to the **Frustrated Lewis Pairs CO₂ Database (FLPCO2DB)** project!

This document will help you get started with the dataset and tools.

## Overview

FLPCO2DB is a curated dataset of Frustrated Lewis Pairs (FLPs) and their interactions with CO₂, derived from the FLPDB database. This repository provides:

- **Organized raw data** from the FLPDB snapshot
- **Registry system** with SMILES, InChI, provenance tracking
- **CLI tools** for data management and exploration
- **Reference materials** including computational chemistry workflows
- **Framework** for students to perform their own calculations

## Quick Start

### 1. Environment Setup

You have two options:

**Option A: Use existing `mml_comp_chem` environment**
```bash
conda activate mml_comp_chem
cd FLPCO2DB
pip install -e .
```

**Option B: Create new `flpco2` environment**
```bash
cd FLPCO2DB
conda env create -f env/environment.yml
conda activate flpco2
pip install -e .
```

**Note**: Python 3.9+ is required (your current Python 3.9.23 works fine!)

### 2. Verify Installation

```bash
flpco2 --help
```

You should see the list of available commands.

### 3. Run the Pipeline

```bash
# Stage and verify raw data
make stage

# Build the registry
make build

# View statistics
make stats

# Validate everything
make validate
```

Or run all at once:
```bash
make all
```

## What's in the Repository?

### Data Structure

```
data/
├── raw/                    # Original FLPDB data (immutable)
│   ├── xyz/               # 133 FLP directories, 927 XYZ files
│   ├── graphs_csv/        # CO2.csv and other substrate CSVs
│   ├── html_pages/        # HTML descriptions
│   └── flpdb_zip/         # Extracted ZIP contents
└── processed/             # Generated registry
    ├── co2_registry.yaml  # Central registry
    └── entries/           # Per-FLP YAML files
```

### Reference Materials

```
reference/
├── mml_studio_07/         # Example notebook & utils.py
│   ├── mml_studio_07.ipynb
│   └── utils.py
├── papers/                # Original research papers
└── original_csvs/         # Pre-extracted CSV files
```

## CLI Commands

### Data Management

```bash
# Stage and verify raw data
flpco2 stage

# Build registry with SMILES generation
flpco2 build-reg

# Build without SMILES validation (faster)
flpco2 build-reg --no-validate

# Build only first 10 entries (testing)
flpco2 build-reg --max-entries 10
```

### Data Exploration

```bash
# Show dataset statistics
flpco2 stats

# Validate registry integrity
flpco2 validate

# Inspect specific FLP entry
flpco2 inspect 108

# Export to CSV for analysis
flpco2 export -o my_export.csv --format csv
```

## Expected Counts

After running `make all`, you should see:

- **133 total FLP entries**
- **132 with CO₂ XYZ files**
- **137 with CO₂ energy data** (from original paper)
- **132 overlap** (entries with both XYZ and energy)

## Using the Registry

### Python API

```python
import yaml
from pathlib import Path

# Load central registry
with open("data/processed/co2_registry.yaml", 'r') as f:
    registry = yaml.safe_load(f)

print(f"Total entries: {registry['counts']['flps_total']}")

# Load specific entry
with open("data/processed/entries/108.yaml", 'r') as f:
    entry = yaml.safe_load(f)

print(f"FLP {entry['flp_id']}")
print(f"SMILES: {entry['structure']['smiles']['flp_bare']['rdkit']}")
print(f"InChIKey: {entry['structure']['smiles']['flp_bare']['inchikey']}")
```

### Command Line

```bash
# View entry as YAML
flpco2 inspect 108

# View entry as JSON
flpco2 inspect 108 --format json
```

## Running Calculations (Future Work)

This repository provides the **data infrastructure** only. Students will implement calculations using:

1. **Reference workflow** in `reference/mml_studio_07/`
   - Example: autoDE, ORCA, xTB setup
   - Utilities for molecular visualization, NBO analysis
   
2. **Your own workflow**
   - Use the XYZ files in `data/raw/xyz/`
   - Follow patterns from `reference/mml_studio_07/utils.py`
   - Document your methods in the registry

## Troubleshooting

### OpenBabel not found

SMILES generation requires OpenBabel:

```bash
conda install -c conda-forge openbabel
```

### Registry build fails

Try building without validation first:

```bash
flpco2 build-reg --no-validate
```

Or test with a small subset:

```bash
flpco2 build-reg --max-entries 10
```

### Import errors

Make sure the package is installed:

```bash
cd FLPCO2DB
pip install -e .
```

## Next Steps

1. **Explore the data**: `make stats`, `flpco2 inspect 108`
2. **Review reference materials**: Check `reference/mml_studio_07/`
3. **Read documentation**: See other files in `docs/`
4. **Plan calculations**: Read `02_COMPUTE_PROTOCOL.md`
5. **Implement workflows**: Use patterns from `mml_studio_07.ipynb`

## Getting Help

- Check other documentation in `docs/`
- Review the reference notebook: `reference/mml_studio_07/mml_studio_07.ipynb`
- Read the original paper: `reference/papers/`
- Inspect registry schema: `docs/01_REGISTRY_SCHEMA.md`

## References

- Original FLPDB paper (see `reference/papers/`)
- autoDE documentation: https://duartegroup.github.io/autodE/
- RDKit documentation: https://www.rdkit.org/docs/

