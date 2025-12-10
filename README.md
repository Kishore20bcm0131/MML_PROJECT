# FLPCO2DB - Frustrated Lewis Pairs CO₂ Database

[![Python](https://img.shields.io/badge/Python-3.11%2B-blue)](https://www.python.org/)
[![License](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)

A curated dataset of Frustrated Lewis Pairs (FLPs) and their interactions with CO₂, derived from the FLPDB database. This repository provides organized data, molecular identifiers, and computational chemistry workflows for Machine Learning in Molecular Sciences (MML) course projects.

## 📋 Overview

**FLPCO2DB** is a student-ready research dataset that includes:

- **133 unique FLP systems** with complete structural data
- **927 XYZ geometry files** (bare FLPs and various small molecule adducts)
- **132 FLP-CO₂ complexes** with binding geometries
- **Recovered thermochemical data** from original DFT calculations
- **SMILES, InChI, InChIKey** molecular identifiers
- **Full provenance tracking** linking all data to sources
- **Reference computational workflows** using autoDE, ORCA, and xTB
- **Extensible registry system** for adding your own calculations

### What Makes This Dataset Unique?

- ✅ **Organized and queryable**: Structured YAML registry with CLI tools
- ✅ **Validated identifiers**: SMILES with round-trip validation
- ✅ **Complete provenance**: Every data point linked to source files
- ✅ **Ready for ML**: Export to CSV/JSON for machine learning pipelines
- ✅ **Extensible**: Framework for adding your own computational results
- ✅ **Educational**: Reference materials and worked examples included

## 🚀 Quick Start

### Prerequisites

- Python 3.9 or higher
- Conda package manager
- OpenBabel (for SMILES generation)
- RDKit (for cheminformatics)

### Installation

```bash
# Clone the repository
git clone https://github.com/digitalmoleculardesign/FLPCO2DB.git
cd FLPCO2DB

# Option A: Use existing mml_comp_chem environment
conda activate mml_comp_chem
pip install -e .

# Option B: Create new environment
conda env create -f env/environment.yml
conda activate flpco2
pip install -e .

# Verify installation
flpco2 --help
```

### Run the Pipeline

```bash
# Full automated pipeline
make all

# Or step-by-step:
make stage      # Extract and verify raw data
make build      # Build registry with SMILES
make validate   # Validate registry integrity
make stats      # Show dataset statistics
```

### Expected Output

After running `make all`, you should see:

```
✓ Data staging verified successfully!
✓ Registry built successfully!
  Entries: 133
  With CO₂ XYZ: 132
  With CO₂ energy: 137
  Overlap: 132
  SMILES validated: ~120
✓ Registry validation passed!
```

## 📓 Student Starter Notebook

For students working on MML course projects, we provide a comprehensive starter notebook:

```bash
cd notebooks/
jupyter notebook flpco2_project_starter.ipynb
```

**What's included:**
- **Part 1: Data Exploration** - Load and explore FLPCO2DB registry
- **Part 2: Molecular Parametrization** - Morgan fingerprints + QM descriptors
- **Part 3: ML Models** - Ridge, Lasso, Bayesian Ridge, Random Forest
- **Part 4: Model Interpretation** - Feature importance and chemical insights
- **Part 5: Candidate Screening** - Rank FLPs for experimental validation
- **Part 6: Project Milestones** - Timeline and deliverables

**6 hands-on exercises** for students to implement key functionality:
- Exercise 1.1: Build dataset DataFrame
- Exercise 2.1: Complete fingerprint pipeline
- Exercise 2.2: Calculate QM descriptors
- Exercise 3.1: Complete ML training pipeline
- Exercise 4.1: Chemical interpretation
- Exercise 5.1: Rank FLP candidates

**Helper utilities in `notebooks/utils.py`:**
- `load_flpco2_registry()` - Load central registry
- `load_flp_entry()` - Load individual entries
- `generate_morgan_fingerprint()` - Fingerprint generation
- `plot_parity()`, `plot_residuals()` - Visualization tools
- `compute_ml_metrics()` - MAE, RMSE, R² metrics

See [`notebooks/README.md`](notebooks/README.md) for complete documentation and workflow tips.

## 📊 Dataset Overview

### Data Coverage

| Category | Count | Description |
|----------|-------|-------------|
| **Total FLPs** | 133 | Unique frustrated Lewis pair systems |
| **XYZ Files** | 927 | Total geometry files (all substrates) |
| **CO₂ Adducts** | 132 | FLP-CO₂ binding geometries |
| **Energy Data** | 137 | Thermochemical data from original study |
| **Complete Entries** | 132 | Entries with both geometry and energy |

### Data Sources

- **Original Database**: FLPDB (Frustrated Lewis Pairs Database)
- **Publication**: *The First Frustrated Lewis Pairs Database* (see `reference/papers/`)
- **Methods**: DFT calculations with various functionals
- **Substrates**: CO₂, H₂, H₂O, CH₃OH, H₂CO, HCOOH

### File Structure

```
FLPCO2DB/
├── data/
│   ├── raw/                     # Original FLPDB data (immutable)
│   │   ├── xyz/                # 133 FLP directories, 927 XYZ files
│   │   ├── graphs_csv/         # CO2.csv with energies
│   │   └── html_pages/         # HTML descriptions
│   └── processed/              # Generated registry
│       ├── co2_registry.yaml   # Central registry
│       └── entries/            # Per-FLP YAML files (133)
├── docs/                        # Complete documentation
│   ├── 00_README.md            # Getting started guide
│   ├── 01_REGISTRY_SCHEMA.md   # Data format documentation
│   ├── 02_COMPUTE_PROTOCOL.md  # Calculation guidance
│   ├── 03_PIPELINE.md          # Pipeline details
│   └── 04_EXAMPLES.md          # Usage examples
├── reference/                   # Reference materials
│   ├── mml_studio_07/          # Example workflows
│   │   ├── mml_studio_07.ipynb # Jupyter notebook examples
│   │   └── utils.py            # Computational chemistry utilities
│   ├── papers/                 # Original research papers
│   └── original_csvs/          # Pre-extracted CSV files
├── src/flpco2/                 # Python package
│   ├── cli.py                  # Command-line interface
│   ├── staging.py              # Data staging utilities
│   ├── smiles_utils.py         # SMILES generation & validation
│   └── registry_builder.py     # Registry generation
├── tests/                       # Test suite
├── notebooks/                   # Student project materials
│   ├── flpco2_project_starter.ipynb  # Starter notebook with exercises
│   ├── utils.py                # Helper functions for ML workflows
│   └── README.md               # Notebook documentation
├── Makefile                     # Automation targets
└── README.md                    # This file
```

## 💻 Usage

### Command-Line Interface

```bash
# Data management
flpco2 stage                    # Stage and verify raw data
flpco2 build-reg                # Build registry with SMILES
flpco2 validate                 # Validate registry integrity

# Data exploration
flpco2 stats                    # Show dataset statistics
flpco2 inspect 108              # View details for FLP 108
flpco2 export -o data.csv       # Export to CSV

# Options
flpco2 build-reg --no-validate  # Skip SMILES validation (faster)
flpco2 build-reg --max-entries 10  # Test with subset
```

### Python API

```python
import yaml
from pathlib import Path

# Load central registry
with open("data/processed/co2_registry.yaml") as f:
    registry = yaml.safe_load(f)

print(f"Total FLPs: {registry['counts']['flps_total']}")

# Load specific entry
with open("data/processed/entries/108.yaml") as f:
    entry = yaml.safe_load(f)

print(f"FLP {entry['flp_id']}")
print(f"SMILES: {entry['structure']['smiles']['flp_bare']['rdkit']}")
print(f"XYZ path: {entry['provenance']['xyz_paths']['flp']}")
```

### With RDKit

```python
from rdkit import Chem
from rdkit.Chem import Descriptors
import yaml

# Load entry and get SMILES
with open("data/processed/entries/108.yaml") as f:
    entry = yaml.safe_load(f)

smiles = entry['structure']['smiles']['flp_bare']['rdkit']
mol = Chem.MolFromSmiles(smiles)

# Calculate properties
mw = Descriptors.MolWt(mol)
logp = Descriptors.MolLogP(mol)
print(f"MW: {mw:.2f}, LogP: {logp:.2f}")
```

### Export for Machine Learning

```bash
# Export to CSV
flpco2 export -o flp_data.csv --format csv

# Use in your ML pipeline
python
>>> import pandas as pd
>>> df = pd.read_csv("flp_data.csv")
>>> df.head()
```

## 🔬 Running Calculations

This repository provides the **data infrastructure only**. Students implement their own calculations using the reference materials.

### Reference Workflow

See `reference/mml_studio_07/` for complete examples using:
- **autoDE**: Automated reaction profile generation
- **ORCA**: DFT calculations
- **xTB**: Fast semi-empirical methods
- **py3Dmol**: 3D molecular visualization
- **NBO**: Natural bond orbital analysis

### Example Calculation Pattern

```python
import autode as ade
from pathlib import Path
import yaml

# Load FLP from registry
with open("data/processed/entries/108.yaml") as f:
    entry = yaml.safe_load(f)

xyz_path = Path(entry['provenance']['xyz_paths']['flp'])

# Create autoDE molecule
mol = ade.Molecule(str(xyz_path))

# Optimize with xTB (fast)
mol.optimise(method=ade.methods.XTB())

# Save result
mol.print_xyz_file("optimized_108.xyz")

# Your analysis here...
```

For complete workflows, see:
- `docs/02_COMPUTE_PROTOCOL.md` - Calculation guidance
- `reference/mml_studio_07/mml_studio_07.ipynb` - Worked examples
- `reference/mml_studio_07/utils.py` - Utility functions

## 📚 Documentation

| Document | Description |
|----------|-------------|
| [`docs/00_README.md`](docs/00_README.md) | Getting started guide |
| [`docs/01_REGISTRY_SCHEMA.md`](docs/01_REGISTRY_SCHEMA.md) | Registry format and fields |
| [`docs/02_COMPUTE_PROTOCOL.md`](docs/02_COMPUTE_PROTOCOL.md) | Computational workflow guidance |
| [`docs/03_PIPELINE.md`](docs/03_PIPELINE.md) | Data pipeline details |
| [`docs/04_EXAMPLES.md`](docs/04_EXAMPLES.md) | Code examples and recipes |
| [`notebooks/README.md`](notebooks/README.md) | Student notebook guide and workflow tips |

## 🎯 Use Cases

### For Students

- **MML Course Projects**: Ready-to-use dataset with starter notebook (`notebooks/flpco2_project_starter.ipynb`)
- **Hands-on Exercises**: 6 exercises covering data exploration, feature engineering, ML models, and interpretation
- **Learning Computational Chemistry**: Reference workflows with autoDE, ORCA, xTB
- **Data Science Practice**: Export to CSV for pandas/scikit-learn workflows
- **Cheminformatics**: Work with SMILES, InChI, molecular descriptors
- **Complete Examples**: Utility functions and plotting tools in `notebooks/utils.py`

### For Researchers

- **FLP Screening**: Identify promising FLP candidates for CO₂ capture
- **Method Benchmarking**: Compare DFT methods on standardized set
- **ML Model Training**: Features + targets for binding energy prediction
- **Descriptor Development**: Test new molecular descriptors on real systems

### Example Projects

1. **Predict CO₂ binding energies** using ML (SMILES → ΔG)
2. **Identify key structural features** controlling FLP reactivity
3. **Benchmark DFT functionals** against original data
4. **Develop new descriptors** for Lewis acidity/basicity
5. **Screen FLP libraries** for optimal CO₂ binding

## 🧪 Registry System

### What's in a Registry Entry?

Each FLP entry includes:

```yaml
flp_id: 108                          # Unique identifier
flp_code: "0BN08011"                 # Original code
provenance:
  xyz_paths:
    flp: "data/raw/xyz/108/108.xyz"  # Bare FLP geometry
    co2: "data/raw/xyz/108/108CO2.xyz"  # CO₂ adduct
structure:
  smiles:
    flp_bare:
      rdkit: "..."                   # Canonical SMILES
      inchi: "..."                   # InChI identifier
      inchikey: "..."                # InChIKey hash
energies_recovered:
  units: "kcal/mol"
  solution:
    G: 6.70                          # Binding free energy
qc_flags:
  has_xyz_co2: true
  has_recovered_energy: true
  smiles_validation_passed: true
```

See `docs/01_REGISTRY_SCHEMA.md` for complete schema documentation.

## 🔧 Development

### Running Tests

```bash
# Install development dependencies
pip install -e ".[dev]"

# Run tests
pytest tests/ -v

# With coverage
pytest tests/ --cov=src/flpco2 --cov-report=html
```

### Code Quality

```bash
# Format code
make format

# Lint code
make lint

# Type checking (if using mypy)
mypy src/flpco2/
```

### Building Documentation

Documentation is in Markdown format in `docs/`. To view locally:

```bash
# Serve with any markdown viewer or:
grip README.md  # GitHub-flavored preview
```

## 🤝 Contributing

This is a course project repository. Students are encouraged to:

1. **Add computational results**: Document your calculations in entry files
2. **Implement descriptors**: Add new molecular features
3. **Improve documentation**: Fix errors, add examples
4. **Report issues**: Found a bug? Open an issue!

See `CONTRIBUTING.md` for guidelines (if created).

## 📖 Citation

If you use this dataset, please cite:

1. **This repository**:
   ```
   FLPCO2DB: Frustrated Lewis Pairs CO₂ Database
   https://github.com/digitalmoleculardesign/FLPCO2DB
   ```

2. **Original FLPDB paper** (see `reference/papers/`):
   ```
   [Citation from original paper]
   ```

## 📄 License

[Specify license - typically MIT or CC-BY-4.0 for datasets]

## 🙏 Acknowledgments

- Original FLPDB authors and database maintainers
- MML Course instructors and teaching team
- autoDE, RDKit, OpenBabel, and ORCA developers

## 📞 Support

- **Documentation**: See `docs/` directory
- **Issues**: GitHub Issues tab
- **Examples**: `docs/04_EXAMPLES.md`
- **Reference code**: `reference/mml_studio_07/`

## 🗺️ Roadmap

Completed:

- [x] Jupyter notebook tutorials (`notebooks/flpco2_project_starter.ipynb`)
- [x] ML model training examples (included in starter notebook)
- [x] Helper utility functions (`notebooks/utils.py`)

Potential future enhancements:

- [ ] Automated descriptor calculation pipeline
- [ ] Web interface for registry browsing
- [ ] Integration with molecular databases (PubChem, ChEMBL)
- [ ] Additional substrate reactions (H₂, H₂O, etc.)
- [ ] Database backend (SQLite/PostgreSQL)
- [ ] Advanced ML examples (neural networks, graph convolutions)

## ⚡ Performance Notes

- **Registry building**: ~2-5 min (without validation), ~10-20 min (with validation)
- **SMILES validation**: 3D coordinate generation is slow; use `--no-validate` for testing
- **Data size**: ~50 MB total (XYZ files + registry)

## 🔗 Related Resources

- **autoDE**: https://duartegroup.github.io/autodE/
- **RDKit**: https://www.rdkit.org/
- **ORCA**: https://orcaforum.kofo.mpg.de/
- **xTB**: https://xtb-docs.readthedocs.io/
- **FLP Chemistry**: See `reference/papers/`

---

**Made with ❤️ for the MML Course**

*For questions about the course or dataset, contact your instructors.*

