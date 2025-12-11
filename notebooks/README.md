# FLPCO2DB Project Notebooks

This directory contains Jupyter notebooks for the FLPCO2DB machine learning project.

## Getting Started

### 1. Activate the Environment

```bash
conda activate mml_comp_chem
```

### 2. Launch Jupyter

```bash
cd notebooks/
jupyter notebook
```

### 3. Start with the Starter Notebook

Open `flpco2_project_starter.ipynb` to begin your project.

## Notebooks

### `flpco2_project_starter.ipynb`

**Purpose:** Introduction and scaffolding for the FLPCO2DB project

**Contents:**
- Project overview and objectives
- Data loading and exploration
- Feature engineering (fingerprints + QM descriptors)
- ML model training and evaluation
- Model interpretation and candidate screening
- Project milestones and deliverables

**What's Provided:**
- Data loading examples
- Workflow templates
- Helper functions in `utils.py`
- References to `mml_studio_07` patterns

**What You Need to Implement:**
- Complete feature extraction pipelines
- Train and optimize ML models
- Interpret results chemically
- Rank FLP candidates for validation

### `utils.py`

**Purpose:** Utility functions for the project

**Key Functions:**
- `load_flpco2_registry()`: Load the central registry
- `load_flp_entry()`: Load individual FLP entries
- `read_xyz_file()`: Parse XYZ coordinates
- `generate_morgan_fingerprint()`: Create fingerprints from SMILES
- `visualize_morgan_bits()`: Visualize fingerprint fragments
- `plot_parity()`: Parity plots for model evaluation
- `plot_residuals()`: Residual analysis
- `plot_feature_importance()`: Feature importance visualization
- `compute_ml_metrics()`: MAE, RMSE, R²

## Reference Materials

### Studio 7 Notebooks

Located in `../reference/mml_studio_07/`:

- `mml_studio_07.ipynb`: Complete reference for:
  - Morgan fingerprint generation and visualization
  - QM descriptor calculation with XTB
  - Conformational ensemble handling
  - ML model training (Ridge, Lasso, Bayesian Ridge, Random Forest)
  - Model interpretation and feature importance
  - Uncertainty quantification
  - Case studies on molecular property prediction

**Use this as a reference for implementation patterns!**

### Project Proposal

`../reference/Project Proposal_MML.md` contains:
- Full project description
- Objectives and milestones
- Timeline
- Key references

## Data Structure

The FLPCO2DB data is organized as follows:

```
../data/
├── raw/                    # Original data from FLPDB
│   ├── xyz/               # XYZ coordinate files (133 FLPs)
│   ├── graphs_csv/        # Energy and metadata CSVs
│   └── html_pages/        # Source HTML pages
└── processed/             # Curated registry
    ├── co2_registry.yaml  # Central registry
    └── entries/           # Per-FLP YAML files (133 entries)
        ├── 1.yaml
        ├── 2.yaml
        └── ...
```

### Registry Structure

Each entry in `data/processed/entries/` contains:

```yaml
flp_id: 1
flp_code: "FLP1"
provenance:
  xyz_paths:
    flp: "path/to/bare_flp.xyz"
    co2: "path/to/co2_adduct.xyz"
  graphs_row_source: "path/to/CO2.csv"
  html_page: "path/to/1.html"
structure:
  la_lb_distance_A: null  # To be computed
  smiles:
    flp_bare:
      rdkit: "SMILES_STRING"
      inchi: "InChI_STRING"
      inchikey: "InChIKey"
energies_recovered:
  units: "kcal/mol"
  gas:
    E: -10.5
    H: -9.8
    G: -7.2
  solution:
    E: -12.1
    H: -11.4
    G: -9.3
qc_flags:
  has_xyz_flp: true
  has_xyz_co2: true
  has_recovered_energy: true
  join_ok: true
  smiles_validation_passed: null
```

## CLI Tools

The FLPCO2DB package provides command-line tools:

```bash
# View statistics
flpco2 stats

# Inspect a specific entry
flpco2 inspect 108

# Validate registry
flpco2 validate

# Export to CSV
flpco2 export --output flp_data.csv --format csv
```

## Workflow Tips

### 1. Start Simple

- Begin with fingerprints only
- Use Ridge regression as baseline
- Establish cross-validation framework
- Add QM features incrementally

### 2. Iterate Systematically

- Compare multiple models
- Track performance metrics
- Visualize predictions (parity plots)
- Analyze errors (residual plots)

### 3. Interpret Results

- Which features are important?
- Do results make chemical sense?
- Can you visualize important fragments?
- How do predictions compare to known FLP chemistry?

### 4. Handle Small Data

With only ~130 samples:
- Use cross-validation (not single train/test split)
- Prefer simpler models (avoid overfitting)
- Consider feature selection
- Leverage domain knowledge

## Getting Help

1. **mml_studio_07**: Reference notebook for implementation patterns
2. **utils.py**: Helper functions for common tasks
3. **CLI tools**: Explore data with `flpco2` commands
4. **Project proposal**: Check objectives and timeline
5. **Literature**: Review FLP chemistry papers in proposal

## Expected Deliverables

By the end of the project:

1. **Complete notebook** with documented workflow
2. **Trained ML models** with performance metrics
3. **Feature importance analysis** with chemical interpretation
4. **Ranked candidate list** for experimental validation
5. **Final presentation** with insights and recommendations

---

**Good luck! Remember: The goal is to discover design principles for CO₂ activation, not just achieve high R².** 🧪🤖
