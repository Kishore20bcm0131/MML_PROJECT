# Computational Protocol Documentation

This document provides guidance for students who want to perform calculations on the FLPCO2DB dataset.

## Overview

**Important**: This repository provides **data infrastructure only**. No calculations are performed automatically. Students will implement their own computational workflows using the reference materials provided.

## Reference Materials

### mml_studio_07 Notebook

Location: `reference/mml_studio_07/`

This folder contains:
- `mml_studio_07.ipynb`: Jupyter notebook with example workflows
- `utils.py`: Utility functions for computational chemistry

**Key features**:
- autoDE setup and configuration
- ORCA and xTB integration
- NBO analysis parsing
- Molecular visualization (py3Dmol)
- Plotting utilities

### Using the Reference Code

```python
# Import utilities from reference
import sys
sys.path.append("reference/mml_studio_07")
from utils import *

# Configure autoDE (example from utils.py)
import autode as ade
ade.Config.n_cores = 6
ade.Config.max_core = 2000

# Your calculations here...
```

## Suggested Computational Workflows

### Option 1: Follow mml_studio_07 Patterns

1. **Load FLP structure** from registry
2. **Use autoDE** for conformer search and optimization
3. **Configure backend** (ORCA, xTB, etc.)
4. **Run calculations** with configurable methods
5. **Store results** in `data/results/` or your own structure

### Option 2: Custom Workflow

Students are encouraged to:
- Use different quantum chemistry packages
- Try different levels of theory
- Implement novel descriptors
- Compare methods

## Proposed Calculation Strategy

### For FLP-CO₂ Binding

Based on the original FLPDB methodology, students might consider:

#### 1. Pre-optimization
- **Method**: xTB (GFN2-xTB) or similar fast method
- **Purpose**: Generate reasonable starting geometries
- **Tool**: Can use xTB directly or via autoDE

#### 2. Final Optimization + Frequencies
- **Method**: DFT (e.g., M06-2X/def2-TZVP) or other configurable backend
- **Purpose**: Accurate energies and thermochemistry
- **Options**:
  - ORCA (with or without UMA)
  - Gaussian
  - Other QC packages

#### 3. Solvent Model
- **Recommended**: SMD (or equivalent)
- **Solvent**: Toluene (as in original study)
- **Purpose**: Solution-phase energetics

#### 4. Thermochemistry
- **Temperature**: 298.15 K
- **Standard state**: 1 M (convert from 1 atm if needed)
- **Entropy**: Quasi-RRHO corrections recommended

### Binding Free Energy Formula

For CO₂ binding:

```
ΔG_bind^sol = G^sol(FLP-CO₂) - G^sol(FLP) - G^sol(CO₂)
```

**Standard state correction** (1 atm → 1 M):
```
ΔG_corr = RT ln(24.46) ≈ 1.89 kcal/mol at 298.15 K
```

Apply this correction to bring all species to 1 M standard state.

## Descriptor Calculations

Students may wish to calculate molecular descriptors for machine learning or analysis:

### Global Descriptors
- **Structural**: LA-LB distance, molecular weight, etc.
- **Electronic**: HOMO, LUMO, gap, dipole moment
- **Reactivity**: χ (electronegativity), η (hardness), S (softness), ω (electrophilicity)

### Local Descriptors
- **Charges**: CM5, NBO, Hirshfeld, etc.
- **Fukui indices**: f⁺, f⁻ (from N±1 calculations)
- **Local reactivity**: Local ω, local S

### Steric Descriptors
- **Buried volume**: %V_bur at various radii
- **Sterimol**: L, B1, B5 parameters
- **Tools**: SEQCROW, ChimeraX, or scripted alternatives

## autoDE Configuration

Example configuration based on `utils.py`:

```python
import autode as ade

# Basic configuration
ade.Config.n_cores = 6  # Adjust for your system
ade.Config.max_core = 2000  # MB

# Set computational methods
ade.Config.hcode = "orca"  # High-level: ORCA
ade.Config.lcode = "xtb"   # Low-level: xTB

# Define method keywords (example)
orca_method = ade.methods.ORCA()
orca_method.keywords.opt = ['M06-2X', 'def2-TZVP', 'TightOPT']
orca_method.keywords.sp = ['M06-2X', 'def2-TZVP']

# Solvent
from autode.solvent.solvents import toluene
```

## Data Storage Recommendations

### Suggested Structure

```
data/
├── calculations/
│   ├── <your_method_name>/
│   │   ├── flp_bare/
│   │   │   ├── 108/
│   │   │   │   ├── input.xyz
│   │   │   │   ├── output.out
│   │   │   │   └── optimized.xyz
│   │   ├── co2_adduct/
│   │   └── results_summary.yaml
```

### Update Registry with Results

Consider adding your results to entry files:

```yaml
compute_plan:
  status: "completed"
  method_profile: "my_DFT_profile_001"
  level_of_theory: "M06-2X/def2-TZVP"
  solvent: "toluene (SMD)"
  date_completed: "2025-11-15"
  
computed_energies:
  units: "kcal/mol"
  method: "M06-2X/def2-TZVP"
  solution:
    G_bind: -8.23
    G_flp: -1234.56
    G_co2: -78.90
    G_adduct: -1321.69
```

## Quality Control

### Check Calculations
- Verify convergence
- Check imaginary frequencies (should be 0 for minima)
- Compare with literature values
- Document all failures/issues

### Reproducibility
- Record all settings/versions
- Save input files
- Document workflow
- Version results

## Comparing with Original Data

Your computed values can be compared with recovered energies:

```python
import yaml

with open("data/processed/entries/108.yaml") as f:
    entry = yaml.safe_load(f)

original_G = entry["energies_recovered"]["solution"]["G"]
your_G = -8.23  # Your computed value

difference = abs(your_G - original_G)
print(f"Difference from original: {difference:.2f} kcal/mol")
```

**Expected agreement**: ±2-3 kcal/mol typical for DFT methods

## Performance Considerations

### Parallelization
- Use multiple cores efficiently
- Don't over-parallelize (diminishing returns)
- Monitor memory usage

### UMA (Universal Machine Learning Assisted) Methods
From memories, the project has experience with:
- ORCA+UMA integration
- GPU-accelerated calculations
- Reduced CPU/memory requirements

See reference materials for UMA configuration if available.

### Batch Processing
- Process multiple FLPs in parallel
- Use job submission systems if available
- Implement restart capability

## Example Workflows

### Minimal Example

```python
from pathlib import Path
import autode as ade

# Load FLP geometry
flp_xyz = Path("data/raw/xyz/108/108.xyz")
mol = ade.Molecule(str(flp_xyz))

# Optimize with xTB
mol.optimise(method=ade.methods.XTB())

# Save result
mol.print_xyz_file("optimized_108.xyz")
```

### Complete Example

See `reference/mml_studio_07/mml_studio_07.ipynb` for full workflows including:
- Conformer generation
- Multi-step optimization
- Frequency calculations
- NBO analysis
- Visualization

## Best Practices

1. **Start small**: Test workflow on 1-2 FLPs first
2. **Document everything**: Methods, settings, versions
3. **Version results**: Tag calculations with date/method
4. **Check convergence**: Always verify calculations succeeded
5. **Compare**: Cross-check with original FLPDB data
6. **Share**: Document your workflow for reproducibility

## Further Resources

- **autoDE docs**: https://duartegroup.github.io/autodE/
- **ORCA manual**: https://orcaforum.kofo.mpg.de/
- **xTB docs**: https://xtb-docs.readthedocs.io/
- **Original FLPDB paper**: See `reference/papers/`
- **mml_studio_07**: See `reference/mml_studio_07/`

## Questions?

- Review `reference/mml_studio_07/utils.py` for working examples
- Check the original FLPDB paper methodology
- Consult course materials and instructors

