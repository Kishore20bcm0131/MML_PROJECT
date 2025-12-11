# Usage Examples

This document provides practical examples of working with FLPCO2DB.

## Command-Line Examples

### Basic Workflow

```bash
# 1. Activate environment
conda activate mml_comp_chem  # or flpco2

# 2. Navigate to project
cd FLPCO2DB

# 3. Run full pipeline
make all

# Or step by step:
make stage      # Stage and verify data
make build      # Build registry
make validate   # Validate registry
make stats      # Show statistics
```

### Inspecting Data

```bash
# Show dataset statistics
flpco2 stats

# Inspect specific FLP
flpco2 inspect 108

# Inspect as JSON
flpco2 inspect 108 --format json

# Validate registry
flpco2 validate

# Export to CSV
flpco2 export -o my_data.csv --format csv
```

### Testing and Development

```bash
# Build only first 10 entries (fast testing)
flpco2 build-reg --max-entries 10

# Build without SMILES validation (faster)
flpco2 build-reg --no-validate

# Clean and rebuild
make clean
make build
```

## Python API Examples

### Loading the Registry

```python
import yaml
from pathlib import Path

# Load central registry
registry_path = Path("data/processed/co2_registry.yaml")
with open(registry_path, 'r') as f:
    registry = yaml.safe_load(f)

# Show statistics
counts = registry['counts']
print(f"Total FLPs: {counts['flps_total']}")
print(f"With CO₂ XYZ: {counts['with_xyz_co2']}")
print(f"Overlap: {counts['overlap']}")

# List all entries
for entry_ref in registry['entries']:
    flp_id = entry_ref['flp_id']
    has_co2 = entry_ref['qc_flags']['has_xyz_co2']
    print(f"FLP {flp_id}: CO₂ XYZ = {has_co2}")
```

### Loading Individual Entries

```python
import yaml

# Load specific entry
with open("data/processed/entries/108.yaml", 'r') as f:
    entry = yaml.safe_load(f)

# Access data
flp_id = entry['flp_id']
flp_code = entry['flp_code']
smiles = entry['structure']['smiles']['flp_bare']['rdkit']
inchikey = entry['structure']['smiles']['flp_bare']['inchikey']

print(f"FLP {flp_id} ({flp_code})")
print(f"SMILES: {smiles}")
print(f"InChIKey: {inchikey}")

# Get XYZ file path
xyz_path = Path(entry['provenance']['xyz_paths']['flp'])
print(f"XYZ file: {xyz_path}")
```

### Batch Processing

```python
import yaml
from pathlib import Path
from tqdm import tqdm

# Load registry
with open("data/processed/co2_registry.yaml", 'r') as f:
    registry = yaml.safe_load(f)

# Process all entries with CO₂ data
results = []
for entry_ref in tqdm(registry['entries']):
    if not entry_ref['qc_flags']['has_xyz_co2']:
        continue
    
    # Load full entry
    entry_file = Path("data/processed") / entry_ref['entry_file']
    with open(entry_file, 'r') as f:
        entry = yaml.safe_load(f)
    
    # Extract data
    flp_id = entry['flp_id']
    smiles = entry['structure']['smiles']['flp_bare']['rdkit']
    
    # Your processing here...
    results.append({
        'flp_id': flp_id,
        'smiles': smiles,
    })

print(f"Processed {len(results)} entries")
```

### Working with RDKit

```python
from rdkit import Chem
from rdkit.Chem import Descriptors
import yaml

# Load entry
with open("data/processed/entries/108.yaml", 'r') as f:
    entry = yaml.safe_load(f)

# Get SMILES
smiles = entry['structure']['smiles']['flp_bare']['rdkit']

# Create RDKit molecule
mol = Chem.MolFromSmiles(smiles)

if mol is not None:
    # Calculate properties
    mw = Descriptors.MolWt(mol)
    logp = Descriptors.MolLogP(mol)
    tpsa = Descriptors.TPSA(mol)
    
    print(f"Molecular Weight: {mw:.2f}")
    print(f"LogP: {logp:.2f}")
    print(f"TPSA: {tpsa:.2f}")
    
    # Count atoms
    n_atoms = mol.GetNumAtoms()
    n_heavy = mol.GetNumHeavyAtoms()
    print(f"Total atoms: {n_atoms}, Heavy: {n_heavy}")
```

### Loading XYZ Files

```python
from pathlib import Path
import numpy as np

def read_xyz(xyz_path):
    """Read XYZ file and return symbols and coordinates."""
    with open(xyz_path, 'r') as f:
        lines = f.readlines()
    
    n_atoms = int(lines[0].strip())
    comment = lines[1].strip()
    
    symbols = []
    coords = []
    
    for line in lines[2:2+n_atoms]:
        parts = line.strip().split()
        symbols.append(parts[0])
        coords.append([float(parts[1]), float(parts[2]), float(parts[3])])
    
    return symbols, np.array(coords), comment

# Load from registry
import yaml
with open("data/processed/entries/108.yaml", 'r') as f:
    entry = yaml.safe_load(f)

xyz_path = Path(entry['provenance']['xyz_paths']['flp'])
symbols, coords, comment = read_xyz(xyz_path)

print(f"Loaded {len(symbols)} atoms")
print(f"Center of mass: {coords.mean(axis=0)}")
```

### Using with autoDE

```python
import autode as ade
from pathlib import Path
import yaml

# Load entry
with open("data/processed/entries/108.yaml", 'r') as f:
    entry = yaml.safe_load(f)

# Get XYZ path
xyz_path = Path(entry['provenance']['xyz_paths']['flp'])

# Create autoDE molecule
mol = ade.Molecule(str(xyz_path))

print(f"Molecule: {mol.name}")
print(f"Charge: {mol.charge}")
print(f"Multiplicity: {mol.mult}")
print(f"Number of atoms: {mol.n_atoms}")

# You can now use mol for calculations
# mol.optimise(method=ade.methods.XTB())
```

### Filtering Entries

```python
import yaml

# Load registry
with open("data/processed/co2_registry.yaml", 'r') as f:
    registry = yaml.safe_load(f)

# Filter: entries with XYZ and energy
complete_entries = [
    e for e in registry['entries']
    if e['qc_flags']['has_xyz_co2'] and e['qc_flags']['has_recovered_energy']
]

print(f"Complete entries: {len(complete_entries)}")

# Filter: entries with validated SMILES
validated_entries = [
    e for e in registry['entries']
    if e['qc_flags'].get('smiles_validation_passed') == True
]

print(f"Validated SMILES: {len(validated_entries)}")

# Filter by FLP code pattern
entries_with_code = [
    e for e in registry['entries']
    if e['flp_code'] and e['flp_code'].startswith('0BN')
]

print(f"FLPs starting with '0BN': {len(entries_with_code)}")
```

### Creating DataFrames for Analysis

```python
import pandas as pd
import yaml
from pathlib import Path

# Load all entries into DataFrame
rows = []

with open("data/processed/co2_registry.yaml", 'r') as f:
    registry = yaml.safe_load(f)

for entry_ref in registry['entries']:
    entry_file = Path("data/processed") / entry_ref['entry_file']
    with open(entry_file, 'r') as f:
        entry = yaml.safe_load(f)
    
    row = {
        'flp_id': entry['flp_id'],
        'flp_code': entry.get('flp_code'),
        'has_xyz_flp': entry['qc_flags']['has_xyz_flp'],
        'has_xyz_co2': entry['qc_flags']['has_xyz_co2'],
        'has_energy': entry['qc_flags']['has_recovered_energy'],
        'smiles': entry['structure']['smiles'].get('flp_bare', {}).get('rdkit'),
        'inchikey': entry['structure']['smiles'].get('flp_bare', {}).get('inchikey'),
    }
    
    # Add energies if available
    if entry['qc_flags']['has_recovered_energy']:
        sol_energies = entry['energies_recovered'].get('solution', {})
        row['G_sol'] = sol_energies.get('G')
    
    rows.append(row)

df = pd.DataFrame(rows)
print(df.head())
print(f"\nShape: {df.shape}")

# Save to CSV
df.to_csv("flp_data.csv", index=False)
```

## Jupyter Notebook Examples

### Minimal Notebook Template

```python
# Cell 1: Imports
import yaml
import pandas as pd
from pathlib import Path
from rdkit import Chem
from rdkit.Chem import Draw
import matplotlib.pyplot as plt

# Cell 2: Load registry
with open("data/processed/co2_registry.yaml", 'r') as f:
    registry = yaml.safe_load(f)

print(f"Total entries: {registry['counts']['flps_total']}")

# Cell 3: Load and visualize a molecule
with open("data/processed/entries/108.yaml", 'r') as f:
    entry = yaml.safe_load(f)

smiles = entry['structure']['smiles']['flp_bare']['rdkit']
mol = Chem.MolFromSmiles(smiles)

# Draw molecule
img = Draw.MolToImage(mol, size=(400, 400))
plt.imshow(img)
plt.axis('off')
plt.show()

# Cell 4: Create DataFrame and analyze
# ... (see DataFrame example above)
```

### 3D Visualization with py3Dmol

```python
import py3Dmol
from pathlib import Path
import yaml

# Load entry
with open("data/processed/entries/108.yaml", 'r') as f:
    entry = yaml.safe_load(f)

# Get XYZ path
xyz_path = Path(entry['provenance']['xyz_paths']['flp'])

# Create 3D view
view = py3Dmol.view(width=600, height=400)
with open(xyz_path, 'r') as f:
    xyz_data = f.read()

view.addModel(xyz_data, 'xyz')
view.setStyle({'stick': {}})
view.setBackgroundColor('white')
view.zoomTo()
view.show()
```

## Advanced Examples

### Compare Recovered vs Computed Energies

```python
import yaml
import matplotlib.pyplot as plt

# Assuming you've computed new energies and stored them
recovered = []
computed = []

for entry_ref in registry['entries'][:10]:  # Example: first 10
    entry_file = Path("data/processed") / entry_ref['entry_file']
    with open(entry_file, 'r') as f:
        entry = yaml.safe_load(f)
    
    if entry['qc_flags']['has_recovered_energy']:
        G_recovered = entry['energies_recovered']['solution'].get('G')
        # G_computed = entry['computed_energies']['solution'].get('G')  # If you added this
        
        # For example only:
        G_computed = G_recovered + np.random.normal(0, 2)  # Simulated
        
        if G_recovered and G_computed:
            recovered.append(G_recovered)
            computed.append(G_computed)

# Plot correlation
plt.figure(figsize=(8, 8))
plt.scatter(recovered, computed, alpha=0.6)
plt.plot([-20, 20], [-20, 20], 'r--', label='Perfect agreement')
plt.xlabel('Recovered Energy (kcal/mol)')
plt.ylabel('Computed Energy (kcal/mol)')
plt.title('Energy Comparison')
plt.legend()
plt.grid(True, alpha=0.3)
plt.show()
```

### Export Specific Subset

```python
import yaml
import pandas as pd
from pathlib import Path

# Load registry
with open("data/processed/co2_registry.yaml", 'r') as f:
    registry = yaml.safe_load(f)

# Filter: only entries with validated SMILES
filtered_rows = []

for entry_ref in registry['entries']:
    if not entry_ref['qc_flags'].get('smiles_validation_passed'):
        continue
    
    entry_file = Path("data/processed") / entry_ref['entry_file']
    with open(entry_file, 'r') as f:
        entry = yaml.safe_load(f)
    
    filtered_rows.append({
        'flp_id': entry['flp_id'],
        'smiles': entry['structure']['smiles']['flp_bare']['rdkit'],
        'inchikey': entry['structure']['smiles']['flp_bare']['inchikey'],
    })

df = pd.DataFrame(filtered_rows)
df.to_csv("validated_flps.csv", index=False)
print(f"Exported {len(df)} validated entries")
```

## Integration with mml_studio_07

```python
# Import utilities from reference
import sys
sys.path.append("reference/mml_studio_07")
from utils import MolTo3DView

# Load FLP XYZ
with open("data/processed/entries/108.yaml", 'r') as f:
    entry = yaml.safe_load(f)

xyz_path = entry['provenance']['xyz_paths']['flp']

# Visualize with utils.py function
viewer = MolTo3DView(xyz_path, size=(600, 400), style="stick")
viewer.show()
```

## Common Tasks

### Task: Get all SMILES

```bash
# Using CLI
flpco2 export -o smiles.csv
# Then filter in spreadsheet or:
# cut -d',' -f1,6 smiles.csv
```

### Task: Find entry by SMILES

```python
import yaml

target_smiles = "CC1=CC=C(C=C1)P(C2=CC=CC=C2)..."

with open("data/processed/co2_registry.yaml", 'r') as f:
    registry = yaml.safe_load(f)

for entry_ref in registry['entries']:
    entry_file = Path("data/processed") / entry_ref['entry_file']
    with open(entry_file, 'r') as f:
        entry = yaml.safe_load(f)
    
    smiles = entry['structure']['smiles'].get('flp_bare', {}).get('rdkit')
    if smiles == target_smiles:
        print(f"Found: FLP {entry['flp_id']}")
        break
```

### Task: Get statistics by subset

```python
import yaml

with open("data/processed/co2_registry.yaml", 'r') as f:
    registry = yaml.safe_load(f)

# Count by QC flag
counts = {
    'has_xyz_flp': 0,
    'has_xyz_co2': 0,
    'has_energy': 0,
    'smiles_validated': 0,
}

for entry_ref in registry['entries']:
    flags = entry_ref['qc_flags']
    if flags['has_xyz_flp']:
        counts['has_xyz_flp'] += 1
    if flags['has_xyz_co2']:
        counts['has_xyz_co2'] += 1
    if flags['has_recovered_energy']:
        counts['has_energy'] += 1
    if flags.get('smiles_validation_passed'):
        counts['smiles_validated'] += 1

for key, value in counts.items():
    print(f"{key}: {value}")
```

##See Also

- `docs/00_README.md` - Getting started
- `docs/01_REGISTRY_SCHEMA.md` - Data structure details
- `docs/02_COMPUTE_PROTOCOL.md` - Running calculations
- `docs/03_PIPELINE.md` - Pipeline details
- `reference/mml_studio_07/` - Working examples

