# Contributing to FLPCO2DB

Thank you for your interest in contributing to the Frustrated Lewis Pairs CO₂ Database!

## For Students in the MML Course

This repository is designed for educational purposes. You are encouraged to:

### 1. Add Your Computational Results

After running calculations, you can update entry files:

```yaml
# In data/processed/entries/108.yaml, add:
computed_energies:
  method: "M06-2X/def2-TZVP"
  date: "2025-11-15"
  computed_by: "Your Name"
  solution:
    G_bind: -8.23
    G_flp: -1234.56
    G_co2: -78.90
    G_adduct: -1321.69
```

### 2. Implement New Features

- Add descriptor calculations in `src/flpco2/descriptors.py`
- Extend CLI with new commands
- Create analysis notebooks in `notebooks/`

### 3. Improve Documentation

- Fix typos or clarify instructions
- Add examples to `docs/04_EXAMPLES.md`
- Document your computational workflows

### 4. Report Issues

Found a bug? Data error? Open an issue with:
- Clear description
- Steps to reproduce
- Expected vs actual behavior
- Your environment (OS, Python version, etc.)

## Development Workflow

### Setup

```bash
# Clone and install
git clone https://github.com/digitalmoleculardesign/FLPCO2DB.git
cd FLPCO2DB
conda activate mml_comp_chem
pip install -e ".[dev]"
```

### Making Changes

1. **Create a branch**
   ```bash
   git checkout -b feature/your-feature-name
   ```

2. **Make your changes**
   - Follow existing code style
   - Add tests if applicable
   - Update documentation

3. **Test your changes**
   ```bash
   pytest tests/ -v
   make lint
   ```

4. **Commit with clear messages**
   ```bash
   git add .
   git commit -m "Add: brief description of changes"
   ```

5. **Push and create Pull Request**
   ```bash
   git push origin feature/your-feature-name
   ```

## Code Style

### Python

- Follow PEP 8 style guide
- Use type hints where appropriate
- Document functions with docstrings
- Maximum line length: 100 characters

```python
def my_function(param: str) -> Dict[str, any]:
    """
    Brief description.
    
    Args:
        param: Parameter description
        
    Returns:
        Dictionary with results
    """
    # Implementation
    pass
```

### Format Code

```bash
# Auto-format with black
black src/flpco2/

# Sort imports
isort src/flpco2/

# Or use Make target
make format
```

## Testing

### Running Tests

```bash
# All tests
pytest tests/ -v

# With coverage
pytest tests/ --cov=src/flpco2

# Specific test
pytest tests/test_registry.py -v
```

### Writing Tests

Place tests in `tests/` directory:

```python
# tests/test_my_feature.py
import pytest
from flpco2 import my_function

def test_my_function():
    result = my_function("test")
    assert result is not None
```

## Documentation

### Markdown Files

- Use clear headings
- Include code examples
- Add links to related docs
- Keep language simple and clear

### Docstrings

Use Google or NumPy style:

```python
def process_entry(entry_id: int) -> Dict[str, any]:
    """
    Process a single FLP entry.
    
    Args:
        entry_id: FLP identifier
        
    Returns:
        Dictionary with processed data
        
    Raises:
        ValueError: If entry_id not found
    """
    pass
```

## Commit Messages

Use clear, descriptive commit messages:

- **Add**: New feature or file
- **Fix**: Bug fix
- **Update**: Modify existing feature
- **Docs**: Documentation changes
- **Test**: Add or modify tests
- **Refactor**: Code restructuring
- **Style**: Formatting changes

Examples:
```
Add: SMILES validation with round-trip check
Fix: Handle missing XYZ files gracefully
Update: Improve registry building performance
Docs: Add examples for autoDE integration
```

## Pull Request Process

1. **Ensure tests pass**: All tests must pass before merging
2. **Update documentation**: Document new features
3. **Clear PR description**: Explain what and why
4. **Link issues**: Reference related issues with `#issue_number`
5. **Request review**: Tag course instructors if needed

## Adding Computational Workflows

If you implement a new calculation pipeline:

1. **Document the method** in `docs/02_COMPUTE_PROTOCOL.md`
2. **Provide example code** in `notebooks/`
3. **Update registry schema** if adding new fields
4. **Test on small subset** before running on full dataset

## Questions?

- Check existing documentation in `docs/`
- Review reference materials in `reference/mml_studio_07/`
- Open an issue with your question
- Ask course instructors

## Code of Conduct

### Be Respectful

- Welcoming and inclusive language
- Respect different viewpoints
- Gracefully accept constructive criticism
- Focus on what's best for the project and learning

### Academic Integrity

- Give credit for code/ideas from others
- Don't submit someone else's work as your own
- Collaborate, but understand what you submit
- Cite sources and references appropriately

## License

By contributing, you agree that your contributions will be licensed under the same license as the project.

## Recognition

Contributors will be acknowledged in:
- Repository contributors list
- Project documentation (if significant contributions)

Thank you for contributing to FLPCO2DB and advancing FLP research! 🎉

