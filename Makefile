.PHONY: help env install stage verify build stats validate export clean test lint format all

# Default target
help:
	@echo "FLPCO2DB Makefile"
	@echo "================="
	@echo ""
	@echo "Available targets:"
	@echo "  env        - Create conda environment"
	@echo "  install    - Install package in editable mode"
	@echo "  stage      - Stage and verify raw data"
	@echo "  verify     - Verify raw data structure"
	@echo "  build      - Build registry from raw data"
	@echo "  stats      - Show dataset statistics"
	@echo "  validate   - Validate registry"
	@echo "  export     - Export registry to CSV"
	@echo "  test       - Run tests"
	@echo "  lint       - Run linting"
	@echo "  format     - Format code with black and isort"
	@echo "  clean      - Clean generated files"
	@echo "  all        - Run full pipeline: stage, build, validate, stats"
	@echo ""
	@echo "Quick start:"
	@echo "  make env && make install && make all"

# Create conda environment
env:
	@echo "Creating conda environment..."
	conda env create -f env/environment.yml

# Install package in editable mode
install:
	@echo "Installing flpco2 package..."
	pip install -e .

# Install with development dependencies
install-dev:
	@echo "Installing flpco2 package with dev dependencies..."
	pip install -e ".[dev]"

# Stage and verify raw data
stage:
	@echo "Staging and verifying raw data..."
	flpco2 stage

# Verify data structure only
verify:
	@echo "Verifying raw data structure..."
	flpco2 stage

# Build registry
build:
	@echo "Building CO₂ registry..."
	flpco2 build-reg

# Build registry without SMILES validation (faster)
build-fast:
	@echo "Building CO₂ registry (no SMILES validation)..."
	flpco2 build-reg --no-validate

# Build registry for testing (first 10 entries)
build-test:
	@echo "Building test registry (10 entries)..."
	flpco2 build-reg --max-entries 10

# Show dataset statistics
stats:
	@echo "Dataset statistics:"
	flpco2 stats

# Validate registry
validate:
	@echo "Validating registry..."
	flpco2 validate

# Export registry to CSV
export:
	@echo "Exporting registry to data/processed/registry_export.csv..."
	flpco2 export -o data/processed/registry_export.csv --format csv

# Inspect a specific entry
inspect:
	@echo "Usage: make inspect FLP_ID=108"
	@echo "Inspecting FLP $(FLP_ID)..."
	flpco2 inspect $(FLP_ID)

# Run tests
test:
	@echo "Running tests..."
	pytest tests/ -v

# Run tests with coverage
test-cov:
	@echo "Running tests with coverage..."
	pytest tests/ -v --cov=src/flpco2 --cov-report=html

# Lint code
lint:
	@echo "Linting code..."
	flake8 src/flpco2/ --max-line-length=100
	@echo "✓ Linting passed"

# Format code
format:
	@echo "Formatting code..."
	black src/flpco2/
	isort src/flpco2/
	@echo "✓ Code formatted"

# Clean generated files
clean:
	@echo "Cleaning generated files..."
	rm -rf data/processed/
	rm -rf data/interim/
	rm -rf src/flpco2/__pycache__
	rm -rf tests/__pycache__
	rm -rf .pytest_cache
	rm -rf htmlcov
	rm -rf dist
	rm -rf build
	rm -rf src/*.egg-info
	find . -type f -name "*.pyc" -delete
	find . -type d -name "__pycache__" -delete
	@echo "✓ Cleaned"

# Full pipeline
all: stage build validate stats
	@echo ""
	@echo "✓ Full pipeline completed!"
	@echo ""
	@echo "Next steps:"
	@echo "  - Explore the registry: make stats"
	@echo "  - Inspect entries: make inspect FLP_ID=108"
	@echo "  - Export data: make export"
	@echo "  - See reference/mml_studio_07/ for calculation examples"

