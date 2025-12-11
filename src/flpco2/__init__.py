"""
FLPCO2DB - Frustrated Lewis Pairs CO₂ Database

A curated dataset of FLP-CO₂ interactions with computational chemistry workflows.
"""

__version__ = "0.1.0"
__author__ = "MML Course Team"

from pathlib import Path

# Project root directory
PROJECT_ROOT = Path(__file__).parent.parent.parent

# Data directories
DATA_DIR = PROJECT_ROOT / "data"
RAW_DATA_DIR = DATA_DIR / "raw"
PROCESSED_DATA_DIR = DATA_DIR / "processed"
REGISTRY_PATH = PROCESSED_DATA_DIR / "co2_registry.yaml"
ENTRIES_DIR = PROCESSED_DATA_DIR / "entries"

# Reference directories
REFERENCE_DIR = PROJECT_ROOT / "reference"

__all__ = [
    "__version__",
    "PROJECT_ROOT",
    "DATA_DIR",
    "RAW_DATA_DIR",
    "PROCESSED_DATA_DIR",
    "REGISTRY_PATH",
    "ENTRIES_DIR",
    "REFERENCE_DIR",
]

