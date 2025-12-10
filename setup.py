"""
Setup script for FLPCO2DB package.
"""

from setuptools import setup, find_packages
from pathlib import Path

# Read README
readme_file = Path(__file__).parent / "README.md"
long_description = readme_file.read_text() if readme_file.exists() else ""

setup(
    name="flpco2",
    version="0.1.0",
    author="MML Course Team",
    author_email="",
    description="Frustrated Lewis Pairs CO₂ Database - Dataset and computational workflows",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/digitalmoleculardesign/FLPCO2DB",
    project_urls={
        "Bug Tracker": "https://github.com/digitalmoleculardesign/FLPCO2DB/issues",
        "Documentation": "https://github.com/digitalmoleculardesign/FLPCO2DB/tree/main/docs",
    },
    packages=find_packages(where="src"),
    package_dir={"": "src"},
    python_requires=">=3.9",
    install_requires=[
        "numpy",
        "scipy",
        "pandas",
        "rdkit",  # Note: install via conda
        "pyyaml",
        "pydantic",
        "typer",
        "rich",
        "tqdm",
    ],
    extras_require={
        "dev": [
            "pytest",
            "pytest-cov",
            "black",
            "isort",
            "flake8",
        ],
        "full": [
            "cctk",
            "autode",
            "pyscf",
            "datamol",
            "matplotlib",
            "seaborn",
            "py3dmol",
            "jupyter",
        ],
    },
    entry_points={
        "console_scripts": [
            "flpco2=flpco2.cli:main",
        ],
    },
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Chemistry",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Programming Language :: Python :: 3.12",
    ],
    keywords="chemistry computational-chemistry frustrated-lewis-pairs co2 database",
)

