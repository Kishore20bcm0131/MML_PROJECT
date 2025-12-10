"""
Basic tests for FLPCO2DB registry functionality.
"""

import pytest
from pathlib import Path
import yaml
from flpco2 import REGISTRY_PATH, ENTRIES_DIR, RAW_DATA_DIR


class TestRegistry:
    """Test registry integrity and structure."""
    
    def test_registry_exists(self):
        """Test that central registry file exists."""
        if not REGISTRY_PATH.exists():
            pytest.skip("Registry not built yet. Run: flpco2 build-reg")
        
        assert REGISTRY_PATH.exists(), "Registry file should exist"
        assert REGISTRY_PATH.suffix == ".yaml", "Registry should be YAML format"
    
    def test_registry_valid_yaml(self):
        """Test that registry is valid YAML."""
        if not REGISTRY_PATH.exists():
            pytest.skip("Registry not built yet")
        
        with open(REGISTRY_PATH, 'r') as f:
            registry = yaml.safe_load(f)
        
        assert isinstance(registry, dict), "Registry should be a dictionary"
    
    def test_registry_required_fields(self):
        """Test that registry has required fields."""
        if not REGISTRY_PATH.exists():
            pytest.skip("Registry not built yet")
        
        with open(REGISTRY_PATH, 'r') as f:
            registry = yaml.safe_load(f)
        
        required_fields = ['generated_at', 'dataset_version', 'counts', 'entries']
        for field in required_fields:
            assert field in registry, f"Registry should have '{field}' field"
    
    def test_registry_counts(self):
        """Test that registry counts are reasonable."""
        if not REGISTRY_PATH.exists():
            pytest.skip("Registry not built yet")
        
        with open(REGISTRY_PATH, 'r') as f:
            registry = yaml.safe_load(f)
        
        counts = registry['counts']
        
        # Expected counts from prompt
        assert counts['flps_total'] >= 130, "Should have at least 130 FLPs"
        assert counts['with_xyz_co2'] >= 130, "Should have at least 130 CO₂ XYZ files"
        assert counts['overlap'] >= 130, "Should have at least 130 overlapping entries"
    
    def test_entry_files_exist(self):
        """Test that referenced entry files exist."""
        if not REGISTRY_PATH.exists():
            pytest.skip("Registry not built yet")
        
        with open(REGISTRY_PATH, 'r') as f:
            registry = yaml.safe_load(f)
        
        entries = registry.get('entries', [])
        assert len(entries) > 0, "Registry should have entries"
        
        # Check first few entries
        for entry_ref in entries[:5]:
            entry_file = ENTRIES_DIR / f"{entry_ref['flp_id']}.yaml"
            assert entry_file.exists(), f"Entry file should exist: {entry_file}"
    
    def test_entry_structure(self):
        """Test that individual entries have correct structure."""
        if not REGISTRY_PATH.exists():
            pytest.skip("Registry not built yet")
        
        with open(REGISTRY_PATH, 'r') as f:
            registry = yaml.safe_load(f)
        
        entries = registry.get('entries', [])
        if not entries:
            pytest.skip("No entries in registry")
        
        # Load first entry
        entry_ref = entries[0]
        entry_file = ENTRIES_DIR / f"{entry_ref['flp_id']}.yaml"
        
        with open(entry_file, 'r') as f:
            entry = yaml.safe_load(f)
        
        # Required top-level fields
        required_fields = ['flp_id', 'provenance', 'structure', 'qc_flags']
        for field in required_fields:
            assert field in entry, f"Entry should have '{field}' field"
        
        # Check QC flags
        assert isinstance(entry['qc_flags'], dict), "qc_flags should be a dictionary"
        assert 'has_xyz_flp' in entry['qc_flags'], "Should have has_xyz_flp flag"


class TestRawData:
    """Test raw data structure."""
    
    def test_raw_data_dir_exists(self):
        """Test that raw data directory exists."""
        assert RAW_DATA_DIR.exists(), "Raw data directory should exist"
    
    def test_xyz_directory_exists(self):
        """Test that XYZ directory exists."""
        xyz_dir = RAW_DATA_DIR / "xyz"
        assert xyz_dir.exists(), "XYZ directory should exist"
    
    def test_csv_files_exist(self):
        """Test that CSV files exist."""
        graphs_csv_dir = RAW_DATA_DIR / "graphs_csv"
        if not graphs_csv_dir.exists():
            pytest.skip("graphs_csv directory not found")
        
        co2_csv = graphs_csv_dir / "CO2.csv"
        assert co2_csv.exists(), "CO2.csv should exist"
    
    def test_xyz_count(self):
        """Test that expected number of XYZ files exist."""
        xyz_dir = RAW_DATA_DIR / "xyz"
        if not xyz_dir.exists():
            pytest.skip("XYZ directory not found")
        
        # Count all XYZ files
        xyz_files = list(xyz_dir.glob("*/*.xyz"))
        assert len(xyz_files) >= 900, f"Should have ~927 XYZ files, found {len(xyz_files)}"
    
    def test_co2_xyz_count(self):
        """Test that expected number of CO2 XYZ files exist."""
        xyz_dir = RAW_DATA_DIR / "xyz"
        if not xyz_dir.exists():
            pytest.skip("XYZ directory not found")
        
        # Count CO2 XYZ files
        co2_files = list(xyz_dir.glob("*/*CO2.xyz"))
        assert len(co2_files) >= 130, f"Should have ~132 CO2 XYZ files, found {len(co2_files)}"


@pytest.fixture
def sample_entry_id():
    """Provide a sample FLP ID for testing."""
    return 108


class TestEntryLoading:
    """Test loading and parsing individual entries."""
    
    def test_load_entry(self, sample_entry_id):
        """Test loading a specific entry."""
        entry_file = ENTRIES_DIR / f"{sample_entry_id}.yaml"
        
        if not entry_file.exists():
            pytest.skip(f"Entry {sample_entry_id} not found")
        
        with open(entry_file, 'r') as f:
            entry = yaml.safe_load(f)
        
        assert entry['flp_id'] == sample_entry_id
        assert 'structure' in entry
        assert 'qc_flags' in entry
    
    def test_smiles_present(self, sample_entry_id):
        """Test that SMILES are present if entry has XYZ."""
        entry_file = ENTRIES_DIR / f"{sample_entry_id}.yaml"
        
        if not entry_file.exists():
            pytest.skip(f"Entry {sample_entry_id} not found")
        
        with open(entry_file, 'r') as f:
            entry = yaml.safe_load(f)
        
        if entry['qc_flags']['has_xyz_flp']:
            smiles = entry['structure']['smiles'].get('flp_bare', {})
            assert smiles.get('rdkit') or smiles.get('openbabel'), \
                "Should have SMILES if XYZ file exists"


if __name__ == "__main__":
    pytest.main([__file__, "-v"])

