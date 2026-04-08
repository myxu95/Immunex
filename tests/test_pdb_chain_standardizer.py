#!/usr/bin/env python3
"""
Unit tests for PDB Chain Standardizer module

Tests cover:
- Chain analysis functionality
- Chain mapping creation
- PDB file standardization
- Validation methods
- Batch processing
- Error handling
"""

import pytest
import tempfile
import shutil
from pathlib import Path
from immunex.analysis.structure.pdb_chain_standardizer import (
    PDBChainStandardizer,
    ChainInfo,
    StandardizationResult
)


@pytest.fixture
def temp_dir():
    """Create temporary directory for test files."""
    tmp = tempfile.mkdtemp()
    yield Path(tmp)
    shutil.rmtree(tmp)


@pytest.fixture
def sample_pdb_5chain(temp_dir):
    """
    Create sample 5-chain PDB file (pHLA-TCR complex).
    Chain order: X(10 res) < Y(100 res) < Z(150 res) < W(220 res) < V(280 res)
    """
    pdb_content = ""
    atom_num = 1

    # Chain X: Peptide (10 residues)
    for res in range(1, 11):
        pdb_content += f"ATOM  {atom_num:5d}  CA  ALA X{res:4d}      10.000  10.000  10.000  1.00 20.00           C\n"
        atom_num += 1

    # Chain Y: HLA-β (100 residues)
    for res in range(1, 101):
        pdb_content += f"ATOM  {atom_num:5d}  CA  ALA Y{res:4d}      20.000  20.000  20.000  1.00 20.00           C\n"
        atom_num += 1

    # Chain Z: TCR-α (150 residues)
    for res in range(1, 151):
        pdb_content += f"ATOM  {atom_num:5d}  CA  ALA Z{res:4d}      30.000  30.000  30.000  1.00 20.00           C\n"
        atom_num += 1

    # Chain W: TCR-β (220 residues)
    for res in range(1, 221):
        pdb_content += f"ATOM  {atom_num:5d}  CA  ALA W{res:4d}      40.000  40.000  40.000  1.00 20.00           C\n"
        atom_num += 1

    # Chain V: HLA-α (280 residues)
    for res in range(1, 281):
        pdb_content += f"ATOM  {atom_num:5d}  CA  ALA V{res:4d}      50.000  50.000  50.000  1.00 20.00           C\n"
        atom_num += 1

    pdb_file = temp_dir / "sample_5chain.pdb"
    pdb_file.write_text(pdb_content)
    return pdb_file


@pytest.fixture
def sample_pdb_standard(temp_dir):
    """Create sample PDB already in standard order (C-B-D-E-A)."""
    pdb_content = ""
    atom_num = 1

    chains = [
        ('C', 10),   # Peptide
        ('B', 100),  # HLA-β
        ('D', 150),  # TCR-α
        ('E', 220),  # TCR-β
        ('A', 280)   # HLA-α
    ]

    for chain_id, num_res in chains:
        for res in range(1, num_res + 1):
            pdb_content += f"ATOM  {atom_num:5d}  CA  ALA {chain_id}{res:4d}      10.000  10.000  10.000  1.00 20.00           C\n"
            atom_num += 1

    pdb_file = temp_dir / "sample_standard.pdb"
    pdb_file.write_text(pdb_content)
    return pdb_file


@pytest.fixture
def sample_pdb_multichain(temp_dir):
    """Create sample PDB with too many chains (7 chains)."""
    pdb_content = ""
    atom_num = 1

    for chain_id in ['A', 'B', 'C', 'D', 'E', 'F', 'G']:
        for res in range(1, 51):
            pdb_content += f"ATOM  {atom_num:5d}  CA  ALA {chain_id}{res:4d}      10.000  10.000  10.000  1.00 20.00           C\n"
            atom_num += 1

    pdb_file = temp_dir / "sample_multichain.pdb"
    pdb_file.write_text(pdb_content)
    return pdb_file


@pytest.fixture
def sample_pdb_insufficient(temp_dir):
    """Create sample PDB with too few chains (3 chains)."""
    pdb_content = ""
    atom_num = 1

    for chain_id in ['A', 'B', 'C']:
        for res in range(1, 101):
            pdb_content += f"ATOM  {atom_num:5d}  CA  ALA {chain_id}{res:4d}      10.000  10.000  10.000  1.00 20.00           C\n"
            atom_num += 1

    pdb_file = temp_dir / "sample_insufficient.pdb"
    pdb_file.write_text(pdb_content)
    return pdb_file


class TestChainInfo:
    """Test ChainInfo dataclass."""

    def test_chain_info_creation(self):
        """Test creating ChainInfo object."""
        chain = ChainInfo(chain_id='A', residue_count=280, atom_count=2240)
        assert chain.chain_id == 'A'
        assert chain.residue_count == 280
        assert chain.atom_count == 2240

    def test_chain_info_repr(self):
        """Test ChainInfo string representation."""
        chain = ChainInfo(chain_id='C', residue_count=10, atom_count=80)
        repr_str = repr(chain)
        assert 'C' in repr_str
        assert '10' in repr_str
        assert '80' in repr_str


class TestStandardizationResult:
    """Test StandardizationResult dataclass."""

    def test_result_creation(self):
        """Test creating StandardizationResult object."""
        result = StandardizationResult(
            task_name="test_task",
            status="OK",
            num_chains=5,
            chain_mapping={'X': 'C', 'Y': 'B'}
        )
        assert result.task_name == "test_task"
        assert result.status == "OK"
        assert result.num_chains == 5
        assert result.chain_mapping == {'X': 'C', 'Y': 'B'}

    def test_result_to_dict(self):
        """Test converting result to dictionary."""
        result = StandardizationResult(
            task_name="test",
            status="OK",
            num_chains=5,
            chain_mapping={'X': 'C', 'Y': 'B'},
            processing_time=1.5
        )
        result_dict = result.to_dict()

        assert result_dict['TaskName'] == "test"
        assert result_dict['Status'] == "OK"
        assert result_dict['NumChains'] == 5
        assert 'X→C' in result_dict['ChainMapping']
        assert result_dict['ProcessingTime_sec'] == 1.5


class TestPDBChainStandardizer:
    """Test PDBChainStandardizer class."""

    def test_initialization_default(self):
        """Test standardizer initialization with defaults."""
        standardizer = PDBChainStandardizer()
        assert standardizer.standard_order == ['C', 'B', 'D', 'E', 'A']
        assert standardizer.expected_chain_count == 5
        assert standardizer.max_residue_threshold == 1000

    def test_initialization_custom(self):
        """Test standardizer initialization with custom parameters."""
        custom_order = ['A', 'B', 'C']
        standardizer = PDBChainStandardizer(
            standard_order=custom_order,
            expected_chain_count=3,
            max_residue_threshold=500
        )
        assert standardizer.standard_order == custom_order
        assert standardizer.expected_chain_count == 3
        assert standardizer.max_residue_threshold == 500

    def test_analyze_chains_5chain(self, sample_pdb_5chain):
        """Test chain analysis on 5-chain PDB."""
        standardizer = PDBChainStandardizer()
        chains = standardizer.analyze_chains(sample_pdb_5chain)

        assert len(chains) == 5
        # Check sorting by residue count
        assert chains[0].residue_count == 10   # X -> smallest
        assert chains[1].residue_count == 100  # Y
        assert chains[2].residue_count == 150  # Z
        assert chains[3].residue_count == 220  # W
        assert chains[4].residue_count == 280  # V -> largest

        # Check chain IDs
        assert chains[0].chain_id == 'X'
        assert chains[4].chain_id == 'V'

    def test_analyze_chains_file_not_found(self, temp_dir):
        """Test chain analysis with non-existent file."""
        standardizer = PDBChainStandardizer()
        nonexistent = temp_dir / "nonexistent.pdb"

        with pytest.raises(FileNotFoundError):
            standardizer.analyze_chains(nonexistent)

    def test_analyze_chains_empty_pdb(self, temp_dir):
        """Test chain analysis with empty PDB file."""
        empty_pdb = temp_dir / "empty.pdb"
        empty_pdb.write_text("HEADER Test\n")

        standardizer = PDBChainStandardizer()

        with pytest.raises(ValueError, match="No protein chains found"):
            standardizer.analyze_chains(empty_pdb)

    def test_create_mapping_ok(self, sample_pdb_5chain):
        """Test creating chain mapping for valid 5-chain complex."""
        standardizer = PDBChainStandardizer()
        chains = standardizer.analyze_chains(sample_pdb_5chain)
        mapping, status = standardizer.create_mapping(chains)

        assert status == "OK"
        assert mapping['X'] == 'C'  # Smallest -> C (peptide)
        assert mapping['Y'] == 'B'  # Second -> B (HLA-β)
        assert mapping['Z'] == 'D'  # Third -> D (TCR-α)
        assert mapping['W'] == 'E'  # Fourth -> E (TCR-β)
        assert mapping['V'] == 'A'  # Largest -> A (HLA-α)

    def test_create_mapping_multichain(self, sample_pdb_multichain):
        """Test creating mapping with too many chains."""
        standardizer = PDBChainStandardizer()
        chains = standardizer.analyze_chains(sample_pdb_multichain)
        mapping, status = standardizer.create_mapping(chains)

        assert status == "MULTICHAIN"
        assert mapping == {}

    def test_create_mapping_insufficient(self, sample_pdb_insufficient):
        """Test creating mapping with too few chains."""
        standardizer = PDBChainStandardizer()
        chains = standardizer.analyze_chains(sample_pdb_insufficient)
        mapping, status = standardizer.create_mapping(chains)

        assert status == "INSUFFICIENT"
        assert mapping == {}

    def test_is_already_standard_true(self, sample_pdb_standard):
        """Test checking if PDB is already in standard order."""
        standardizer = PDBChainStandardizer()
        chains = standardizer.analyze_chains(sample_pdb_standard)
        mapping, _ = standardizer.create_mapping(chains)

        is_standard = standardizer.is_already_standard(chains, mapping)
        assert is_standard is True

    def test_is_already_standard_false(self, sample_pdb_5chain):
        """Test checking if PDB is not in standard order."""
        standardizer = PDBChainStandardizer()
        chains = standardizer.analyze_chains(sample_pdb_5chain)
        mapping, _ = standardizer.create_mapping(chains)

        is_standard = standardizer.is_already_standard(chains, mapping)
        assert is_standard is False

    def test_standardize_file(self, sample_pdb_5chain, temp_dir):
        """Test standardizing a PDB file."""
        standardizer = PDBChainStandardizer()
        chains = standardizer.analyze_chains(sample_pdb_5chain)
        mapping, _ = standardizer.create_mapping(chains)

        output_pdb = temp_dir / "standardized.pdb"
        standardizer.standardize_file(sample_pdb_5chain, output_pdb, mapping)

        # Verify output file exists
        assert output_pdb.exists()

        # Verify chains were remapped
        output_content = output_pdb.read_text()
        assert 'C' in output_content  # Peptide chain
        assert 'A' in output_content  # HLA-α chain

        # Count chain occurrences
        c_count = output_content.count(' C ')
        assert c_count == 10  # Peptide has 10 residues

    def test_process_single_ok(self, sample_pdb_5chain, temp_dir):
        """Test processing a single PDB file successfully."""
        standardizer = PDBChainStandardizer()
        output_pdb = temp_dir / "output.pdb"

        result = standardizer.process_single(
            input_pdb=sample_pdb_5chain,
            output_pdb=output_pdb,
            task_name="test_task"
        )

        assert result.status == "OK"
        assert result.task_name == "test_task"
        assert result.num_chains == 5
        assert len(result.chain_mapping) == 5
        assert output_pdb.exists()

    def test_process_single_already_standard(self, sample_pdb_standard, temp_dir):
        """Test processing PDB that's already in standard order."""
        standardizer = PDBChainStandardizer()
        output_pdb = temp_dir / "output.pdb"

        result = standardizer.process_single(
            input_pdb=sample_pdb_standard,
            output_pdb=output_pdb,
            task_name="test_standard",
            skip_if_standard=True
        )

        assert result.status == "ALREADY_STANDARD"
        assert result.num_chains == 5

    def test_process_single_no_pdb(self, temp_dir):
        """Test processing non-existent PDB file."""
        standardizer = PDBChainStandardizer()
        nonexistent = temp_dir / "nonexistent.pdb"
        output_pdb = temp_dir / "output.pdb"

        result = standardizer.process_single(
            input_pdb=nonexistent,
            output_pdb=output_pdb,
            task_name="test_missing"
        )

        assert result.status == "NO_PDB"
        assert result.error_message is not None

    def test_process_single_multichain(self, sample_pdb_multichain, temp_dir):
        """Test processing PDB with too many chains."""
        standardizer = PDBChainStandardizer()
        output_pdb = temp_dir / "output.pdb"

        result = standardizer.process_single(
            input_pdb=sample_pdb_multichain,
            output_pdb=output_pdb,
            task_name="test_multi"
        )

        assert result.status == "MULTICHAIN"
        assert result.num_chains == 7
        assert result.error_message is not None

    def test_batch_process(self, sample_pdb_5chain, sample_pdb_standard, temp_dir):
        """Test batch processing multiple PDB files."""
        standardizer = PDBChainStandardizer()

        pairs = [
            (sample_pdb_5chain, temp_dir / "out1.pdb", "task1"),
            (sample_pdb_standard, temp_dir / "out2.pdb", "task2")
        ]

        results = standardizer.batch_process(pairs, n_processes=1)

        assert len(results) == 2
        assert results[0].status == "OK"
        assert results[1].status == "ALREADY_STANDARD"

    def test_save_report(self, temp_dir):
        """Test saving results to CSV report."""
        standardizer = PDBChainStandardizer()

        results = [
            StandardizationResult(
                task_name="task1",
                status="OK",
                num_chains=5,
                chain_mapping={'X': 'C', 'Y': 'B', 'Z': 'D', 'W': 'E', 'V': 'A'},
                processing_time=1.2
            ),
            StandardizationResult(
                task_name="task2",
                status="ALREADY_STANDARD",
                num_chains=5,
                chain_mapping={'C': 'C', 'B': 'B', 'D': 'D', 'E': 'E', 'A': 'A'},
                processing_time=0.8
            )
        ]

        csv_file = temp_dir / "report.csv"
        standardizer.save_report(results, csv_file)

        assert csv_file.exists()

        # Check CSV content
        content = csv_file.read_text()
        assert "TaskName" in content
        assert "task1" in content
        assert "task2" in content
        assert "OK" in content

    def test_validate_standardized_pdb_valid(self, sample_pdb_standard):
        """Test validating a correctly standardized PDB."""
        standardizer = PDBChainStandardizer()
        is_valid, message = standardizer.validate_standardized_pdb(sample_pdb_standard)

        assert is_valid is True
        assert "Valid" in message

    def test_validate_standardized_pdb_invalid(self, sample_pdb_5chain):
        """Test validating a non-standard PDB."""
        standardizer = PDBChainStandardizer()
        is_valid, message = standardizer.validate_standardized_pdb(sample_pdb_5chain)

        assert is_valid is False
        assert "does not match" in message

    def test_validate_standardized_pdb_wrong_count(self, sample_pdb_insufficient):
        """Test validating PDB with wrong chain count."""
        standardizer = PDBChainStandardizer()
        is_valid, message = standardizer.validate_standardized_pdb(sample_pdb_insufficient)

        assert is_valid is False
        assert "Expected 5 chains" in message


class TestIntegration:
    """Integration tests for complete workflows."""

    def test_complete_workflow(self, sample_pdb_5chain, temp_dir):
        """Test complete standardization workflow."""
        standardizer = PDBChainStandardizer()
        output_pdb = temp_dir / "final.pdb"

        # Process
        result = standardizer.process_single(
            input_pdb=sample_pdb_5chain,
            output_pdb=output_pdb,
            task_name="integration_test"
        )

        assert result.status == "OK"

        # Validate output
        is_valid, message = standardizer.validate_standardized_pdb(output_pdb)
        assert is_valid is True

        # Re-process should detect already standard
        result2 = standardizer.process_single(
            input_pdb=output_pdb,
            output_pdb=temp_dir / "final2.pdb",
            task_name="integration_test_2",
            skip_if_standard=True
        )

        assert result2.status == "ALREADY_STANDARD"

    def test_batch_with_report(self, sample_pdb_5chain, sample_pdb_multichain,
                                sample_pdb_insufficient, temp_dir):
        """Test batch processing with report generation."""
        standardizer = PDBChainStandardizer()

        pairs = [
            (sample_pdb_5chain, temp_dir / "out1.pdb", "task1"),
            (sample_pdb_multichain, temp_dir / "out2.pdb", "task2"),
            (sample_pdb_insufficient, temp_dir / "out3.pdb", "task3")
        ]

        results = standardizer.batch_process(pairs, n_processes=1)

        # Save report
        report_file = temp_dir / "batch_report.csv"
        standardizer.save_report(results, report_file)

        assert report_file.exists()

        # Verify results
        ok_results = [r for r in results if r.status == "OK"]
        multi_results = [r for r in results if r.status == "MULTICHAIN"]
        insuf_results = [r for r in results if r.status == "INSUFFICIENT"]

        assert len(ok_results) == 1
        assert len(multi_results) == 1
        assert len(insuf_results) == 1


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
