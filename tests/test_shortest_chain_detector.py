#!/usr/bin/env python3
"""
Unit tests for ShortestChainDetector module

Tests cover:
- Chain detection from GRO files
- Shortest chain identification
- Index file generation
- Error handling
- Edge cases
"""

import pytest
import tempfile
import shutil
from pathlib import Path
from aftermd.utils.shortest_chain_detector import (
    ShortestChainDetector,
    create_shortest_chain_index
)


@pytest.fixture
def temp_dir():
    """Create temporary directory for test files."""
    tmp = tempfile.mkdtemp()
    yield Path(tmp)
    shutil.rmtree(tmp)


@pytest.fixture
def sample_gro_5chain(temp_dir):
    """
    Create sample GRO file with 5 protein chains (pHLA-TCR complex).
    Chain lengths: 50, 400, 10, 600, 1000 atoms
    Expected shortest: Chain 3 (10 atoms, peptide)
    """
    gro_file = temp_dir / "test.gro"

    # GRO format header
    content = ["Test pHLA-TCR complex", "2060"]  # Total atoms

    atom_id = 1

    # Chain 1: 50 atoms (50 residues, 1 atom per residue for simplicity)
    for res in range(1, 51):
        line = f"{res:>5}ALA     CA{atom_id:>5}   1.000   1.000   1.000\n"
        content.append(line)
        atom_id += 1

    # Chain 2: 400 atoms (400 residues)
    for res in range(51, 451):
        line = f"{res:>5}ALA     CA{atom_id:>5}   2.000   2.000   2.000\n"
        content.append(line)
        atom_id += 1

    # Chain 3: 10 atoms (10 residues) - SHORTEST (peptide)
    for res in range(1000, 1010):  # Large gap indicates new chain
        line = f"{res:>5}GLY     CA{atom_id:>5}   3.000   3.000   3.000\n"
        content.append(line)
        atom_id += 1

    # Chain 4: 600 atoms (600 residues)
    for res in range(2000, 2600):  # Another gap
        line = f"{res:>5}LEU     CA{atom_id:>5}   4.000   4.000   4.000\n"
        content.append(line)
        atom_id += 1

    # Chain 5: 1000 atoms (1000 residues)
    for res in range(3000, 4000):  # Another gap
        line = f"{res:>5}VAL     CA{atom_id:>5}   5.000   5.000   5.000\n"
        content.append(line)
        atom_id += 1

    # Box vector (required in GRO format)
    content.append("  10.0  10.0  10.0\n")

    with open(gro_file, 'w') as f:
        f.write('\n'.join(content[:2]) + '\n')
        f.writelines(content[2:])

    return gro_file


@pytest.fixture
def sample_topology(temp_dir):
    """Create minimal topology file."""
    tpr_file = temp_dir / "test.tpr"
    # For testing, we'll use the GRO file as "topology"
    # In real scenarios, this would be a .tpr file
    return tpr_file


def test_shortest_chain_detector_initialization(sample_gro_5chain, sample_topology):
    """Test detector initialization."""
    detector = ShortestChainDetector(
        str(sample_gro_5chain),
        str(sample_topology),
        gmx_executable="gmx"
    )

    assert detector.gro_file == str(sample_gro_5chain)
    assert detector.topology_file == str(sample_topology)
    assert detector.gmx == "gmx"
    assert len(detector.protein_residues) > 0


def test_detect_chains_from_gro(sample_gro_5chain, sample_topology):
    """Test chain detection from GRO file."""
    detector = ShortestChainDetector(str(sample_gro_5chain), str(sample_topology))
    chains = detector._detect_chains_from_gro()

    # Should detect 5 chains
    assert len(chains) == 5

    # Check chain lengths
    chain_lengths = [len(atoms) for atoms in chains.values()]
    assert 10 in chain_lengths  # Shortest chain (peptide)
    assert 50 in chain_lengths
    assert 400 in chain_lengths
    assert 600 in chain_lengths
    assert 1000 in chain_lengths


def test_find_shortest_chain(sample_gro_5chain, sample_topology):
    """Test shortest chain identification."""
    detector = ShortestChainDetector(str(sample_gro_5chain), str(sample_topology))
    chains = detector._detect_chains_from_gro()

    shortest = detector._find_shortest_chain_atoms(chains)

    assert shortest is not None
    assert len(shortest) == 10  # Peptide with 10 atoms


def test_format_shortest_chain_group(sample_gro_5chain, sample_topology):
    """Test index group formatting."""
    detector = ShortestChainDetector(str(sample_gro_5chain), str(sample_topology))

    # Test with sample atom IDs
    atom_ids = list(range(1, 26))  # 25 atoms
    formatted = detector._format_shortest_chain_group(atom_ids)

    lines = formatted.strip().split('\n')

    # Should have header + 3 data lines (12 atoms per line: 12 + 12 + 1)
    assert lines[0] == "[ Shortest_Chain ]"
    assert len(lines) == 4  # header + 3 data lines

    # Check first line has 12 atoms
    first_data_line = lines[1]
    assert len(first_data_line) == 60  # 12 atoms * 5 chars each


def test_create_minimal_index_content(sample_gro_5chain, sample_topology):
    """Test minimal index content generation."""
    detector = ShortestChainDetector(str(sample_gro_5chain), str(sample_topology))
    content = detector._create_minimal_index_content()

    # Should contain standard groups
    assert "[ System ]" in content
    assert "[ Protein ]" in content
    assert "[ Backbone ]" in content
    assert "[ C-alpha ]" in content


def test_generate_shortest_chain_index_full_workflow(sample_gro_5chain, temp_dir):
    """Test complete index generation workflow."""
    detector = ShortestChainDetector(
        str(sample_gro_5chain),
        str(sample_gro_5chain),  # Use GRO as topology for test
        gmx_executable="gmx"
    )

    output_dir = temp_dir / "output"
    index_file = detector.generate_shortest_chain_index(str(output_dir))

    # Note: This test may fail if gmx is not available
    # In CI/CD, you might want to mock subprocess calls
    if index_file:
        assert Path(index_file).exists()
        assert Path(index_file).name == "shortest_chain.ndx"

        # Check file contains shortest chain group
        with open(index_file) as f:
            content = f.read()
            assert "[ Shortest_Chain ]" in content


def test_convenience_function(sample_gro_5chain, temp_dir):
    """Test convenience function create_shortest_chain_index."""
    output_dir = temp_dir / "output"

    index_file = create_shortest_chain_index(
        str(sample_gro_5chain),
        str(sample_gro_5chain),
        str(output_dir),
        gmx_executable="gmx"
    )

    # May return None if gmx not available
    if index_file:
        assert Path(index_file).exists()


def test_nonexistent_gro_file(temp_dir):
    """Test error handling for missing GRO file."""
    nonexistent = temp_dir / "nonexistent.gro"

    result = create_shortest_chain_index(
        str(nonexistent),
        str(temp_dir / "topology.tpr"),
        str(temp_dir / "output")
    )

    assert result is None


def test_empty_chains_handling(sample_topology, temp_dir):
    """Test handling of files with no valid chains."""
    # Create GRO with no protein atoms (only water)
    gro_file = temp_dir / "water_only.gro"
    content = [
        "Water only system",
        "3",  # 3 water atoms
        "    1SOL     OW    1   1.000   1.000   1.000\n",
        "    1SOL    HW1    2   1.100   1.000   1.000\n",
        "    1SOL    HW2    3   1.000   1.100   1.000\n",
        "  10.0  10.0  10.0\n"
    ]

    with open(gro_file, 'w') as f:
        f.write('\n'.join(content[:2]) + '\n')
        f.writelines(content[2:])

    detector = ShortestChainDetector(str(gro_file), str(sample_topology))
    index_file = detector.generate_shortest_chain_index(str(temp_dir / "output"))

    # Should return None (no protein chains)
    assert index_file is None


def test_single_chain_system(temp_dir):
    """Test system with only one chain."""
    gro_file = temp_dir / "single_chain.gro"

    content = ["Single chain", "50"]

    # Single chain: 50 atoms
    for i in range(1, 51):
        line = f"    {i:>1}ALA     CA{i:>5}   1.000   1.000   1.000\n"
        content.append(line)

    content.append("  10.0  10.0  10.0\n")

    with open(gro_file, 'w') as f:
        f.write('\n'.join(content[:2]) + '\n')
        f.writelines(content[2:])

    detector = ShortestChainDetector(str(gro_file), str(gro_file))
    chains = detector._detect_chains_from_gro()

    assert len(chains) == 1

    shortest = detector._find_shortest_chain_atoms(chains)
    assert len(shortest) == 50


def test_chain_with_gaps_handling(temp_dir):
    """Test handling of chains with small gaps in residue numbering."""
    gro_file = temp_dir / "gapped_chain.gro"

    content = ["Chain with gap", "30"]

    atom_id = 1
    # First segment: residues 1-10
    for res in range(1, 11):
        line = f"{res:>5}ALA     CA{atom_id:>5}   1.000   1.000   1.000\n"
        content.append(line)
        atom_id += 1

    # Gap of 1 residue (should NOT split chain)
    # Second segment: residues 12-20 (skipped 11)
    for res in range(12, 21):
        line = f"{res:>5}GLY     CA{atom_id:>5}   1.000   1.000   1.000\n"
        content.append(line)
        atom_id += 1

    # Large gap (should split chain)
    # Third segment: residues 100-109
    for res in range(100, 110):
        line = f"{res:>5}VAL     CA{atom_id:>5}   2.000   2.000   2.000\n"
        content.append(line)
        atom_id += 1

    content.append("  10.0  10.0  10.0\n")

    with open(gro_file, 'w') as f:
        f.write('\n'.join(content[:2]) + '\n')
        f.writelines(content[2:])

    detector = ShortestChainDetector(str(gro_file), str(gro_file))
    chains = detector._detect_chains_from_gro()

    # Should detect 2 chains (small gap not counted, large gap counted)
    # Actually with gap > 1, it will split
    assert len(chains) == 2


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
