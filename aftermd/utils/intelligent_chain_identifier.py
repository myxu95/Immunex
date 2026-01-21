"""
Intelligent Chain Identifier for pHLA-TCR Complexes

Uses ANARCI to identify TCR alpha/beta chains, combined with length-based
heuristics for peptide and beta2-microglobulin identification.

Author: AfterMD Development Team
"""

import logging
from pathlib import Path
from typing import Dict, List, Optional, Tuple
from dataclasses import dataclass

from .pdb_sequence_extractor import PDBSequenceExtractor
from .cdr_manager import ANARCIWrapper

logger = logging.getLogger(__name__)


@dataclass
class ChainIdentification:
    """Result of chain identification."""
    chain_id: str
    length: int
    sequence: str
    chain_type: str  # 'peptide', 'beta2m', 'TCR_alpha', 'TCR_beta', 'HLA_alpha'
    confidence: float  # 0.0 to 1.0
    anarci_result: Optional[Dict] = None


class IntelligentChainIdentifier:
    """
    Intelligent chain identifier using ANARCI and length heuristics.

    Identification strategy:
    1. Identify peptide by length (shortest, typically 5-22 aa)
    2. Identify beta2-microglobulin by length (~99-100 aa)
    3. Use ANARCI to identify remaining 3 chains:
       - 2 TCR chains (alpha and beta)
       - 1 HLA-alpha chain (longest)

    Standard chain assignment:
    - Chain A: HLA-alpha (longest)
    - Chain B: beta2-microglobulin (~100 aa)
    - Chain C: Peptide (shortest)
    - Chain D: TCR-alpha
    - Chain E: TCR-beta
    """

    # Length thresholds for identification
    PEPTIDE_MAX_LENGTH = 25
    BETA2M_MIN_LENGTH = 95
    BETA2M_MAX_LENGTH = 105
    HLA_ALPHA_MIN_LENGTH = 260

    def __init__(self, use_anarci: bool = True):
        """
        Initialize intelligent chain identifier.

        Args:
            use_anarci: Whether to use ANARCI for TCR identification
        """
        self.use_anarci = use_anarci
        self.sequence_extractor = PDBSequenceExtractor()

        if use_anarci:
            try:
                self.anarci_wrapper = ANARCIWrapper(allow_fallback=False)
                logger.info("ANARCI available for TCR chain identification")
            except RuntimeError as e:
                logger.warning(f"ANARCI not available: {e}")
                logger.warning("Will use fallback length-based identification")
                self.use_anarci = False
                self.anarci_wrapper = None
        else:
            self.anarci_wrapper = None

    def identify_chains(
        self,
        pdb_file: str
    ) -> Dict[str, ChainIdentification]:
        """
        Identify all chains in PDB file using intelligent strategy.

        Args:
            pdb_file: Path to PDB file

        Returns:
            Dictionary mapping original chain_id to ChainIdentification

        Example:
            >>> identifier = IntelligentChainIdentifier()
            >>> results = identifier.identify_chains("1ao7.pdb")
            >>> for chain_id, info in results.items():
            ...     print(f"{chain_id}: {info.chain_type} ({info.length} aa)")
        """
        # Step 1: Extract sequences
        chains_data = self.sequence_extractor.extract_sequences_from_pdb(pdb_file)

        if len(chains_data) != 5:
            logger.warning(
                f"Expected 5 chains in pHLA-TCR complex, found {len(chains_data)}"
            )

        # Step 2: Sort by length
        sorted_chains = sorted(
            chains_data.items(),
            key=lambda x: x[1]['length']
        )

        identifications = {}

        # Step 3: Identify peptide (shortest)
        peptide_chain_id, peptide_data = sorted_chains[0]
        if peptide_data['length'] <= self.PEPTIDE_MAX_LENGTH:
            identifications[peptide_chain_id] = ChainIdentification(
                chain_id=peptide_chain_id,
                length=peptide_data['length'],
                sequence=peptide_data['sequence'],
                chain_type='peptide',
                confidence=1.0
            )
            logger.info(f"Identified {peptide_chain_id} as peptide ({peptide_data['length']} aa)")
        else:
            logger.warning(
                f"Shortest chain {peptide_chain_id} has {peptide_data['length']} aa, "
                f"exceeds peptide threshold ({self.PEPTIDE_MAX_LENGTH})"
            )
            identifications[peptide_chain_id] = ChainIdentification(
                chain_id=peptide_chain_id,
                length=peptide_data['length'],
                sequence=peptide_data['sequence'],
                chain_type='unknown',
                confidence=0.5
            )

        # Step 4: Identify beta2-microglobulin (~100 aa)
        beta2m_found = False
        for chain_id, chain_data in sorted_chains[1:]:
            if self.BETA2M_MIN_LENGTH <= chain_data['length'] <= self.BETA2M_MAX_LENGTH:
                identifications[chain_id] = ChainIdentification(
                    chain_id=chain_id,
                    length=chain_data['length'],
                    sequence=chain_data['sequence'],
                    chain_type='beta2m',
                    confidence=0.95
                )
                logger.info(f"Identified {chain_id} as beta2m ({chain_data['length']} aa)")
                beta2m_found = True
                break

        if not beta2m_found:
            logger.warning("Could not identify beta2m chain by length heuristic")

        # Step 5: Identify remaining chains using ANARCI
        remaining_chains = {
            cid: cdata for cid, cdata in chains_data.items()
            if cid not in identifications
        }

        if len(remaining_chains) == 3:
            tcr_hla_identifications = self._identify_tcr_and_hla(remaining_chains)
            identifications.update(tcr_hla_identifications)
        else:
            logger.warning(
                f"Expected 3 remaining chains (2 TCR + 1 HLA-alpha), "
                f"found {len(remaining_chains)}"
            )
            # Fallback: mark as unknown
            for chain_id, chain_data in remaining_chains.items():
                identifications[chain_id] = ChainIdentification(
                    chain_id=chain_id,
                    length=chain_data['length'],
                    sequence=chain_data['sequence'],
                    chain_type='unknown',
                    confidence=0.0
                )

        return identifications

    def _identify_tcr_and_hla(
        self,
        chains: Dict[str, Dict]
    ) -> Dict[str, ChainIdentification]:
        """
        Identify TCR alpha, TCR beta, and HLA-alpha from 3 remaining chains.

        Strategy:
        1. Use ANARCI to identify TCR chains
        2. Longest chain is HLA-alpha
        3. Distinguish TCR-alpha from TCR-beta using ANARCI output

        Args:
            chains: Dictionary of 3 remaining chains

        Returns:
            Dictionary of ChainIdentification results
        """
        identifications = {}

        if not self.use_anarci or self.anarci_wrapper is None:
            # Fallback: length-based identification
            return self._identify_by_length_fallback(chains)

        # Run ANARCI on all 3 chains
        anarci_results = {}
        for chain_id, chain_data in chains.items():
            try:
                result = self.anarci_wrapper.run_anarci(
                    sequence=chain_data['sequence'],
                    chain_type='TCR'
                )
                anarci_results[chain_id] = result
                logger.info(
                    f"ANARCI result for {chain_id}: "
                    f"chain_type={result.get('chain_type', 'unknown')}"
                )
            except Exception as e:
                logger.warning(f"ANARCI failed for chain {chain_id}: {e}")
                anarci_results[chain_id] = None

        # Identify chains based on ANARCI results
        tcr_alpha_candidates = []
        tcr_beta_candidates = []
        hla_alpha_candidates = []

        for chain_id, result in anarci_results.items():
            chain_data = chains[chain_id]

            if result and 'chain_type' in result:
                chain_type_anarci = result['chain_type'].lower()

                if 'alpha' in chain_type_anarci and 'tcr' in chain_type_anarci:
                    tcr_alpha_candidates.append((chain_id, chain_data, result))
                elif 'beta' in chain_type_anarci and 'tcr' in chain_type_anarci:
                    tcr_beta_candidates.append((chain_id, chain_data, result))
                else:
                    # Likely HLA-alpha (ANARCI doesn't recognize it as TCR)
                    hla_alpha_candidates.append((chain_id, chain_data, result))
            else:
                # No TCR recognition - likely HLA-alpha
                hla_alpha_candidates.append((chain_id, chain_data, None))

        # Assign TCR-alpha
        if len(tcr_alpha_candidates) == 1:
            chain_id, chain_data, result = tcr_alpha_candidates[0]
            identifications[chain_id] = ChainIdentification(
                chain_id=chain_id,
                length=chain_data['length'],
                sequence=chain_data['sequence'],
                chain_type='TCR_alpha',
                confidence=0.95,
                anarci_result=result
            )
            logger.info(f"Identified {chain_id} as TCR-alpha (ANARCI)")
        elif len(tcr_alpha_candidates) > 1:
            logger.warning(f"Multiple TCR-alpha candidates: {[c[0] for c in tcr_alpha_candidates]}")
            # Use shortest as alpha
            tcr_alpha_candidates.sort(key=lambda x: x[1]['length'])
            chain_id, chain_data, result = tcr_alpha_candidates[0]
            identifications[chain_id] = ChainIdentification(
                chain_id=chain_id,
                length=chain_data['length'],
                sequence=chain_data['sequence'],
                chain_type='TCR_alpha',
                confidence=0.7,
                anarci_result=result
            )

        # Assign TCR-beta
        if len(tcr_beta_candidates) == 1:
            chain_id, chain_data, result = tcr_beta_candidates[0]
            identifications[chain_id] = ChainIdentification(
                chain_id=chain_id,
                length=chain_data['length'],
                sequence=chain_data['sequence'],
                chain_type='TCR_beta',
                confidence=0.95,
                anarci_result=result
            )
            logger.info(f"Identified {chain_id} as TCR-beta (ANARCI)")
        elif len(tcr_beta_candidates) > 1:
            logger.warning(f"Multiple TCR-beta candidates: {[c[0] for c in tcr_beta_candidates]}")
            # Use longest as beta
            tcr_beta_candidates.sort(key=lambda x: x[1]['length'], reverse=True)
            chain_id, chain_data, result = tcr_beta_candidates[0]
            identifications[chain_id] = ChainIdentification(
                chain_id=chain_id,
                length=chain_data['length'],
                sequence=chain_data['sequence'],
                chain_type='TCR_beta',
                confidence=0.7,
                anarci_result=result
            )

        # Assign HLA-alpha (remaining chain, should be longest)
        remaining = [
            (cid, cdata) for cid, cdata in chains.items()
            if cid not in identifications
        ]

        if len(remaining) == 1:
            chain_id, chain_data = remaining[0]
            identifications[chain_id] = ChainIdentification(
                chain_id=chain_id,
                length=chain_data['length'],
                sequence=chain_data['sequence'],
                chain_type='HLA_alpha',
                confidence=0.9
            )
            logger.info(f"Identified {chain_id} as HLA-alpha (longest remaining)")
        elif len(remaining) > 1:
            # Multiple remaining - assign longest as HLA-alpha
            remaining.sort(key=lambda x: x[1]['length'], reverse=True)
            chain_id, chain_data = remaining[0]
            identifications[chain_id] = ChainIdentification(
                chain_id=chain_id,
                length=chain_data['length'],
                sequence=chain_data['sequence'],
                chain_type='HLA_alpha',
                confidence=0.7
            )
            logger.warning(f"Multiple unidentified chains, assigned longest {chain_id} as HLA-alpha")

            # Mark others as unknown
            for chain_id, chain_data in remaining[1:]:
                identifications[chain_id] = ChainIdentification(
                    chain_id=chain_id,
                    length=chain_data['length'],
                    sequence=chain_data['sequence'],
                    chain_type='unknown',
                    confidence=0.3
                )

        return identifications

    def _identify_by_length_fallback(
        self,
        chains: Dict[str, Dict]
    ) -> Dict[str, ChainIdentification]:
        """
        Fallback identification based on length alone (without ANARCI).

        Assumes:
        - Longest chain: HLA-alpha
        - Middle 2 chains: TCR (shorter=alpha, longer=beta)

        Args:
            chains: Dictionary of 3 chains

        Returns:
            Dictionary of ChainIdentification results
        """
        identifications = {}
        sorted_chains = sorted(chains.items(), key=lambda x: x[1]['length'])

        if len(sorted_chains) != 3:
            logger.error(f"Expected 3 chains, got {len(sorted_chains)}")
            return identifications

        # Shortest of the 3: TCR-alpha (usually)
        chain_id, chain_data = sorted_chains[0]
        identifications[chain_id] = ChainIdentification(
            chain_id=chain_id,
            length=chain_data['length'],
            sequence=chain_data['sequence'],
            chain_type='TCR_alpha',
            confidence=0.6  # Lower confidence without ANARCI
        )
        logger.info(f"Identified {chain_id} as TCR-alpha (length fallback)")

        # Middle: TCR-beta (usually)
        chain_id, chain_data = sorted_chains[1]
        identifications[chain_id] = ChainIdentification(
            chain_id=chain_id,
            length=chain_data['length'],
            sequence=chain_data['sequence'],
            chain_type='TCR_beta',
            confidence=0.6
        )
        logger.info(f"Identified {chain_id} as TCR-beta (length fallback)")

        # Longest: HLA-alpha
        chain_id, chain_data = sorted_chains[2]
        identifications[chain_id] = ChainIdentification(
            chain_id=chain_id,
            length=chain_data['length'],
            sequence=chain_data['sequence'],
            chain_type='HLA_alpha',
            confidence=0.8
        )
        logger.info(f"Identified {chain_id} as HLA-alpha (length fallback)")

        return identifications

    def create_standardization_mapping(
        self,
        identifications: Dict[str, ChainIdentification]
    ) -> Dict[str, str]:
        """
        Create chain ID mapping for standardization.

        Standard ordering:
        - Chain A: HLA-alpha
        - Chain B: beta2m
        - Chain C: Peptide
        - Chain D: TCR-alpha
        - Chain E: TCR-beta

        Args:
            identifications: Dictionary from identify_chains()

        Returns:
            Dictionary mapping old_chain_id -> new_chain_id

        Example:
            >>> mapping = identifier.create_standardization_mapping(identifications)
            >>> print(mapping)  # {'X': 'C', 'Y': 'B', 'Z': 'D', 'W': 'E', 'V': 'A'}
        """
        mapping = {}

        for chain_id, info in identifications.items():
            if info.chain_type == 'peptide':
                mapping[chain_id] = 'C'
            elif info.chain_type == 'beta2m':
                mapping[chain_id] = 'B'
            elif info.chain_type == 'TCR_alpha':
                mapping[chain_id] = 'D'
            elif info.chain_type == 'TCR_beta':
                mapping[chain_id] = 'E'
            elif info.chain_type == 'HLA_alpha':
                mapping[chain_id] = 'A'
            else:
                logger.warning(f"Unknown chain type for {chain_id}: {info.chain_type}")
                mapping[chain_id] = chain_id  # Keep original

        return mapping
