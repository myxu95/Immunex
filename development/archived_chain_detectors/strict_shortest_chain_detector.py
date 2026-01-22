#!/usr/bin/env python3
"""
Strict Shortest Chain Detector

For complex systems, ONLY the shortest chain (peptide) can be used as center group.
No fallback strategies allowed - if shortest chain detection fails, processing must stop.
"""

import re
import logging
import subprocess
from pathlib import Path
from typing import Dict, List, Optional, Tuple
import tempfile

logger = logging.getLogger(__name__)


class StrictShortestChainDetector:
    """
    Strict detector that ONLY identifies shortest chain for centering.
    No fallback strategies - either find shortest chain or fail safely.
    """

    def __init__(self, topology_file: str, gro_file: Optional[str] = None, gmx_executable: str = "gmx"):
        self.topology_file = topology_file
        self.gro_file = gro_file
        self.gmx = gmx_executable

        # Extended protein residue types
        self.protein_residues = {
            # Standard amino acids
            'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE',
            'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL',
            # Modified residues
            'ACE', 'NME', 'HIE', 'HID', 'HIP', 'CYX', 'ASH', 'GLH', 'LYN',
            # Alternative names and caps
            'HSD', 'HSE', 'HSP', 'CSS', 'CYM', 'NTER', 'CTER'
        }

    def detect_shortest_chain_or_fail(self) -> Tuple[str, str]:
        """
        Simplified shortest chain detection using gmx make_ndx + splitch.

        Process:
        1. Use 'gmx make_ndx' to load standard groups
        2. Execute 'splitch 1' to split Protein group by chains
        3. Find the shortest newly created chain group (peptide)
        4. Return the group number for use in trjconv -center

        Returns:
            Tuple of (group_id, index_file_path)

        Raises:
            RuntimeError: If shortest chain cannot be detected
        """
        logger.info("开始简化的严格最短链检测 - 使用gmx make_ndx + splitch")

        try:
            # Execute simplified detection using gmx make_ndx + splitch
            result = self._detect_shortest_chain_splitch()

            if result:
                shortest_group_id, index_file = result
                logger.info(f"✅ 检测到最短链(peptide)")
                logger.info(f"   组序号: {shortest_group_id}")
                logger.info(f"   Index文件: {index_file}")
                return shortest_group_id, index_file
            else:
                raise RuntimeError("pHLA-TCR体系splitch检测失败，无法识别peptide链")

        except Exception as e:
            error_msg = (
                f"❌ pHLA-TCR体系最短链检测失败: {e}\n"
                "   pHLA-TCR复合物应包含5条链: HLA重链+β2m+peptide+TCR(αβ)\n"
                "   请检查:\n"
                "   1. TPR文件是否完整包含所有链信息\n"
                "   2. GROMACS版本是否支持splitch命令\n"
                "   3. 链之间是否有明确的分界标识\n"
                "   4. 可能需要手动验证gmx make_ndx输出"
            )
            logger.error(error_msg)
            raise RuntimeError("pHLA-TCR体系peptide检测失败，停止处理以保证科学正确性")

    def _detect_shortest_chain_splitch(self) -> Optional[Tuple[str, str]]:
        """
        使用gmx make_ndx + splitch 1检测最短链

        Returns:
            Tuple of (group_id, index_file_path) or None if failed
        """
        try:
            # Create temporary index file
            temp_dir = tempfile.mkdtemp(prefix="aftermd_splitch_")
            index_file = Path(temp_dir) / "chains.ndx"

            logger.info(f"创建临时index文件: {index_file}")

            # Use confirmed splitch syntax for the server's GROMACS version
            commands = "splitch 1\nq\n"

            # Execute gmx make_ndx
            logger.debug(f"执行命令: {self.gmx} make_ndx -f {self.topology_file} -o {index_file}")
            logger.debug(f"输入命令: {repr(commands)}")

            result = subprocess.run(
                [self.gmx, "make_ndx", "-f", self.topology_file, "-o", str(index_file)],
                input=commands,
                text=True,
                capture_output=True,
                timeout=60
            )

            logger.info(f"命令返回码: {result.returncode}")

            # Log the complete output for debugging pHLA-TCR systems
            logger.info("=== COMPLETE GMX MAKE_NDX OUTPUT ===")
            logger.info("STDERR:")
            logger.info(result.stderr)
            logger.info("STDOUT:")
            logger.info(result.stdout)
            logger.info("=== END OUTPUT ===")

            if result.returncode != 0:
                logger.warning(f"gmx make_ndx returned non-zero code {result.returncode}, but continuing to parse output")

            # Try to parse the output - check both stderr and stdout
            # In GROMACS 2023.2, chain info appears in stdout after splitch command

            # Parse stderr first (traditional location)
            chain_groups = self._parse_splitch_output(result.stderr)

            # If not found in stderr, try stdout (GROMACS 2023.2 behavior)
            if not chain_groups:
                logger.info("未在stderr中找到链组，尝试解析stdout")
                chain_groups = self._parse_splitch_output(result.stdout)

            if not chain_groups:
                logger.error("未找到splitch生成的链组")
                return None

            # Find shortest chain group
            shortest_group = self._find_shortest_chain_group(chain_groups)

            if not shortest_group:
                logger.error("无法确定最短链组")
                return None

            group_id, atom_count = shortest_group
            logger.info(f"检测到最短链: 组{group_id}, {atom_count}个原子")

            # Validate the index file exists and is readable
            if not index_file.exists():
                logger.error("Index文件未成功创建")
                return None

            return str(group_id), str(index_file)

        except Exception as e:
            logger.error(f"splitch检测失败: {e}")
            return None

    def _detect_from_gro_multiformat(self) -> Optional[Tuple[Dict[int, List[int]], List[int]]]:
        """多格式GRO文件解析"""
        if not self.gro_file or not Path(self.gro_file).exists():
            logger.warning("GRO文件不存在，跳过GRO检测")
            return None

        try:
            # Method 1: Standard format parsing
            result = self._parse_gro_standard_format()
            if result and self._is_valid_chain_detection(result):
                return result

            # Method 2: Flexible format parsing
            result = self._parse_gro_flexible_format()
            if result and self._is_valid_chain_detection(result):
                return result

            # Method 3: Residue-gap based parsing
            result = self._parse_gro_residue_gaps()
            if result and self._is_valid_chain_detection(result):
                return result

            return None

        except Exception as e:
            logger.error(f"GRO多格式解析失败: {e}")
            return None

    def _detect_from_combined_analysis(self) -> Optional[Tuple[Dict[int, List[int]], List[int]]]:
        """结合TPR和GRO的分析"""
        try:
            # Combine information from both sources
            # This would implement a more sophisticated approach
            # that cross-validates findings from TPR and GRO
            logger.info("组合分析方法尚未实现")
            return None

        except Exception as e:
            logger.error(f"组合分析失败: {e}")
            return None

    def _parse_gro_standard_format(self) -> Optional[Tuple[Dict[int, List[int]], List[int]]]:
        """标准GRO格式解析 - 更严格的格式检查"""
        try:
            with open(self.gro_file, 'r') as f:
                lines = f.readlines()

            if len(lines) < 3:
                logger.warning("GRO文件格式不正确")
                return None

            # Parse header
            title = lines[0].strip()
            try:
                atom_count = int(lines[1].strip())
            except ValueError:
                logger.warning("无法解析原子数量")
                return None

            # Ensure we have the right number of lines
            expected_lines = atom_count + 3  # title + count + atoms + box
            if len(lines) < expected_lines:
                logger.warning(f"GRO文件行数不足: 期望{expected_lines}, 实际{len(lines)}")
                return None

            atom_lines = lines[2:2+atom_count]
            chains = {}
            current_chain_id = 0
            current_atoms = []
            last_resid = None

            for line_num, line in enumerate(atom_lines, 1):
                try:
                    atom_info = self._parse_gro_line_strict(line, line_num)
                    if not atom_info:
                        continue

                    resid, resname, atom_id = atom_info

                    if resname in self.protein_residues:
                        # Chain boundary detection
                        if last_resid is not None and resid != last_resid + 1:
                            # New chain detected
                            if len(current_atoms) >= 5:  # Minimum peptide size
                                chains[current_chain_id] = current_atoms
                                current_chain_id += 1
                            current_atoms = []

                        current_atoms.append(atom_id)
                        last_resid = resid

                except Exception as e:
                    logger.debug(f"解析第{line_num}行失败: {e}")
                    continue

            # Add final chain
            if len(current_atoms) >= 5:
                chains[current_chain_id] = current_atoms

            if len(chains) < 2:
                logger.warning(f"检测到的链数太少: {len(chains)}")
                return None

            # Find shortest chain
            shortest_chain_atoms = min(chains.values(), key=len)

            logger.info(f"标准格式解析: 检测到{len(chains)}条链")
            for i, atoms in chains.items():
                logger.info(f"  链{i}: {len(atoms)}个原子")

            return chains, shortest_chain_atoms

        except Exception as e:
            logger.error(f"标准GRO解析失败: {e}")
            return None

    def _parse_gro_line_strict(self, line: str, line_num: int) -> Optional[Tuple[int, str, int]]:
        """严格的GRO行解析，支持多种格式变体"""

        if len(line.strip()) < 30:
            return None

        try:
            # GRO format: positions are fixed
            # %5d%-5s%5s%5d%8.3f%8.3f%8.3f
            # resid resname atomname atomid x y z

            # Method 1: Fixed position parsing (most common)
            if len(line) >= 44:
                resid_resname = line[0:5].strip()
                atom_name = line[5:10].strip()
                atom_id_str = line[10:15].strip()

                # Parse residue ID and name
                patterns = [
                    r'^(\d+)([A-Z]{3,4})$',      # 1ALA, 123PRO
                    r'^(\d+)\s+([A-Z]{3,4})$',   # "1 ALA"
                    r'^(\d+)([A-Z]+)$'           # 1ACE, 1NTER
                ]

                for pattern in patterns:
                    match = re.match(pattern, resid_resname)
                    if match:
                        resid = int(match.group(1))
                        resname = match.group(2)
                        atom_id = int(atom_id_str)
                        return resid, resname, atom_id

            # Method 2: Space-separated parsing (fallback)
            parts = line.split()
            if len(parts) >= 6:
                # Try to parse first field as "123ALA" format
                first_field = parts[0]
                match = re.match(r'^(\d+)([A-Z]+)', first_field)
                if match:
                    resid = int(match.group(1))
                    resname = match.group(2)
                    # Atom ID usually in 3rd position for space-separated
                    atom_id = int(parts[2])
                    return resid, resname, atom_id

            return None

        except (ValueError, IndexError) as e:
            logger.debug(f"行{line_num}解析失败: {e}")
            return None

    def _parse_gro_flexible_format(self) -> Optional[Tuple[Dict[int, List[int]], List[int]]]:
        """灵活格式GRO解析 - 处理格式变体"""
        # Implementation similar to standard but with more flexible parsing
        logger.info("灵活格式解析尚未完全实现")
        return None

    def _parse_gro_residue_gaps(self) -> Optional[Tuple[Dict[int, List[int]], List[int]]]:
        """基于残基间隔的链检测"""
        # Implementation that looks for larger gaps in residue numbering
        logger.info("残基间隔检测尚未完全实现")
        return None

    def _parse_splitch_output(self, stderr_output: str) -> List[Tuple[int, int]]:
        """
        解析gmx make_ndx splitch命令的输出，提取链组信息

        Args:
            stderr_output: gmx make_ndx的stderr输出

        Returns:
            List of (group_id, atom_count) tuples for newly created chain groups
        """
        chain_groups = []

        try:
            lines = stderr_output.split('\n')

            # Focus on parsing - detailed output already logged above

            # Multiple patterns to try, ordered by specificity
            patterns = [
                # Pattern 1: Actual GROMACS 2023.2 splitch output format
                r'^(\d+):\s*(\d+)\s+atoms\s*\(\d+\s+to\s+\d+\)',
                # Pattern 2: Standard splitch output with "Protein_chain_"
                r'^\s*(\d+)\s+Protein_chain_\w+\s*:\s*(\d+)\s+atoms',
                # Pattern 3: General chain pattern (more flexible)
                r'^\s*(\d+)\s+\w*[Cc]hain\w*\s*:\s*(\d+)\s+atoms',
                # Pattern 4: Very flexible pattern for any group >= 18
                r'^\s*(\d+)\s+(\S+)\s*:\s*(\d+)\s+atoms'
            ]

            # Try each pattern
            for pattern_num, pattern in enumerate(patterns, 1):
                pattern_chains = []

                for line in lines:
                    match = re.search(pattern, line)
                    if match:
                        if pattern_num == 1:  # Pattern 1: GROMACS 2023.2 format "1: 4343 atoms"
                            chain_id = int(match.group(1))
                            atom_count = int(match.group(2))
                            # Convert chain ID to group ID (add 17 to get group number)
                            group_id = chain_id + 17
                            pattern_chains.append((group_id, atom_count))
                            logger.debug(f"Pattern {pattern_num} 匹配: 链{chain_id}->组{group_id}, {atom_count}个原子")

                        elif pattern_num <= 3:  # Patterns 2-3: direct group_id, atom_count
                            group_id = int(match.group(1))
                            atom_count = int(match.group(2))
                            # Only accept newly created groups (>= 18)
                            if group_id >= 18:
                                pattern_chains.append((group_id, atom_count))
                                logger.debug(f"Pattern {pattern_num} 匹配: 组{group_id}, {atom_count}个原子")

                        else:  # Pattern 4: group_id, group_name, atom_count
                            group_id = int(match.group(1))
                            group_name = match.group(2)
                            atom_count = int(match.group(3))

                            # For pattern 4, only accept groups >= 18 with chain-like names
                            if group_id < 18:
                                continue
                            if not any(keyword in group_name.lower() for keyword in ['chain', 'prot']):
                                continue

                            pattern_chains.append((group_id, atom_count))
                            logger.debug(f"Pattern {pattern_num} 匹配: 组{group_id}, {atom_count}个原子")

                if pattern_chains:
                    chain_groups = pattern_chains
                    logger.info(f"使用Pattern {pattern_num}解析到 {len(chain_groups)} 个链组")
                    break

            # If still no chains found, show diagnostic info
            if not chain_groups:
                logger.warning("未找到任何链组，显示诊断信息:")

                # Look for "Splitting" message
                splitting_found = False
                for line in lines:
                    if 'splitting' in line.lower() or 'splitch' in line.lower():
                        logger.warning(f"找到分割消息: {line.strip()}")
                        splitting_found = True

                if not splitting_found:
                    logger.warning("未找到splitting消息，splitch可能未执行")

                # Show any lines with numbers that might be groups
                logger.warning("可能的组行:")
                for line in lines:
                    if re.search(r'^\s*\d+\s+\S+.*atoms', line):
                        logger.warning(f"  {line.strip()}")

            logger.info(f"最终解析到 {len(chain_groups)} 个链组")
            return chain_groups

        except Exception as e:
            logger.error(f"解析splitch输出失败: {e}")
            return []

    def _find_shortest_chain_group(self, chain_groups: List[Tuple[int, int]]) -> Optional[Tuple[int, int]]:
        """
        从链组列表中找到最短的链（peptide）

        Args:
            chain_groups: List of (group_id, atom_count) tuples

        Returns:
            (group_id, atom_count) of shortest chain, or None if invalid
        """
        if not chain_groups:
            logger.error("没有链组可供选择")
            return None

        if len(chain_groups) < 2:
            logger.error(f"检测到的链数不足: {len(chain_groups)} < 2，不是复合物系统")
            return None

        # Sort by atom count to find shortest
        sorted_chains = sorted(chain_groups, key=lambda x: x[1])
        shortest_group_id, shortest_atom_count = sorted_chains[0]

        # Validate shortest chain is reasonable for peptide (≤20 AAs = ~300 atoms)
        if shortest_atom_count < 10:
            logger.error(f"最短链太短: {shortest_atom_count} < 10原子，可能不是有效peptide")
            return None

        if shortest_atom_count > 300:
            logger.error(f"最短链太长: {shortest_atom_count} > 300原子，超出peptide范围(≤20AAs)")
            return None

        # Check size distribution - ensure meaningful difference
        longest_atom_count = sorted_chains[-1][1]
        if shortest_atom_count * 2 > longest_atom_count:
            logger.warning(f"链尺寸差异不明显: 最短{shortest_atom_count}, 最长{longest_atom_count}")

        logger.info(f"选择最短链作为peptide: 组{shortest_group_id}, {shortest_atom_count}个原子")
        for group_id, atom_count in sorted_chains:
            logger.info(f"  链组{group_id}: {atom_count}个原子")

        return shortest_group_id, shortest_atom_count

    def _is_valid_chain_detection(self, result: Tuple[Dict[int, List[int]], List[int]]) -> bool:
        """验证链检测结果是否有效"""
        chains, shortest_chain_atoms = result

        # Must have at least 2 chains (complex system)
        if len(chains) < 2:
            logger.warning(f"检测到的链数不足: {len(chains)} < 2")
            return False

        # Shortest chain must be reasonable peptide size
        if len(shortest_chain_atoms) < 5:
            logger.warning(f"最短链太短: {len(shortest_chain_atoms)} < 5")
            return False

        # Shortest chain shouldn't be too long (likely not a peptide)
        if len(shortest_chain_atoms) > 1000:
            logger.warning(f"最短链太长: {len(shortest_chain_atoms)} > 1000")
            return False

        # Check chain size distribution
        chain_sizes = [len(atoms) for atoms in chains.values()]
        min_size = min(chain_sizes)
        max_size = max(chain_sizes)

        # Ensure there's meaningful size difference
        if min_size * 3 > max_size:
            logger.warning(f"链尺寸差异不明显: 最小{min_size}, 最大{max_size}")
            return False

        logger.info(f"✅ 链检测结果验证通过")
        logger.info(f"   链数: {len(chains)}")
        logger.info(f"   最短链: {min_size} 原子")
        logger.info(f"   最长链: {max_size} 原子")

        return True

    def _validate_shortest_chain_result(self, result: Tuple[Dict[int, List[int]], List[int]]) -> bool:
        """验证最短链检测结果"""
        return self._is_valid_chain_detection(result)


    def _validate_index_file(self, index_file: str, group_id: str) -> bool:
        """验证index文件是否包含指定的组"""
        try:
            with open(index_file, 'r') as f:
                content = f.read()

            # Check if the group ID exists in the file
            # Look for group headers like "[ Protein_chain_A ]" or similar
            group_pattern = rf'\[\s*\w*\s*\]'
            groups_found = re.findall(group_pattern, content)

            if len(groups_found) >= int(group_id):
                logger.info(f"Index文件验证通过: 包含组{group_id}")
                return True
            else:
                logger.error(f"Index文件验证失败: 未找到组{group_id}")
                return False

        except Exception as e:
            logger.error(f"Index文件验证失败: {e}")
            return False


def detect_shortest_chain_strict(topology_file: str,
                               gro_file: Optional[str] = None,
                               gmx_executable: str = "gmx") -> Tuple[str, str]:
    """
    简化的严格最短链检测 - 使用gmx make_ndx + splitch。

    Process:
    1. 执行 'gmx make_ndx -f topology.tpr' 加载标准组
    2. 执行 'splitch 1' 命令按链拆分Protein组
    3. 解析输出找到新生成的链组中最短的一个(peptide)
    4. 返回该组的序号，用于后续 'gmx trjconv -center' 步骤

    Args:
        topology_file: TPR文件路径
        gro_file: GRO文件路径 (当前未使用，保留接口兼容性)
        gmx_executable: GROMACS可执行文件路径

    Returns:
        Tuple of (group_id, index_file_path)
        - group_id: 最短链组的序号 (如 "18", "19" 等)
        - index_file_path: 生成的index文件路径

    Raises:
        RuntimeError: 如果无法检测到最短链

    Example:
        group_id, index_file = detect_shortest_chain_strict("system.tpr")
        # 然后可以在trjconv中使用: gmx trjconv -center -pbc mol -n index_file
        # 选择组group_id进行居中
    """
    detector = StrictShortestChainDetector(topology_file, gro_file, gmx_executable)
    return detector.detect_shortest_chain_or_fail()