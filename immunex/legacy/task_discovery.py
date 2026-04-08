"""
Task Discovery - Automatic task detection from directory structures.

This module provides tools to automatically discover MD simulation tasks
from directory structures and create PipelineContext instances.
"""

from pathlib import Path
import warnings
from typing import List, Dict, Optional, Callable
import logging

from ..core.context import PipelineContext

logger = logging.getLogger(__name__)


class TaskDiscovery:
    """
    Legacy task discovery tool returning PipelineContext objects.

    This tool follows the principle: "Convention over Configuration".
    It uses sensible defaults for common directory structures but allows
    custom discovery rules.

    Example:
        >>> discovery = TaskDiscovery()
        >>> tasks = discovery.discover("/data/simulations")
        >>> print(f"Found {len(tasks)} tasks")
        >>> for task in tasks:
        ...     print(f"  - {task.system_id}")
    """

    def __init__(self):
        """Initialize legacy task discovery."""
        warnings.warn(
            "TaskDiscovery is legacy. Prefer immunex.core.discover_tasks() for the public discovery contract, or LegacyTaskDiscovery when PipelineContext output is explicitly required.",
            DeprecationWarning,
            stacklevel=2,
        )
        self.discovery_rules = []

    def add_rule(self, rule: Callable[[Path], Optional[Dict]]):
        """
        Add custom discovery rule.

        A discovery rule is a function that takes a Path and returns
        a dict with task information, or None if not a valid task.

        Args:
            rule: Discovery rule function

        Example:
            >>> def my_rule(task_dir):
            ...     if (task_dir / "custom.xtc").exists():
            ...         return {
            ...             'topology': str(task_dir / "custom.tpr"),
            ...             'trajectory_raw': str(task_dir / "custom.xtc")
            ...         }
            ...     return None
            >>>
            >>> discovery.add_rule(my_rule)
        """
        self.discovery_rules.append(rule)

    def discover(self,
                base_directory: str,
                pattern: str = "*",
                structure_type: str = "auto",
                require_files: Optional[List[str]] = None,
                output_base_dir: Optional[str] = None) -> List[PipelineContext]:
        """
        Discover tasks in a base directory.

        Supports multiple directory structures:
        - "auto": Automatically detect structure type
        - "nested": {pdb_id}/{variant}/ structure
        - "flat": {task_name}/ flat structure

        Args:
            base_directory: Base directory to search
            pattern: Subdirectory pattern (glob style)
            structure_type: Directory structure type ("auto", "nested", "flat")
            require_files: Additional required files
            output_base_dir: Base directory for outputs (default: ./results)

        Returns:
            List of PipelineContext instances

        Example:
            >>> # Auto-detect structure
            >>> tasks = discovery.discover("/data/simulations")

            >>> # Nested structure: {pdb_id}/{variant}/
            >>> tasks = discovery.discover("/data/simulations", structure_type="nested")

            >>> # Flat structure: {task_name}/
            >>> tasks = discovery.discover("/data/simulations", structure_type="flat")

            >>> # Custom output directory
            >>> tasks = discovery.discover("/data/sims", output_base_dir="/results/batch1")
        """
        base_path = Path(base_directory)

        if not base_path.exists():
            raise ValueError(f"Base directory not found: {base_directory}")

        if not base_path.is_dir():
            raise ValueError(f"Base directory is not a directory: {base_directory}")

        tasks = []

        if structure_type == "auto":
            structure_type = self._detect_structure_type(base_path)
            logger.info(f"Auto-detected structure type: {structure_type}")

        if structure_type == "nested":
            tasks = self._discover_nested(base_path, pattern, require_files, output_base_dir)
        else:  # flat
            tasks = self._discover_flat(base_path, pattern, require_files, output_base_dir)

        logger.info(f"Discovered {len(tasks)} tasks in {base_directory}")
        return tasks

    def _detect_structure_type(self, base_path: Path) -> str:
        """
        Auto-detect directory structure type.

        Checks if subdirectories contain MD files directly (flat)
        or have further nested directories (nested).

        Args:
            base_path: Base directory path

        Returns:
            "nested" or "flat"
        """
        for subdir in list(base_path.iterdir())[:5]:  # Check first 5 dirs
            if not subdir.is_dir():
                continue

            # Check if MD files exist directly
            if (subdir / "md.xtc").exists() or (subdir / "prod" / "md.xtc").exists():
                return "flat"

            # Check if has nested directories with MD files
            for nested_dir in subdir.iterdir():
                if nested_dir.is_dir():
                    if (nested_dir / "md.xtc").exists():
                        return "nested"

        return "flat"  # Default to flat

    def _discover_nested(self,
                        base_path: Path,
                        pattern: str,
                        require_files: Optional[List[str]],
                        output_base_dir: Optional[str]) -> List[PipelineContext]:
        """
        Discover tasks in nested structure: {pdb_id}/{variant}/

        Args:
            base_path: Base directory
            pattern: Directory pattern
            require_files: Required files
            output_base_dir: Output base directory

        Returns:
            List of PipelineContext instances
        """
        tasks = []

        for pdb_dir in sorted(base_path.glob(pattern)):
            if not pdb_dir.is_dir():
                continue

            pdb_id = pdb_dir.name

            # Scan variant directories
            for variant_dir in sorted(pdb_dir.iterdir()):
                if not variant_dir.is_dir():
                    continue

                variant_name = variant_dir.name
                system_id = f"{pdb_id}_{variant_name}"

                logger.debug(f"Checking nested task: {system_id}")

                # Try custom rules first
                task_data = None
                if self.discovery_rules:
                    for rule in self.discovery_rules:
                        task_data = rule(variant_dir)
                        if task_data:
                            logger.debug(f"Custom rule matched for {system_id}")
                            break

                # Fall back to default discovery
                if not task_data:
                    task_data = self._default_discovery(variant_dir, require_files)

                # Create context if valid task found
                if task_data:
                    task_data['pdb_id'] = pdb_id
                    task_data['variant'] = variant_name
                    context = self._create_context(system_id, task_data, output_base_dir)
                    tasks.append(context)
                    logger.info(f"Discovered nested task: {system_id}")

        return tasks

    def _discover_flat(self,
                      base_path: Path,
                      pattern: str,
                      require_files: Optional[List[str]],
                      output_base_dir: Optional[str]) -> List[PipelineContext]:
        """
        Discover tasks in flat structure: {task_name}/

        Args:
            base_path: Base directory
            pattern: Directory pattern
            require_files: Required files
            output_base_dir: Output base directory

        Returns:
            List of PipelineContext instances
        """
        tasks = []

        for task_dir in sorted(base_path.glob(pattern)):
            if not task_dir.is_dir():
                continue

            task_name = task_dir.name
            logger.debug(f"Checking flat task: {task_name}")

            # Try custom rules first
            task_data = None
            if self.discovery_rules:
                for rule in self.discovery_rules:
                    task_data = rule(task_dir)
                    if task_data:
                        logger.debug(f"Custom rule matched for {task_name}")
                        break

            # Fall back to default discovery
            if not task_data:
                task_data = self._default_discovery(task_dir, require_files)

            # Create context if valid task found
            if task_data:
                context = self._create_context(task_name, task_data, output_base_dir)
                tasks.append(context)
                logger.info(f"Discovered flat task: {task_name}")

        return tasks

    def _default_discovery(self,
                          task_dir: Path,
                          require_files: Optional[List[str]] = None) -> Optional[Dict]:
        """
        Default discovery rule: find md.xtc and md.tpr.

        Searches in two locations:
        1. Directly in task_dir
        2. In task_dir/prod subdirectory

        Args:
            task_dir: Task directory path
            require_files: Additional required files

        Returns:
            Dict with task data, or None if not found
        """
        search_locations = [
            task_dir,
            task_dir / "prod"
        ]

        for location in search_locations:
            if not location.exists():
                continue

            trajectory = location / "md.xtc"
            topology = location / "md.tpr"

            if not (trajectory.exists() and topology.exists()):
                continue

            # Check additional required files
            if require_files:
                missing = [f for f in require_files if not (location / f).exists()]
                if missing:
                    logger.debug(f"Missing required files in {location}: {missing}")
                    continue

            # Try to find PDB structure
            structure_pdb = None
            for pdb_name in ["structure.pdb", "protein.pdb", "complex.pdb"]:
                pdb_path = location / pdb_name
                if pdb_path.exists():
                    structure_pdb = str(pdb_path)
                    break

            logger.debug(f"Found MD files in {location}")
            result = {
                'topology': str(topology),
                'trajectory_raw': str(trajectory),
                'task_dir': str(task_dir)
            }

            if structure_pdb:
                result['structure_pdb'] = structure_pdb

            return result

        return None

    def _create_context(self,
                       task_name: str,
                       task_data: Dict,
                       output_base_dir: Optional[str] = None) -> PipelineContext:
        """
        Create PipelineContext from task data.

        Args:
            task_name: Task name (used as system_id)
            task_data: Task data dict
            output_base_dir: Base directory for outputs

        Returns:
            PipelineContext instance
        """
        # Determine output directory
        if output_base_dir:
            output_dir = str(Path(output_base_dir) / task_name)
        else:
            output_dir = f"./results/{task_name}"

        # Extract metadata
        metadata = {
            'task_dir': task_data.get('task_dir', ''),
            'discovery_method': 'default'
        }

        # Add PDB ID and variant if available
        if 'pdb_id' in task_data:
            metadata['pdb_id'] = task_data['pdb_id']
        if 'variant' in task_data:
            metadata['variant'] = task_data['variant']

        return PipelineContext(
            system_id=task_name,
            topology=task_data['topology'],
            trajectory_raw=task_data['trajectory_raw'],
            structure_pdb=task_data.get('structure_pdb'),
            output_dir=output_dir,
            metadata=metadata
        )

    def discover_from_list(self, task_list: List[Dict]) -> List[PipelineContext]:
        """
        Create contexts from a list of task dictionaries.

        Useful when you have a custom list of tasks (e.g., from a CSV file).

        Args:
            task_list: List of task dictionaries

        Returns:
            List of PipelineContext instances

        Example:
            >>> tasks = [
            ...     {
            ...         'system_id': '1ao7',
            ...         'topology': 'data/1ao7/md.tpr',
            ...         'trajectory_raw': 'data/1ao7/md.xtc'
            ...     },
            ...     {
            ...         'system_id': '1bd2',
            ...         'topology': 'data/1bd2/md.tpr',
            ...         'trajectory_raw': 'data/1bd2/md.xtc'
            ...     }
            ... ]
            >>> contexts = discovery.discover_from_list(tasks)
        """
        contexts = []

        for task_data in task_list:
            if 'system_id' not in task_data:
                logger.warning(f"Skipping task without system_id: {task_data}")
                continue

            if 'topology' not in task_data or 'trajectory_raw' not in task_data:
                logger.warning(f"Skipping task without required fields: {task_data['system_id']}")
                continue

            context = PipelineContext(**task_data)
            contexts.append(context)

        logger.info(f"Created {len(contexts)} contexts from task list")
        return contexts

    def __repr__(self) -> str:
        return f"TaskDiscovery(custom_rules={len(self.discovery_rules)})"
