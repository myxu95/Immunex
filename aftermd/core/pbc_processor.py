import subprocess
import os
from pathlib import Path
from typing import List, Optional, Dict, Any
import logging
from ..utils.group_selector import GroupSelector
from ..utils.cleanup_manager import CleanupManager

logger = logging.getLogger(__name__)


class PBCProcessor:
    """PBC (Periodic Boundary Conditions) processing for MD trajectories."""
    
    def __init__(self, gmx_executable: str = "gmx", keep_temp_files: bool = False):
        """
        Initialize PBCProcessor.
        
        Args:
            gmx_executable: GROMACS executable command
            keep_temp_files: Whether to keep temporary files for debugging (default: False)
        """
        self.gmx = gmx_executable
        self.group_selector = None  # Will be initialized when topology is provided
        self.keep_temp_files = keep_temp_files
        self.cleanup_manager = CleanupManager(keep_temp_files=keep_temp_files, verbose=True)
        self._check_gromacs()
    
    def _check_gromacs(self):
        """Check if GROMACS is available."""
        try:
            result = subprocess.run([self.gmx, "--version"], 
                                  capture_output=True, text=True, check=True)
            logger.info("GROMACS found and accessible")
        except (subprocess.CalledProcessError, FileNotFoundError):
            logger.warning("GROMACS not found or not accessible")
    
    def _register_temp_file(self, filepath: str) -> str:
        """Register a temporary file for later cleanup."""
        return self.cleanup_manager.register_file(filepath)
    
    def _cleanup_temp_files(self, output_dir: Optional[str] = None):
        """Clean up all registered temporary files and apply pattern-based cleanup."""
        # Clean registered files
        self.cleanup_manager.cleanup_registered_files()
        
        # Apply pattern-based cleanup in the output directory if provided
        if output_dir and not self.keep_temp_files:
            self.cleanup_manager.cleanup_by_patterns(
                directory=output_dir,
                categories=['gromacs_temp', 'index_files']
            )
    
    def _run_gmx_command(self, 
                        command: List[str], 
                        stdin_input: Optional[str] = None,
                        cwd: Optional[str] = None) -> subprocess.CompletedProcess:
        """
        Run GROMACS command with error handling.
        
        Args:
            command: GROMACS command list
            stdin_input: Input to pass to stdin
            cwd: Working directory
            
        Returns:
            CompletedProcess result
        """
        full_command = [self.gmx] + command
        logger.info(f"Running: {' '.join(full_command)}")
        
        try:
            result = subprocess.run(
                full_command,
                input=stdin_input,
                text=True,
                capture_output=True,
                check=True,
                cwd=cwd
            )
            return result
        except subprocess.CalledProcessError as e:
            logger.error(f"GROMACS command failed: {e}")
            logger.error(f"stdout: {e.stdout}")
            logger.error(f"stderr: {e.stderr}")
            raise
    
    def remove_pbc_2step(self,
                         trajectory: str,
                         topology: str,
                         output: str,
                         fit_group: Optional[str] = None,
                         output_group: Optional[str] = None,
                         dt: Optional[float] = None) -> str:
        """
        Remove PBC using simplified 2-step process: nojump -> fit rot+trans.

        This is the recommended method for large complexes (e.g., TCR-pMHC) based on
        empirical experience. It avoids potential issues with centering and whole operations.

        Args:
            trajectory: Input trajectory file (.xtc or .trr)
            topology: Topology file (.tpr or .gro)
            output: Output trajectory file
            fit_group: Group for fitting rot+trans (default: Backbone)
            output_group: Group to output (default: System)
            dt: Time interval for frame sampling in ps (default: no downsampling)

        Returns:
            Path to output file
        """
        # Initialize group selector for this topology
        if self.group_selector is None:
            self.group_selector = GroupSelector(topology, self.gmx)

        # Auto-select groups if not provided
        if fit_group is None:
            fit_group = self.group_selector.select_group("fit_backbone")
        if output_group is None:
            output_group = self.group_selector.select_group("output_system")

        # Create and register temporary file
        temp_nojump = self._register_temp_file(str(Path(output).with_suffix('.temp_nojump.xtc')))

        try:
            # Step 1: PBC nojump
            logger.info("Step 1: Applying PBC nojump to prevent atom jumps")
            nojump_cmd = [
                "trjconv",
                "-f", trajectory,
                "-s", topology,
                "-o", temp_nojump,
                "-pbc", "nojump"
            ]

            if dt is not None:
                nojump_cmd.extend(["-dt", str(dt)])
                logger.info(f"Using dt={dt} ps for trajectory downsampling")

            nojump_input = f"{output_group}\n"
            self._run_gmx_command(nojump_cmd, stdin_input=nojump_input)
            logger.info("PBC nojump completed")

            # Step 2: Fit rot+trans
            logger.info("Step 2: Fitting rotational and translational motion")
            fit_cmd = [
                "trjconv",
                "-f", temp_nojump,
                "-s", topology,
                "-o", output,
                "-fit", "rot+trans"
            ]

            fit_input = f"{fit_group}\n{output_group}\n"
            self._run_gmx_command(fit_cmd, stdin_input=fit_input)
            logger.info(f"Applied rot+trans fitting using group: {fit_group}")

            logger.info(f"2-step PBC processing completed: {output}")
            logger.info(f"Used groups - Fit: {fit_group}, Output: {output_group}")
            return output

        finally:
            # Clean up temporary files
            output_dir = str(Path(output).parent)
            self._cleanup_temp_files(output_dir)

    def remove_pbc(self,
                   trajectory: str,
                   topology: str,
                   output: str,
                   center_group: Optional[str] = None,
                   output_group: Optional[str] = None,
                   fit_group: Optional[str] = None,
                   dt: Optional[float] = None,
                   use_nojump: bool = False) -> str:
        """
        Remove PBC using three-step process: center -> pbc whole -> fit rot+trans.

        Note: The 2-step method (remove_pbc_2step) is recommended for most cases.
        Use this 3-step method only when you specifically need centering and whole operations.

        Args:
            trajectory: Input trajectory file (.xtc or .trr)
            topology: Topology file (.tpr or .gro)
            output: Output trajectory file
            center_group: Group to center (default: Protein)
            output_group: Group to output (default: System)
            fit_group: Group for fitting rot+trans (default: Protein)
            dt: Time interval for frame sampling in ps (default: no downsampling)
            use_nojump: Enable -pbc nojump to prevent atom jumps across boundaries

        Returns:
            Path to output file
        """
        # Initialize group selector for this topology
        if self.group_selector is None:
            self.group_selector = GroupSelector(topology, self.gmx)
        
        # Auto-select groups if not provided
        if center_group is None:
            center_group = self.group_selector.select_group("center_protein")
        if output_group is None:
            output_group = self.group_selector.select_group("output_system")
        if fit_group is None:
            fit_group = self.group_selector.select_group("fit_backbone")
        
        # Create and register temporary files
        temp_centered = self._register_temp_file(str(Path(output).with_suffix('.temp_centered.xtc')))
        temp_whole = self._register_temp_file(str(Path(output).with_suffix('.temp_whole.xtc')))
        
        try:
            # Step 1: Center trajectory with PBC
            pbc_method = "nojump" if use_nojump else "atom"
            if dt is not None:
                logger.info(f"Step 1: Centering trajectory (center pbc {pbc_method}) with dt={dt} ps")
            else:
                logger.info(f"Step 1: Centering trajectory (center pbc {pbc_method})")

            center_cmd = [
                "trjconv",
                "-f", trajectory,
                "-s", topology,
                "-o", temp_centered,
                "-center",
                "-pbc", pbc_method
            ]
            
            # Add index file if shortest chain is available
            if self.group_selector and self.group_selector.has_shortest_chain_group():
                index_file = self.group_selector.get_shortest_chain_index_file()
                if index_file:
                    center_cmd.extend(["-n", index_file])
                    logger.info(f"Using custom index file for centering: {index_file}")
            
            # Add dt option if specified
            if dt is not None:
                center_cmd.extend(["-dt", str(dt)])
            
            center_input = f"{center_group}\n{output_group}\n"
            self._run_gmx_command(center_cmd, stdin_input=center_input)
            logger.info(f"Centered trajectory using group: {center_group}")
            
            # Step 2: Apply PBC whole
            logger.info("Step 2: Applying PBC whole")
            whole_cmd = [
                "trjconv",
                "-f", temp_centered,
                "-s", topology,
                "-o", temp_whole,
                "-pbc", "whole"
            ]
            
            # Note: dt not needed for step 2 since we're processing the already sampled trajectory
            
            whole_input = f"{output_group}\n"
            self._run_gmx_command(whole_cmd, stdin_input=whole_input)
            logger.info("Applied PBC whole")
            
            # Step 3: Fit rot+trans
            logger.info("Step 3: Fitting rotational and translational motion")
            fit_cmd = [
                "trjconv",
                "-f", temp_whole,
                "-s", topology,
                "-o", output,
                "-fit", "rot+trans"
            ]
            
            # Note: dt not needed for step 3 since we're processing the already sampled trajectory
            
            fit_input = f"{fit_group}\n{output_group}\n"
            self._run_gmx_command(fit_cmd, stdin_input=fit_input)
            logger.info(f"Applied rot+trans fitting using group: {fit_group}")
            
            logger.info(f"PBC processing completed: {output}")
            logger.info(f"Used groups - Center: {center_group}, Output: {output_group}, Fit: {fit_group}")
            return output
            
        finally:
            # Register shortest chain index file for cleanup if it exists
            if self.group_selector and self.group_selector.has_shortest_chain_group():
                index_file = self.group_selector.get_shortest_chain_index_file()
                if index_file:
                    self._register_temp_file(index_file)
            
            # Clean up all registered temporary files and apply pattern cleanup
            output_dir = str(Path(output).parent)
            self._cleanup_temp_files(output_dir)

    def comprehensive_pbc_process(self,
                                trajectory: str,
                                topology: str,
                                output_dir: str,
                                method: str = "2step",
                                center_group: Optional[str] = None,
                                fit_group: Optional[str] = None,
                                dt: Optional[float] = None,
                                auto_dt: bool = True,
                                use_nojump: bool = True) -> Dict[str, str]:
        """
        Comprehensive PBC processing workflow.

        Args:
            trajectory: Input trajectory file
            topology: Topology file
            output_dir: Output directory
            method: PBC processing method (default: "2step")
                - "2step": nojump -> fit (recommended, based on empirical experience)
                - "3step": center -> whole -> fit (standard GROMACS workflow)
            center_group: Group for centering (only for 3step, default: auto-select)
            fit_group: Group for fitting (default: auto-select)
            dt: Time interval for frame sampling in ps (default: 100 ps)
            auto_dt: Automatically use default dt=100ps (default: True)
            use_nojump: Enable -pbc nojump in 3step method (default: True, only for 3step)

        Returns:
            Dictionary with paths to output files

        Note:
            The 2-step method is recommended for large complexes (e.g., TCR-pMHC)
            based on empirical experience. It avoids potential issues with centering
            and whole operations while effectively removing PBC artifacts.

            Standard: 100 ns trajectory -> 1000 frames (dt = 100 ps)
        """
        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)

        base_name = Path(trajectory).stem

        # Initialize group selector
        if self.group_selector is None:
            self.group_selector = GroupSelector(topology, self.gmx)

        # Validate method
        if method not in ["2step", "3step"]:
            raise ValueError(f"Invalid method: {method}. Must be '2step' or '3step'")

        # Auto-calculate dt if not provided
        if dt is None and auto_dt:
            dt = 100.0
            logger.info(f"Using default dt={dt} ps for downsampling (targets ~1000 frames per 100 ns)")
            logger.info(f"  To disable downsampling, use --no-auto-dt")
            logger.info(f"  To customize, use --dt <value_in_ps>")

        # Auto-select fit group if not provided
        if fit_group is None:
            fit_group = self.group_selector.select_group("fit_backbone")

        # Prepare output file
        processed_file = output_path / f"{base_name}_processed.xtc"

        logger.info(f"Starting comprehensive PBC processing for {trajectory}")
        logger.info(f"PBC processing method: {method}")
        if dt is not None:
            logger.info(f"Using dt={dt} ps for trajectory downsampling")

        # Execute processing based on method
        if method == "2step":
            logger.info("Using 2-step method: nojump -> fit (recommended)")
            logger.info(f"Using groups - Fit: {fit_group}")

            final_trajectory = self.remove_pbc_2step(
                trajectory=trajectory,
                topology=topology,
                output=str(processed_file),
                fit_group=fit_group,
                dt=dt
            )
        else:  # 3step
            # Auto-select center group for 3step
            if center_group is None:
                center_group = self.group_selector.select_group("center_protein")

            logger.info("Using 3-step method: center -> whole -> fit")
            logger.info(f"Using groups - Center: {center_group}, Fit: {fit_group}")

            final_trajectory = self.remove_pbc(
                trajectory=trajectory,
                topology=topology,
                output=str(processed_file),
                center_group=center_group,
                fit_group=fit_group,
                dt=dt,
                use_nojump=use_nojump
            )

        # Copy reference files (md.tpr and md.gro) for RMSD calculation and visualization
        copied_files = self._copy_reference_files(trajectory, topology, output_path)

        results = {
            "processed": final_trajectory,
            "method": method,
            "fit_group": fit_group,
            "output_directory": str(output_path),
            "reference_structures": copied_files["gro_files"],
            "reference_topology": copied_files["tpr_file"]
        }

        # Add center_group for 3step method
        if method == "3step":
            results["center_group"] = center_group

        logger.info(f"Comprehensive PBC processing completed in: {output_dir}")
        return results
    
    def _copy_reference_files(self, trajectory: str, topology: str, output_dir: Path) -> Dict[str, Any]:
        """
        Copy md.tpr and md.gro files from source directory to output directory.

        Args:
            trajectory: Path to trajectory file
            topology: Path to topology file
            output_dir: Output directory path

        Returns:
            Dictionary with 'tpr_file' and 'gro_files' lists
        """
        import shutil

        # Get the source directory (where the trajectory file is located)
        traj_path = Path(trajectory)
        source_dir = traj_path.parent

        # Also check the topology file's directory (in case files are in different locations)
        topo_path = Path(topology)
        search_dirs = [source_dir, topo_path.parent]

        # Remove duplicates while preserving order
        search_dirs = list(dict.fromkeys(search_dirs))

        copied_gro_files = []
        copied_tpr_file = None

        # Common .gro file names to look for
        gro_candidates = [
            "md.gro",
            "prod.gro",
            "production.gro",
            f"{traj_path.stem}.gro"
        ]

        # Common .tpr file names to look for
        tpr_candidates = [
            "md.tpr",
            "prod.tpr",
            "production.tpr",
            f"{traj_path.stem}.tpr"
        ]

        # Copy .tpr file (highest priority for RMSD reference)
        for search_dir in search_dirs:
            for tpr_name in tpr_candidates:
                tpr_file = search_dir / tpr_name

                if tpr_file.exists() and tpr_file.is_file():
                    try:
                        dest_name = "md.tpr"
                        dest_path = output_dir / dest_name

                        if not dest_path.exists():
                            shutil.copy2(tpr_file, dest_path)
                            copied_tpr_file = str(dest_path)
                            logger.info(f"Copied TPR file: {tpr_file} -> {dest_path}")
                            break

                    except Exception as e:
                        logger.warning(f"Failed to copy {tpr_file}: {e}")

            if copied_tpr_file:
                break

        # Copy .gro file
        for search_dir in search_dirs:
            for gro_name in gro_candidates:
                gro_file = search_dir / gro_name

                if gro_file.exists() and gro_file.is_file():
                    try:
                        if gro_name == "md.gro":
                            dest_name = "md.gro"
                        else:
                            dest_name = f"reference_{gro_name}"

                        dest_path = output_dir / dest_name

                        if dest_path.exists():
                            continue

                        shutil.copy2(gro_file, dest_path)
                        copied_gro_files.append(str(dest_path))
                        logger.info(f"Copied GRO file: {gro_file} -> {dest_path}")

                        if gro_name == "md.gro":
                            break

                    except Exception as e:
                        logger.warning(f"Failed to copy {gro_file}: {e}")

            if copied_gro_files and copied_gro_files[0].endswith("md.gro"):
                break

        # Log summary
        if copied_tpr_file:
            logger.info(f"Successfully copied TPR reference file for RMSD calculation")
        else:
            logger.warning("No .tpr file found - RMSD will use last frame as reference")
            for search_dir in search_dirs:
                logger.info(f"  Searched in: {search_dir}")

        if copied_gro_files:
            logger.info(f"Successfully copied {len(copied_gro_files)} GRO structure file(s)")
        else:
            logger.info("No .gro structure file found for trajectory visualization")

        return {
            "tpr_file": copied_tpr_file,
            "gro_files": copied_gro_files
        }