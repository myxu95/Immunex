import os
import glob
from pathlib import Path
from typing import List, Callable, Any, Dict, Optional, Tuple
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor
import logging

logger = logging.getLogger(__name__)


# Global function for multiprocessing compatibility
def _process_md_task_global(task_info: Tuple[str, Tuple[str, str]],
                           output_base_dir: str, dt: Optional[float],
                           method: str, gmx_executable: str) -> Dict[str, Any]:
    """
    Global function to process a single MD task - required for pickle compatibility.
    """
    from ..preprocessing.pbc_processor import PBCProcessor

    task_name, (trajectory, topology) = task_info

    try:
        logger.info(f"Processing task: {task_name}")

        # Create task-specific output directory
        task_output_dir = Path(output_base_dir) / task_name
        task_output_dir.mkdir(parents=True, exist_ok=True)

        # Initialize PBC processor
        pbc_processor = PBCProcessor(gmx_executable)

        # Run comprehensive PBC processing
        result = pbc_processor.comprehensive_pbc_process(
            trajectory=trajectory,
            topology=topology,
            output_dir=str(task_output_dir),
            method=method,
            dt=dt
        )

        result["task_name"] = task_name
        result["status"] = "success"
        logger.info(f"Successfully processed task: {task_name}")
        return result

    except Exception as e:
        logger.error(f"Failed to process task {task_name}: {e}")
        return {
            "task_name": task_name,
            "status": "failed",
            "error": str(e)
        }


# Global function for CDR recognition - required for pickle compatibility
def _process_cdr_task_global(task_info: Tuple[str, str],
                            output_base_dir: str,
                            chain_selections: Dict[str, str]) -> Dict[str, Any]:
    """
    Global function to process a single CDR recognition task.

    Args:
        task_info: Tuple of (task_name, topology_file)
        output_base_dir: Base output directory
        chain_selections: Chain selection dictionary

    Returns:
        Task processing result dictionary
    """
    from .cdr_manager import CDRManager
    from datetime import datetime
    import traceback

    task_name, topology_file = task_info
    start_time = datetime.now()

    try:
        logger.info(f"Processing CDR recognition: {task_name}")

        # Create task-specific output directory
        task_output_dir = Path(output_base_dir) / task_name
        task_output_dir.mkdir(parents=True, exist_ok=True)

        # Initialize CDR manager
        manager = CDRManager(
            topology_file=topology_file,
            allow_fallback=True
        )

        # Detect CDR regions
        cdr_data = manager.detect_all_cdr_regions(chains=chain_selections)

        # Generate index files
        # CDR3 only
        cdr3_index = task_output_dir / "cdr3_only.ndx"
        manager.generate_index_file(
            output_file=str(cdr3_index),
            include_cdrs=[3]
        )

        # All CDRs
        all_cdr_index = task_output_dir / "cdr_regions.ndx"
        manager.generate_index_file(
            output_file=str(all_cdr_index),
            include_cdrs=[1, 2, 3]
        )

        # Save metadata
        metadata_file = task_output_dir / "cdr_metadata.json"
        manager.save_metadata(
            output_file=str(metadata_file),
            task_name=task_name
        )

        # Get summary
        summary = manager.get_cdr_summary()

        processing_time = (datetime.now() - start_time).total_seconds()

        result = {
            "task_name": task_name,
            "status": "success",
            "processing_time_sec": processing_time,
            "detection_method": summary['detection_method'],
            "chains_processed": summary['total_chains'],
            "output_dir": str(task_output_dir)
        }

        # Add CDR3 sequences if detected
        for chain_name, chain_summary in summary['chains'].items():
            if chain_summary.get('status') == 'success':
                result[f"{chain_name}_cdr3"] = chain_summary.get('cdr3_sequence', 'N/A')
                result[f"{chain_name}_cdrs"] = ','.join(chain_summary.get('cdrs_detected', []))

        logger.info(f"Successfully processed CDR recognition: {task_name} ({processing_time:.1f}s)")
        return result

    except Exception as e:
        processing_time = (datetime.now() - start_time).total_seconds()
        logger.error(f"Failed to process CDR recognition {task_name}: {e}")
        logger.debug(traceback.format_exc())
        return {
            "task_name": task_name,
            "status": "failed",
            "error": str(e),
            "processing_time_sec": processing_time
        }


class BatchProcessor:
    """Batch processing utility for handling multiple files efficiently."""
    
    def __init__(self, max_workers: int = None):
        """
        Initialize BatchProcessor.
        
        Args:
            max_workers: Maximum number of worker processes/threads
        """
        self.max_workers = max_workers
        
    def find_files(self, directory: str, pattern: str) -> List[str]:
        """
        Find files matching pattern in directory.
        
        Args:
            directory: Directory path to search
            pattern: File pattern (e.g., "*.xtc", "*.pdb")
            
        Returns:
            List of matching file paths
        """
        search_path = os.path.join(directory, pattern)
        files = glob.glob(search_path, recursive=True)
        return sorted(files)
    
    def find_md_input_files(self, task_path: str) -> Optional[Tuple[str, str]]:
        """
        Find MD input files (md.xtc and md.tpr) using specific discovery rules.
        
        Discovery rules:
        1. Look directly in the task path
        2. Look in "prod" subfolder under task path
        
        Args:
            task_path: Path to the task directory
            
        Returns:
            Tuple of (trajectory_file, topology_file) if both found, None otherwise
        """
        task_dir = Path(task_path)
        
        if not task_dir.exists():
            logger.warning(f"Task directory does not exist: {task_path}")
            return None
        
        # Define search locations in priority order
        search_locations = [
            task_dir,              # Rule 1: directly in task path
            task_dir / "prod"      # Rule 2: in prod subfolder
        ]
        
        for location in search_locations:
            if not location.exists():
                continue
                
            xtc_file = location / "md.xtc"
            tpr_file = location / "md.tpr"
            
            if xtc_file.exists() and tpr_file.exists():
                logger.info(f"Found MD files in {location}: md.xtc, md.tpr")
                return str(xtc_file), str(tpr_file)
            else:
                missing = []
                if not xtc_file.exists():
                    missing.append("md.xtc")
                if not tpr_file.exists():
                    missing.append("md.tpr")
                logger.debug(f"Location {location}: missing {', '.join(missing)}")
        
        logger.warning(f"No complete MD input files (md.xtc + md.tpr) found in task: {task_path}")
        return None
    
    def discover_batch_tasks(self, base_directory: str) -> Dict[str, Tuple[str, str]]:
        """
        Discover all valid MD tasks in a base directory.
        
        Args:
            base_directory: Base directory containing multiple task folders
            
        Returns:
            Dictionary mapping task names to (trajectory, topology) file pairs
        """
        base_path = Path(base_directory)
        
        if not base_path.exists():
            logger.error(f"Base directory does not exist: {base_directory}")
            return {}
        
        discovered_tasks = {}
        
        # Look for subdirectories that might contain MD files
        for item in base_path.iterdir():
            if item.is_dir():
                task_name = item.name
                md_files = self.find_md_input_files(str(item))

                if md_files:
                    trajectory, topology = md_files
                    discovered_tasks[task_name] = (trajectory, topology)
                    logger.info(f"Task '{task_name}': found MD files")
                else:
                    logger.debug(f"Task '{task_name}': no MD files found")
        
        logger.info(f"Discovered {len(discovered_tasks)} valid MD tasks in {base_directory}")
        return discovered_tasks
    
    def process_files(self, 
                     files: List[str], 
                     process_func: Callable,
                     use_multiprocessing: bool = True,
                     **kwargs) -> List[Any]:
        """
        Process multiple files using parallel execution.
        
        Args:
            files: List of file paths to process
            process_func: Function to apply to each file
            use_multiprocessing: Use ProcessPoolExecutor if True, ThreadPoolExecutor if False
            **kwargs: Additional arguments to pass to process_func
            
        Returns:
            List of results from processing each file
        """
        executor_class = ProcessPoolExecutor if use_multiprocessing else ThreadPoolExecutor
        
        with executor_class(max_workers=self.max_workers) as executor:
            futures = [executor.submit(process_func, file_path, **kwargs) for file_path in files]
            results = []
            
            for future in futures:
                try:
                    result = future.result()
                    results.append(result)
                except Exception as e:
                    logger.error(f"Error processing file: {e}")
                    results.append(None)
                    
        return results
    
    def batch_analyze(self, 
                     input_pattern: str,
                     output_dir: str,
                     analysis_func: Callable,
                     **analysis_kwargs) -> Dict[str, Any]:
        """
        Perform batch analysis on files matching pattern.
        
        Args:
            input_pattern: Pattern to find input files (e.g., "/path/to/*.xtc")
            output_dir: Directory to save results
            analysis_func: Analysis function to apply
            **analysis_kwargs: Additional arguments for analysis function
            
        Returns:
            Dictionary with analysis results summary
        """
        files = glob.glob(input_pattern)
        if not files:
            raise ValueError(f"No files found matching pattern: {input_pattern}")
            
        Path(output_dir).mkdir(parents=True, exist_ok=True)
        
        results = self.process_files(files, analysis_func, 
                                   output_dir=output_dir, **analysis_kwargs)
        
        successful = sum(1 for r in results if r is not None)
        total = len(files)
        
        summary = {
            "total_files": total,
            "successful": successful,
            "failed": total - successful,
            "results": results
        }
        
        logger.info(f"Batch analysis completed: {successful}/{total} files processed successfully")
        return summary
    
    def batch_pbc_processing(self,
                           base_directory: str,
                           output_base_dir: str,
                           method: str = "2step",
                           dt: Optional[float] = None,
                           gmx_executable: str = "gmx") -> Dict[str, Any]:
        """
        Perform batch PBC processing with shortest chain detection on multiple MD tasks.

        Args:
            base_directory: Base directory containing task folders
            output_base_dir: Base directory for output files
            method: PBC processing method ("2step" or "3step", default: "2step")
            dt: Time interval for frame sampling in ps (optional)
            gmx_executable: GROMACS executable command

        Returns:
            Dictionary with batch processing results summary
        """
        from ..preprocessing.pbc_processor import PBCProcessor

        # Discover all valid MD tasks
        discovered_tasks = self.discover_batch_tasks(base_directory)

        if not discovered_tasks:
            logger.error("No valid MD tasks found for batch processing")
            return {"total_tasks": 0, "successful": 0, "failed": 0, "results": []}

        # Create output base directory
        Path(output_base_dir).mkdir(parents=True, exist_ok=True)

        # Prepare task list for parallel processing
        task_list = list(discovered_tasks.items())

        # Use custom parallel processing for tasks with global function
        executor_class = ProcessPoolExecutor
        with executor_class(max_workers=self.max_workers) as executor:
            futures = [
                executor.submit(_process_md_task_global, task_info, output_base_dir, dt, method, gmx_executable)
                for task_info in task_list
            ]
            results = []

            for future in futures:
                try:
                    result = future.result()
                    results.append(result)
                except Exception as e:
                    logger.error(f"Error in batch processing: {e}")
                    results.append({"status": "failed", "error": str(e)})
        
        # Summarize results
        successful = sum(1 for r in results if r.get("status") == "success")
        failed = len(results) - successful
        
        summary = {
            "total_tasks": len(discovered_tasks),
            "successful": successful,
            "failed": failed,
            "results": results,
            "output_directory": output_base_dir
        }
        
        logger.info(f"Batch PBC processing completed: {successful}/{len(discovered_tasks)} tasks processed successfully")

        # Log task details
        for result in results:
            task_name = result.get("task_name", "unknown")
            status = result.get("status", "unknown")
            if status == "success":
                logger.info(f"✅ {task_name}: PBC processing completed")
            else:
                error = result.get("error", "unknown error")
                logger.error(f"❌ {task_name}: {error}")

        return summary

    def batch_cdr_recognition(self,
                            base_directory: str,
                            output_dir: str,
                            chain_selections: Optional[Dict[str, str]] = None,
                            num_workers: Optional[int] = None) -> Dict[str, Any]:
        """
        Batch CDR recognition for multiple MD tasks.

        Args:
            base_directory: Base directory containing task folders with topology files
            output_dir: Output directory for CDR analysis results
            chain_selections: Dictionary mapping chain names to selection strings
                            Default: {'TCR_alpha': 'chainID D', 'TCR_beta': 'chainID E'}
            num_workers: Number of parallel workers (overrides self.max_workers)

        Returns:
            Dictionary with batch processing results summary
        """
        from datetime import datetime
        import pandas as pd

        # Default chain selections for pHLA-TCR complexes
        if chain_selections is None:
            chain_selections = {
                'TCR_alpha': 'chainID D and protein',
                'TCR_beta': 'chainID E and protein'
            }

        logger.info("="*60)
        logger.info("Batch CDR Recognition")
        logger.info("="*60)
        logger.info(f"Base directory: {base_directory}")
        logger.info(f"Output directory: {output_dir}")
        logger.info(f"Chain selections: {chain_selections}")
        logger.info("="*60)

        # Create output directory
        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)

        # Discover tasks - look for topology files
        base_path = Path(base_directory)

        if not base_path.exists():
            logger.error(f"Base directory does not exist: {base_directory}")
            return {"total_tasks": 0, "successful": 0, "failed": 0, "results": []}

        discovered_tasks = {}

        # Look for task directories
        for item in sorted(base_path.iterdir()):
            if not item.is_dir():
                continue

            task_name = item.name

            # Look for topology file (.tpr, .pdb, or .gro)
            topology_file = None
            for ext in ['.tpr', '.pdb', '.gro']:
                # Try different patterns
                candidates = [
                    item / f"md{ext}",
                    item / f"{task_name}{ext}",
                    item / f"system{ext}"
                ]

                for candidate in candidates:
                    if candidate.exists():
                        topology_file = str(candidate)
                        break

                if topology_file:
                    break

            if topology_file:
                discovered_tasks[task_name] = topology_file
                logger.info(f"Task '{task_name}': found topology {Path(topology_file).name}")
            else:
                logger.debug(f"Task '{task_name}': no topology file found")

        if not discovered_tasks:
            logger.error("No valid tasks found for batch CDR recognition")
            return {"total_tasks": 0, "successful": 0, "failed": 0, "results": []}

        logger.info(f"\nDiscovered {len(discovered_tasks)} tasks for processing")

        # Prepare task list
        task_list = list(discovered_tasks.items())

        # Parallel processing
        workers = num_workers if num_workers else self.max_workers
        start_time = datetime.now()

        with ProcessPoolExecutor(max_workers=workers) as executor:
            futures = [
                executor.submit(_process_cdr_task_global, task_info, output_dir, chain_selections)
                for task_info in task_list
            ]
            results = []

            for future in futures:
                try:
                    result = future.result()
                    results.append(result)
                except Exception as e:
                    logger.error(f"Error in batch CDR processing: {e}")
                    results.append({"status": "failed", "error": str(e)})

        # Calculate statistics
        successful = sum(1 for r in results if r.get("status") == "success")
        failed = len(results) - successful
        total_time = (datetime.now() - start_time).total_seconds()

        # Create summary CSV
        summary_csv = output_path / "batch_summary.csv"
        try:
            df = pd.DataFrame(results)
            df.to_csv(summary_csv, index=False)
            logger.info(f"\nBatch summary saved: {summary_csv}")
        except Exception as e:
            logger.warning(f"Failed to save summary CSV: {e}")
            summary_csv = None

        # Create failed tasks log
        failed_tasks = [r for r in results if r.get("status") == "failed"]
        if failed_tasks:
            failed_log = output_path / "failed_tasks.log"
            with open(failed_log, 'w') as f:
                f.write(f"Failed CDR Recognition Tasks\n")
                f.write(f"Total failed: {len(failed_tasks)}\n")
                f.write(f"="*60 + "\n\n")
                for task in failed_tasks:
                    f.write(f"Task: {task.get('task_name', 'unknown')}\n")
                    f.write(f"Error: {task.get('error', 'unknown')}\n")
                    f.write("-"*60 + "\n\n")
            logger.info(f"Failed tasks log: {failed_log}")

        # Print summary
        logger.info("\n" + "="*60)
        logger.info("Batch CDR Recognition Summary")
        logger.info("="*60)
        logger.info(f"Total tasks: {len(discovered_tasks)}")
        logger.info(f"Successful: {successful}")
        logger.info(f"Failed: {failed}")
        logger.info(f"Success rate: {successful/len(discovered_tasks)*100:.1f}%")
        logger.info(f"Total time: {total_time:.1f}s")
        logger.info(f"Average time per task: {total_time/len(discovered_tasks):.1f}s")
        logger.info("="*60)

        # Log task details
        for result in results:
            task_name = result.get("task_name", "unknown")
            status = result.get("status", "unknown")
            if status == "success":
                time = result.get("processing_time_sec", 0)
                method = result.get("detection_method", "unknown")
                logger.info(f"✓ {task_name}: {method} ({time:.1f}s)")
            else:
                error = result.get("error", "unknown")
                logger.error(f"✗ {task_name}: {error}")

        return {
            "total_tasks": len(discovered_tasks),
            "successful": successful,
            "failed": failed,
            "results": results,
            "output_directory": output_dir,
            "summary_csv": str(summary_csv) if summary_csv else None,
            "total_time_sec": total_time
        }
