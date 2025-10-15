"""
Multiprocessing wrapper for merging AMR Finder Plus results into parquet files.

This script coordinates parallel processing of multiple parquet files using
multiprocessing to efficiently handle ~145 files.
"""

import multiprocessing
from pathlib import Path
import logging
from typing import List, Tuple
import sys
from datetime import datetime

# Import the core processing function
from merge_amr_data import process_single_parquet

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(processName)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


def find_all_parquet_files(base_dir: Path) -> List[Path]:
    """
    Recursively find all .parquet files in the base directory.
    
    Args:
        base_dir: Base directory to search
    
    Returns:
        List of Path objects for all parquet files
    """
    parquet_files = list(base_dir.rglob('*.parquet'))
    logger.info(f"Found {len(parquet_files)} parquet files in {base_dir}")
    return sorted(parquet_files)


def process_file_wrapper(args: Tuple[Path, Path, Path]) -> Tuple[Path, bool]:
    """
    Wrapper function for multiprocessing.
    
    Args:
        args: Tuple of (parquet_path, amrf_results_dir, output_dir)
    
    Returns:
        Tuple of (parquet_path, success_flag)
    """
    parquet_path, amrf_results_dir, output_dir = args
    try:
        success = process_single_parquet(parquet_path, amrf_results_dir, output_dir)
        return (parquet_path, success)
    except Exception as e:
        logger.error(f"Exception in worker processing {parquet_path}: {e}")
        return (parquet_path, False)


def main():
    """
    Main function to coordinate multiprocessing pipeline.
    """
    # Define paths
    input_base_dir = Path(
        '/home/dca36/rds/hpc-work/data/BacFormer/processed/'
        'ecoli_embeddings_with_ast_labels_flattened'
    )
    amrf_results_dir = Path('/home/dca36/rds/hpc-work/data/amr_finder_plus/data/amrf_results')
    output_base_dir = Path(
        '/home/dca36/rds/hpc-work/data/BacFormer/processed/'
        'ecoli_embeddings_with_amrf_labels_flattened'
    )
    
    # Validate input directory exists
    if not input_base_dir.exists():
        logger.error(f"Input directory does not exist: {input_base_dir}")
        sys.exit(1)
    
    if not amrf_results_dir.exists():
        logger.error(f"AMR results directory does not exist: {amrf_results_dir}")
        sys.exit(1)
    
    # Create output base directory
    output_base_dir.mkdir(parents=True, exist_ok=True)
    logger.info(f"Output directory: {output_base_dir}")
    
    # Find all parquet files
    parquet_files = find_all_parquet_files(input_base_dir)
    
    if not parquet_files:
        logger.warning("No parquet files found!")
        sys.exit(0)
    
    # Detect available CPUs
    n_cpus = multiprocessing.cpu_count()
    # Use n_cpus - 1 to leave one core free, with minimum of 1
    n_workers = max(1, n_cpus - 1)
    logger.info(f"Detected {n_cpus} CPUs, using {n_workers} workers")
    
    # Prepare arguments for each file
    args_list = [
        (pq_file, amrf_results_dir, output_base_dir)
        for pq_file in parquet_files
    ]
    
    # Process files in parallel
    logger.info(f"Starting parallel processing of {len(parquet_files)} files...")
    start_time = datetime.now()
    
    # Use multiprocessing Pool
    with multiprocessing.Pool(processes=n_workers) as pool:
        # Use imap_unordered for better performance and progress tracking
        results = []
        for i, result in enumerate(pool.imap_unordered(process_file_wrapper, args_list), 1):
            results.append(result)
            pq_path, success = result
            status = "SUCCESS" if success else "FAILED"
            logger.info(f"[{i}/{len(parquet_files)}] {status}: {pq_path.name}")
    
    # Summary statistics
    end_time = datetime.now()
    elapsed = end_time - start_time
    
    successes = sum(1 for _, success in results if success)
    failures = len(results) - successes
    
    logger.info("=" * 70)
    logger.info("PROCESSING COMPLETE")
    logger.info("=" * 70)
    logger.info(f"Total files: {len(results)}")
    logger.info(f"Successful: {successes}")
    logger.info(f"Failed: {failures}")
    logger.info(f"Time elapsed: {elapsed}")
    logger.info("=" * 70)
    
    if failures > 0:
        logger.warning("Failed files:")
        for pq_path, success in results:
            if not success:
                logger.warning(f"  - {pq_path}")
    
    # Exit with error code if any failures
    sys.exit(0 if failures == 0 else 1)


if __name__ == "__main__":
    # Set start method for multiprocessing (important for some systems)
    multiprocessing.set_start_method('spawn', force=True)
    main()

