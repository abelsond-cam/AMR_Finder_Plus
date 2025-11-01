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
import argparse
import pandas as pd

# Import the core processing function
from merge_amr_data import process_single_parquet

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(processName)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


def find_all_parquet_files(base_dir: Path, exclude_dirs: List[str] = None) -> List[Path]:
    """
    Recursively find all .parquet files in the base directory, excluding specified directories.
    
    Args:
        base_dir: Base directory to search
        exclude_dirs: List of directory names to exclude from search
    
    Returns:
        List of Path objects for all parquet files
    """
    if exclude_dirs is None:
        exclude_dirs = []
    
    parquet_files = []
    for parquet_file in base_dir.rglob('*.parquet'):
        # Check if file is in any excluded directory
        exclude_file = False
        for exclude_dir in exclude_dirs:
            if exclude_dir in parquet_file.parts:
                exclude_file = True
                break
        
        if not exclude_file:
            parquet_files.append(parquet_file)
    
    logger.info(f"Found {len(parquet_files)} parquet files in {base_dir} (excluding {exclude_dirs})")
    return sorted(parquet_files)


def process_file_wrapper(args: Tuple[Path, Path, Path, Path, int]) -> Tuple[Path, bool, str]:
    """
    Wrapper function for multiprocessing.
    
    Args:
        args: Tuple of (parquet_path, amrf_results_dir, output_dir, input_base_dir, timeout_seconds)
    
    Returns:
        Tuple of (parquet_path, success_flag, status_message)
    """
    import time
    start_time = time.time()
    parquet_path, amrf_results_dir, output_dir, input_base_dir, timeout_seconds = args
    
    try:
        success = process_single_parquet(
            parquet_path, 
            amrf_results_dir, 
            output_dir, 
            input_base_dir,
            timeout_seconds
        )
        elapsed = time.time() - start_time
        status = f"completed in {elapsed:.2f}s"
        return (parquet_path, success, status)
    except Exception as e:
        elapsed = time.time() - start_time
        error_msg = f"failed after {elapsed:.2f}s: {str(e)}"
        logger.error(f"Exception in worker processing {parquet_path}: {e}")
        return (parquet_path, False, error_msg)


def parse_args():
    """
    Parse command-line arguments.
    
    Returns:
        argparse.Namespace object with parsed arguments
    """
    # Define default paths
    default_input_dir = '/home/dca36/rds/hpc-work/data/BacFormer/processed/ast_esm_embeddings'
    default_amrf_dir = '/home/dca36/rds/hpc-work/data/amr_finder_plus/data/amrf_results'
    
    parser = argparse.ArgumentParser(
        description='Merge AMR Finder Plus results into parquet files using multiprocessing'
    )
    
    parser.add_argument(
        '--input-dir',
        type=str,
        default=default_input_dir,
        help=f'Input directory containing parquet files (default: {default_input_dir})'
    )
    
    parser.add_argument(
        '--amrf-dir',
        type=str,
        default=default_amrf_dir,
        help=f'Directory containing AMR result TSV files (default: {default_amrf_dir})'
    )
    
    parser.add_argument(
        '--output-dir',
        type=str,
        help='Output directory for processed files (default: same as input-dir)'
    )
    
    parser.add_argument(
        '--exclude-dirs',
        type=str,
        nargs='+',
        default=['ecoli'],
        help='Directory names to exclude from processing (default: ecoli)'
    )
    
    parser.add_argument(
        '--timeout',
        type=int,
        default=300,  # 5 minutes
        help='Maximum time in seconds to process a single file (default: 300s)'
    )
    
    parser.add_argument(
        '--resume',
        action='store_true',
        help='Resume processing - skip files that already have AMR data'
    )
    
    args = parser.parse_args()
    
    # Convert to Path objects
    args.input_dir = Path(args.input_dir)
    args.amrf_dir = Path(args.amrf_dir)
    
    # If output-dir not specified, use input-dir (overwrite mode)
    if args.output_dir is None:
        args.output_dir = args.input_dir
    else:
        args.output_dir = Path(args.output_dir)
    
    return args


def main():
    """
    Main function to coordinate multiprocessing pipeline.
    """
    # Parse command-line arguments
    args = parse_args()
    
    input_base_dir = args.input_dir
    amrf_results_dir = args.amrf_dir
    output_base_dir = args.output_dir
    exclude_dirs = args.exclude_dirs
    timeout_seconds = args.timeout
    resume_mode = args.resume
    
    logger.info("=" * 70)
    logger.info("MERGE AMR DATA - MULTIPROCESSING")
    logger.info("=" * 70)
    logger.info(f"Input directory: {input_base_dir}")
    logger.info(f"AMR results directory: {amrf_results_dir}")
    logger.info(f"Output directory: {output_base_dir}")
    logger.info(f"Excluding directories: {exclude_dirs}")
    logger.info(f"File processing timeout: {timeout_seconds}s")
    logger.info(f"Resume mode: {resume_mode}")
    logger.info("=" * 70)
    
    # Validate input directory exists
    if not input_base_dir.exists():
        logger.error(f"Input directory does not exist: {input_base_dir}")
        sys.exit(1)
    
    if not amrf_results_dir.exists():
        logger.error(f"AMR results directory does not exist: {amrf_results_dir}")
        sys.exit(1)
    
    # Create output base directory if different from input
    if output_base_dir != input_base_dir:
        output_base_dir.mkdir(parents=True, exist_ok=True)
    logger.info(f"Output directory: {output_base_dir}")
    
    # Find all parquet files
    parquet_files = find_all_parquet_files(input_base_dir, exclude_dirs=exclude_dirs)
    
    if not parquet_files:
        logger.warning("No parquet files found!")
        sys.exit(0)
    
    # Detect available CPUs
    n_cpus = multiprocessing.cpu_count()
    # Use n_cpus - 1 to leave one core free, with minimum of 1
    n_workers = max(1, n_cpus - 1)
    logger.info(f"Detected {n_cpus} CPUs, using {n_workers} workers")
    
    # Check if we need to filter files in resume mode
    if resume_mode:
        logger.info("Resume mode enabled - checking for already processed files...")
        filtered_files = []
        for pq_file in parquet_files:
            try:
                # Quick check if file has AMR columns already
                df = pd.read_parquet(pq_file, columns=[], engine='pyarrow')
                if 'element_symbol' in df.columns and 'element_name' in df.columns:
                    logger.info(f"Skipping already processed file: {pq_file.name}")
                    continue
                filtered_files.append(pq_file)
            except Exception as e:
                logger.warning(f"Error checking file {pq_file}: {e}")
                filtered_files.append(pq_file)
                
        skipped_count = len(parquet_files) - len(filtered_files)
        logger.info(f"Resume mode: {skipped_count} files already processed, {len(filtered_files)} files remaining")
        parquet_files = filtered_files
    
    # Prepare arguments for each file
    args_list = [
        (pq_file, amrf_results_dir, output_base_dir, input_base_dir, timeout_seconds)
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
            pq_path, success, status_msg = result
            status = "SUCCESS" if success else "FAILED"
            logger.info(f"[{i}/{len(parquet_files)}] {status}: {pq_path.name} - {status_msg}")
    
    # Summary statistics
    end_time = datetime.now()
    elapsed = end_time - start_time
    
    successes = sum(1 for _, success, _ in results if success)
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
        for pq_path, success, status_msg in results:
            if not success:
                logger.warning(f"  - {pq_path}: {status_msg}")
    
    # Exit with error code if any failures
    sys.exit(0 if failures == 0 else 1)


if __name__ == "__main__":
    # Set start method for multiprocessing (important for some systems)
    multiprocessing.set_start_method('spawn', force=True)
    main()

