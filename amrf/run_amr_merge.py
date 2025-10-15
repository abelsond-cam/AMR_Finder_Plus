"""
Single-threaded wrapper for merging AMR Finder Plus results into parquet files.

This script processes parquet files sequentially (no multiprocessing).
Suitable for running on login nodes or when multiprocessing isn't needed.
"""

from pathlib import Path
import logging
from typing import List
import sys
from datetime import datetime

# Import the core processing function
from merge_amr_data import process_single_parquet

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
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


def main():
    """
    Main function to process files sequentially.
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
    
    logger.info(f"Processing {len(parquet_files)} files sequentially (single-threaded)...")
    
    # Process files one by one
    start_time = datetime.now()
    results = []
    
    for i, pq_file in enumerate(parquet_files, 1):
        logger.info(f"[{i}/{len(parquet_files)}] Processing: {pq_file.name}")
        
        try:
            success = process_single_parquet(pq_file, amrf_results_dir, output_base_dir)
            results.append((pq_file, success))
            
            status = "SUCCESS" if success else "FAILED"
            logger.info(f"[{i}/{len(parquet_files)}] {status}: {pq_file.name}")
            
        except Exception as e:
            logger.error(f"Exception processing {pq_file}: {e}", exc_info=True)
            results.append((pq_file, False))
    
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
    logger.info(f"Average time per file: {elapsed / len(results) if results else 'N/A'}")
    logger.info("=" * 70)
    
    if failures > 0:
        logger.warning("Failed files:")
        for pq_path, success in results:
            if not success:
                logger.warning(f"  - {pq_path}")
    
    # Exit with error code if any failures
    sys.exit(0 if failures == 0 else 1)


if __name__ == "__main__":
    main()

