#!/usr/bin/env python3
"""
Simplified single-threaded AMR data processor for parquet files.

This script processes parquet files for a specific bacterial species,
adding AMR data and creating flattened protein IDs as needed.
Removes multiprocessing complexity to avoid hanging issues.
"""

import pandas as pd
import numpy as np
from pathlib import Path
import logging
import argparse
import sys
import traceback
from typing import List, Dict
import time

# Import functions from existing module
from amrf.merge_amr_data import (
    flatten_protein_ids,
    load_amr_data,
    map_amr_to_proteins,
    validate_amr_data
)

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


def find_parquet_files(species_dir: Path) -> List[Path]:
    """
    Recursively find all .parquet files in the species directory.
    
    Args:
        species_dir: Directory path for the bacterial species
        
    Returns:
        List of Path objects for all parquet files
    """
    if not species_dir.exists():
        logger.error(f"Species directory does not exist: {species_dir}")
        return []
    
    parquet_files = list(species_dir.rglob('*.parquet'))
    logger.info(f"Found {len(parquet_files)} parquet files in {species_dir}")
    
    return sorted(parquet_files)


def check_existing_columns(df: pd.DataFrame) -> Dict[str, bool]:
    """
    Check which required columns exist in the dataframe.
    
    Args:
        df: Input dataframe
        
    Returns:
        Dictionary with column existence status
    """
    checks = {
        'has_amr_data': all(col in df.columns for col in [
            'element_symbol', 'element_name', 'amrf_type', 'amrf_subtype',
            'amrf_class', 'amrf_subclass', 'pct_identity_to_reference'
        ]),
        'has_protein_id_flattened': 'protein_id_flattened' in df.columns,
        'has_gene_names': 'gene_names' in df.columns,
        'has_gene_products': 'gene_products' in df.columns,
        'has_protein_id': 'protein_id' in df.columns
    }
    
    return checks


def create_protein_id_flattened(df: pd.DataFrame) -> pd.DataFrame:
    """
    Create flattened protein ID column from nested protein_id column.
    
    Args:
        df: Input dataframe with nested protein_id column
        
    Returns:
        DataFrame with added protein_id_flattened column
    """
    logger.info("Creating protein_id_flattened column...")
    
    protein_ids_flat_list = []
    
    for idx, row in df.iterrows():
        genome_name = row['genome_name']
        protein_ids_nested = row['protein_id']
        
        # Use existing flatten function
        protein_ids_flat = flatten_protein_ids(protein_ids_nested, genome_name)
        protein_ids_flat_list.append(protein_ids_flat)
    
    df['protein_id_flattened'] = protein_ids_flat_list
    logger.info(f"Created protein_id_flattened for {len(df)} genomes")
    
    return df


def process_single_file(
    parquet_path: Path,
    amrf_results_dir: Path,
    base_dir: Path
) -> bool:
    """
    Process a single parquet file, adding AMR data and flattened protein IDs as needed.
    
    Args:
        parquet_path: Path to the parquet file to process
        amrf_results_dir: Directory containing AMR result TSV files
        base_dir: Base directory for relative path calculations
        
    Returns:
        True if successful, False otherwise
    """
    start_time = time.time()
    
    try:
        logger.info(f"Processing file: {parquet_path}")
        logger.info(f"File size: {parquet_path.stat().st_size / (1024*1024):.2f} MB")
        
        # Read parquet file
        try:
            df = pd.read_parquet(parquet_path, engine='pyarrow')
            logger.info(f"Loaded dataframe with shape: {df.shape}")
            logger.info(f"Columns: {list(df.columns)}")
        except Exception as e:
            logger.error(f"CRITICAL: Failed to read parquet file {parquet_path}")
            logger.error(f"Error: {type(e).__name__}: {e}")
            logger.error(f"Traceback: {traceback.format_exc()}")
            return False
        
        # Check existing columns
        column_checks = check_existing_columns(df)
        logger.info(f"Column status: {column_checks}")
        
        # Check if AMR data already exists
        if column_checks['has_amr_data']:
            logger.info("AMR data columns already present - skipping AMR merge")
            if column_checks['has_protein_id_flattened']:
                logger.info("All required columns present - no processing needed")
                elapsed_time = time.time() - start_time
                logger.info(f"File check completed in {elapsed_time:.2f}s")
                return True
        
        # Verify required columns exist
        if not column_checks['has_protein_id']:
            logger.error(f"CRITICAL: Missing required 'protein_id' column in {parquet_path}")
            return False
        
        if not column_checks['has_gene_names']:
            logger.warning(f"Missing 'gene_names' column in {parquet_path}")
        
        if not column_checks['has_gene_products']:
            logger.warning(f"Missing 'gene_products' column in {parquet_path}")
        
        # Track what modifications we make
        modifications_made = []
        
        # Create protein_id_flattened if missing
        if not column_checks['has_protein_id_flattened']:
            logger.info("Creating protein_id_flattened column")
            df = create_protein_id_flattened(df)
            modifications_made.append("created protein_id_flattened")
        
        # Process AMR data if not already present
        if not column_checks['has_amr_data']:
            logger.info("Adding AMR data...")
            
            # Initialize AMR columns
            amr_columns = [
                'element_symbol', 'element_name', 'amrf_type', 'amrf_subtype',
                'amrf_class', 'amrf_subclass', 'pct_identity_to_reference'
            ]
            
            # Collect AMR data for all genomes
            all_amr_data = {col: [] for col in amr_columns}
            
            for idx, row in df.iterrows():
                genome_name = row['genome_name']
                protein_ids_flat = row['protein_id_flattened']
                
                logger.info(f"Processing AMR for genome {idx+1}/{len(df)}: {genome_name}")
                
                # Load AMR data
                try:
                    amr_df = load_amr_data(genome_name, amrf_results_dir)
                except SystemExit:
                    # load_amr_data calls sys.exit on critical errors
                    logger.error(f"Critical AMR loading error for {genome_name}")
                    return False
                
                # Map AMR data to proteins
                amr_data = map_amr_to_proteins(protein_ids_flat, amr_df, genome_name)
                
                # Validate mapping
                if not validate_amr_data(protein_ids_flat, amr_data, amr_df, genome_name):
                    logger.error(f"AMR validation failed for {genome_name}")
                    return False
                
                # Collect data for this genome
                for col_name, arr in amr_data.items():
                    all_amr_data[col_name].append(arr)
            
            # Add AMR columns to dataframe
            for col_name, data_list in all_amr_data.items():
                df[col_name] = data_list
            
            modifications_made.append("added AMR data")
        
        # Write back to file if modifications were made
        if modifications_made:
            logger.info(f"Writing modified file back to: {parquet_path}")
            logger.info(f"Modifications made: {', '.join(modifications_made)}")
            
            try:
                df.to_parquet(parquet_path, engine='pyarrow', index=False)
                logger.info(f"Successfully wrote file with new shape: {df.shape}")
            except Exception as e:
                logger.error(f"CRITICAL: Failed to write parquet file {parquet_path}")
                logger.error(f"Error: {type(e).__name__}: {e}")
                logger.error(f"Traceback: {traceback.format_exc()}")
                return False
        else:
            logger.info("No modifications needed")
        
        elapsed_time = time.time() - start_time
        logger.info(f"Successfully processed {parquet_path.name} in {elapsed_time:.2f}s")
        return True
        
    except Exception as e:
        elapsed_time = time.time() - start_time
        logger.error(f"CRITICAL EXCEPTION processing {parquet_path} after {elapsed_time:.2f}s")
        logger.error(f"Exception type: {type(e).__name__}")
        logger.error(f"Exception message: {e}")
        logger.error(f"Full traceback:")
        logger.error(traceback.format_exc())
        return False


def parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description='Process parquet files for a specific bacterial species, adding AMR data'
    )
    
    parser.add_argument(
        '--species',
        type=str,
        required=True,
        choices=['acinetobacter_baumannii', 'campylobacter_jejuni'],
        help='Bacterial species to process'
    )
    
    parser.add_argument(
        '--base-dir',
        type=str,
        default='/home/dca36/rds/hpc-work/data/BacFormer/processed/ast_esm_embeddings',
        help='Base directory containing species subdirectories'
    )
    
    parser.add_argument(
        '--amrf-dir',
        type=str,
        default='/home/dca36/rds/hpc-work/data/amr_finder_plus/data/amrf_results',
        help='Directory containing AMR result TSV files'
    )
    
    parser.add_argument(
        '--dry-run',
        action='store_true',
        help='List files that would be processed without actually processing them'
    )
    
    return parser.parse_args()


def main():
    """Main processing function."""
    args = parse_args()
    
    # Convert to Path objects
    base_dir = Path(args.base_dir)
    amrf_dir = Path(args.amrf_dir)
    species_dir = base_dir / args.species
    
    logger.info("=" * 80)
    logger.info("SIMPLIFIED AMR DATA PROCESSOR")
    logger.info("=" * 80)
    logger.info(f"Species: {args.species}")
    logger.info(f"Species directory: {species_dir}")
    logger.info(f"AMR results directory: {amrf_dir}")
    logger.info(f"Dry run mode: {args.dry_run}")
    logger.info("=" * 80)
    
    # Recommend running AMR coverage check first
    logger.info("IMPORTANT: Before processing, ensure all genomes have AMR data by running:")
    logger.info(f"  python check_amr_coverage.py --species {args.species} --check-only")
    logger.info("If missing AMR files are found, generate them with:")
    logger.info(f"  python check_amr_coverage.py --species {args.species}")
    logger.info("=" * 80)
    
    # Validate directories
    if not base_dir.exists():
        logger.error(f"Base directory does not exist: {base_dir}")
        sys.exit(1)
    
    if not species_dir.exists():
        logger.error(f"Species directory does not exist: {species_dir}")
        sys.exit(1)
    
    if not amrf_dir.exists():
        logger.error(f"AMR results directory does not exist: {amrf_dir}")
        sys.exit(1)
    
    # Find all parquet files
    parquet_files = find_parquet_files(species_dir)
    
    if not parquet_files:
        logger.warning("No parquet files found!")
        sys.exit(0)
    
    logger.info(f"Found {len(parquet_files)} files to process")
    
    # Dry run mode - just list files
    if args.dry_run:
        logger.info("DRY RUN MODE - Files that would be processed:")
        for i, pq_file in enumerate(parquet_files, 1):
            rel_path = pq_file.relative_to(species_dir)
            logger.info(f"  {i:3d}. {rel_path}")
        logger.info(f"Total: {len(parquet_files)} files")
        sys.exit(0)
    
    # Process files sequentially
    start_time = time.time()
    success_count = 0
    failure_count = 0
    failed_files = []
    
    for i, parquet_file in enumerate(parquet_files, 1):
        logger.info(f"\n[{i}/{len(parquet_files)}] Processing: {parquet_file.name}")
        
        success = process_single_file(parquet_file, amrf_dir, base_dir)
        
        if success:
            success_count += 1
            logger.info(f"✓ SUCCESS: {parquet_file.name}")
        else:
            failure_count += 1
            failed_files.append(parquet_file)
            logger.error(f"✗ FAILED: {parquet_file.name}")
            
            # Exit immediately on first failure for debugging
            logger.error("STOPPING PROCESSING DUE TO FAILURE")
            break
    
    # Final summary
    end_time = time.time()
    elapsed = end_time - start_time
    
    logger.info("=" * 80)
    logger.info("PROCESSING SUMMARY")
    logger.info("=" * 80)
    logger.info(f"Total files processed: {success_count + failure_count}")
    logger.info(f"Successful: {success_count}")
    logger.info(f"Failed: {failure_count}")
    logger.info(f"Total time: {elapsed:.2f} seconds")
    logger.info(f"Average time per file: {elapsed/(success_count + failure_count):.2f} seconds")
    
    if failed_files:
        logger.error("Failed files:")
        for failed_file in failed_files:
            logger.error(f"  - {failed_file}")
    
    logger.info("=" * 80)
    
    # Exit with error code if any failures
    sys.exit(0 if failure_count == 0 else 1)


if __name__ == "__main__":
    main()
