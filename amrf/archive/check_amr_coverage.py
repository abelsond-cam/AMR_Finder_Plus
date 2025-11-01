#!/usr/bin/env python3
"""
AMR Coverage Checker - Debug missing AMR result files.

This script checks if all genomes in the raw protein sequence files have 
corresponding AMR result files. Enhanced --check-only mode provides detailed
debugging to validate species and locate genomes in parquet files.
"""

import pandas as pd
from pathlib import Path
import logging
import argparse
import sys
from typing import List, Set, Dict
from tqdm import tqdm
from amrf.run_amrf_genome_file import process_parquet_file

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


def find_raw_parquet_files(species_dir: Path) -> List[Path]:
    """
    Find all parquet files in the raw protein sequences directory.
    
    Args:
        species_dir: Path to species directory in raw/ast/protein_sequences/
        
    Returns:
        List of parquet file paths
    """
    if not species_dir.exists():
        logger.error(f"Species directory does not exist: {species_dir}")
        return []
    
    parquet_files = list(species_dir.glob('*.parquet'))
    logger.info(f"Found {len(parquet_files)} raw parquet files in {species_dir}")
    
    return sorted(parquet_files)


def extract_genome_names_from_parquet(parquet_path: Path) -> Set[str]:
    """
    Extract all unique genome names from a parquet file.
    
    Args:
        parquet_path: Path to parquet file
        
    Returns:
        Set of genome names found in the file
    """
    try:
        # Read only the genome_name column for efficiency
        df = pd.read_parquet(parquet_path, columns=['genome_name'], engine='pyarrow')
        genome_names = set(df['genome_name'].unique())
        logger.debug(f"Extracted {len(genome_names)} unique genomes from {parquet_path.name}")
        return genome_names
    
    except Exception as e:
        logger.error(f"Failed to read {parquet_path}: {e}")
        return set()


def get_all_genome_names(raw_parquet_files: List[Path]) -> Set[str]:
    """
    Extract all unique genome names from all parquet files.
    
    Args:
        raw_parquet_files: List of parquet file paths
        
    Returns:
        Set of all unique genome names
    """
    all_genomes = set()
    
    for i, parquet_file in enumerate(raw_parquet_files, 1):
        logger.info(f"Processing file {i}/{len(raw_parquet_files)}: {parquet_file.name}")
        genomes = extract_genome_names_from_parquet(parquet_file)
        all_genomes.update(genomes)
        logger.info(f"  - Found {len(genomes)} genomes, total unique: {len(all_genomes)}")
    
    logger.info(f"Total unique genomes across all files: {len(all_genomes)}")
    return all_genomes


def check_amr_files_exist(genome_names: Set[str], amrf_results_dir: Path) -> Dict[str, bool]:
    """
    Check which genomes have corresponding AMR result files.
    
    Args:
        genome_names: Set of genome names to check
        amrf_results_dir: Directory containing AMR result files
        
    Returns:
        Dictionary mapping genome_name -> file_exists (bool)
    """
    logger.info(f"Checking AMR file availability for {len(genome_names)} genomes...")
    
    amr_status = {}
    missing_count = 0
    
    for genome_name in sorted(genome_names):
        amr_file_path = amrf_results_dir / f"{genome_name}_amr_results.tsv"
        exists = amr_file_path.exists()
        amr_status[genome_name] = exists
        
        if not exists:
            missing_count += 1
            if missing_count <= 10:  # Show first 10 missing files
                logger.warning(f"Missing AMR file: {amr_file_path}")
            elif missing_count == 11:
                logger.warning("... (additional missing files not shown)")
    
    available_count = len(genome_names) - missing_count
    logger.info(f"AMR file status: {available_count} available, {missing_count} missing")
    
    return amr_status


def validate_species_in_binary_labels(missing_genomes: List[str], expected_species: str) -> Dict[str, str]:
    """
    Check species for missing genomes in binary_labels.csv file.
    
    Args:
        missing_genomes: List of genome names missing AMR files
        expected_species: Expected species name (e.g., 'Campylobacter jejuni')
        
    Returns:
        Dictionary mapping genome_name -> actual_species
    """
    binary_labels_file = Path('/home/dca36/rds/hpc-work/data/BacFormer/raw/ast/binary_labels.csv')
    
    if not binary_labels_file.exists():
        logger.error(f"Binary labels file not found: {binary_labels_file}")
        return {}
    
    logger.info(f"Checking species for {len(missing_genomes)} missing genomes in binary_labels.csv...")
    
    try:
        # Read the binary labels file
        df = pd.read_csv(binary_labels_file)
        logger.info(f"Loaded binary_labels.csv with {len(df)} rows and columns: {df.columns.tolist()}")
        
        # Check if Species column exists
        if 'Species' not in df.columns:
            logger.error(f"'Species' column not found in binary_labels.csv. Available columns: {df.columns.tolist()}")
            return {}
        
        # Filter for missing genomes
        genome_species = {}
        found_count = 0
        not_found_count = 0
        wrong_species_count = 0
        
        for genome_name in missing_genomes:
            # Look for the genome in the dataframe
            genome_rows = df[df['genome_name'] == genome_name]
            
            if len(genome_rows) == 0:
                logger.warning(f"Genome {genome_name} not found in binary_labels.csv")
                genome_species[genome_name] = "NOT_FOUND"
                not_found_count += 1
            else:
                species = genome_rows.iloc[0]['Species']
                genome_species[genome_name] = species
                found_count += 1
                
                if species != expected_species:
                    wrong_species_count += 1
                    logger.warning(f"Genome {genome_name} has species '{species}', expected '{expected_species}'")
        
        logger.info("Species validation results:")
        logger.info(f"  - Found in binary_labels.csv: {found_count}")
        logger.info(f"  - Not found in binary_labels.csv: {not_found_count}")
        logger.info(f"  - Wrong species: {wrong_species_count}")
        logger.info(f"  - Correct species: {found_count - wrong_species_count}")
        
        return genome_species
        
    except Exception as e:
        logger.error(f"Error reading binary_labels.csv: {e}")
        return {}


def check_parquet_amrf_columns(raw_parquet_files: List[Path]) -> Dict[Path, bool]:
    """
    Check which parquet files have AMRF columns computed.
    
    Args:
        raw_parquet_files: List of parquet file paths to check
        
    Returns:
        Dictionary mapping parquet_file_path -> has_amrf_columns (bool)
    """
    logger.info(f"Checking AMRF columns in {len(raw_parquet_files)} parquet files...")
    
    # The 7 AMRF columns that indicate computation
    amrf_columns = [
        'element_symbol', 'element_name', 'amrf_type', 'amrf_subtype',
        'amrf_class', 'amrf_subclass', 'pct_identity_to_reference'
    ]
    
    file_status = {}
    
    for i, parquet_file in enumerate(raw_parquet_files, 1):
        try:
            # Efficiently read only the schema/columns without loading data
            # Read a single row to get column information
            df_sample = pd.read_parquet(parquet_file, engine='pyarrow', nrows=1)
            file_columns = set(df_sample.columns)
            
            # Check if all AMRF columns are present
            has_all_columns = all(col in file_columns for col in amrf_columns)
            file_status[parquet_file] = has_all_columns
            
            if i % 50 == 0:  # Progress update every 50 files
                logger.debug(f"Checked {i}/{len(raw_parquet_files)} files...")
                
        except Exception as e:
            logger.error(f"Error checking columns in {parquet_file.name}: {e}")
            file_status[parquet_file] = False
            continue
    
    with_columns = sum(1 for has_cols in file_status.values() if has_cols)
    without_columns = len(file_status) - with_columns
    
    logger.info(f"AMRF column check complete: {with_columns} files with columns, {without_columns} files without columns")
    
    return file_status


def print_amrf_column_examples(
    file_status: Dict[Path, bool],
    num_examples: int = 5
) -> None:
    """
    Print examples of files with and without AMRF columns.
    
    Args:
        file_status: Dictionary mapping parquet_file_path -> has_amrf_columns (bool)
        num_examples: Number of examples to show for each category (default: 5)
    """
    files_with_columns = [f for f, has_cols in file_status.items() if has_cols]
    files_without_columns = [f for f, has_cols in file_status.items() if not has_cols]
    
    logger.info("=" * 80)
    logger.info("AMRF COLUMN EXAMPLES")
    logger.info("=" * 80)
    
    # Files with AMRF columns
    if files_with_columns:
        logger.info(f"Files WITH AMRF columns ({len(files_with_columns)} total):")
        for i, parquet_file in enumerate(files_with_columns[:num_examples], 1):
            logger.info(f"  {i}. {parquet_file.resolve()}")
        if len(files_with_columns) > num_examples:
            logger.info(f"  ... and {len(files_with_columns) - num_examples} more files with columns")
    else:
        logger.warning("No files found with AMRF columns")
    
    logger.info("")
    
    # Files without AMRF columns
    if files_without_columns:
        logger.info(f"Files WITHOUT AMRF columns ({len(files_without_columns)} total):")
        for i, parquet_file in enumerate(files_without_columns[:num_examples], 1):
            logger.info(f"  {i}. {parquet_file.resolve()}")
        if len(files_without_columns) > num_examples:
            logger.info(f"  ... and {len(files_without_columns) - num_examples} more files without columns")
    else:
        logger.info("All files have AMRF columns")
    
    logger.info("=" * 80)


def locate_genomes_in_parquet_files(
    missing_genomes: List[str], 
    raw_parquet_files: List[Path]
) -> Dict[str, List[str]]:
    """
    Find which parquet files contain the missing genomes.
    
    Args:
        missing_genomes: List of genome names to locate
        raw_parquet_files: List of raw parquet files to search
        
    Returns:
        Dictionary mapping genome_name -> list_of_files_containing_genome
    """
    logger.info(f"Searching for {len(missing_genomes)} missing genomes in {len(raw_parquet_files)} parquet files...")
    
    genome_locations = {genome: [] for genome in missing_genomes}
    total_found = 0
    
    for i, parquet_file in enumerate(raw_parquet_files, 1):
        logger.info(f"Searching file {i}/{len(raw_parquet_files)}: {parquet_file.name}")
        
        try:
            # Read only genome_name column for efficiency
            df = pd.read_parquet(parquet_file, columns=['genome_name'], engine='pyarrow')
            file_genomes = set(df['genome_name'].unique())
            logger.info(f"  - File contains {len(file_genomes)} unique genomes")
            
            # Check which missing genomes are in this file
            found_in_file = []
            for genome_name in missing_genomes:
                if genome_name in file_genomes:
                    genome_locations[genome_name].append(parquet_file.name)
                    found_in_file.append(genome_name)
            
            if found_in_file:
                logger.info(f"  - Found {len(found_in_file)} missing genomes: {found_in_file[:5]}{'...' if len(found_in_file) > 5 else ''}")
                total_found += len(found_in_file)
            
        except Exception as e:
            logger.error(f"Error reading {parquet_file}: {e}")
            continue
    
    # Summary
    found_genomes = [g for g, files in genome_locations.items() if files]
    not_found_genomes = [g for g, files in genome_locations.items() if not files]
    
    logger.info("=" * 60)
    logger.info("GENOME LOCATION SUMMARY")
    logger.info("=" * 60)
    logger.info(f"Missing genomes searched for: {len(missing_genomes)}")
    logger.info(f"Found in parquet files: {len(found_genomes)}")
    logger.info(f"Not found in any parquet file: {len(not_found_genomes)}")
    
    if not_found_genomes:
        logger.warning("Genomes NOT found in any parquet file:")
        for genome in not_found_genomes[:10]:  # Show first 10
            logger.warning(f"  - {genome}")
        if len(not_found_genomes) > 10:
            logger.warning(f"  ... and {len(not_found_genomes) - 10} more")
    
    if found_genomes:
        logger.info("Sample of found genomes:")
        for genome in found_genomes[:5]:
            files = genome_locations[genome]
            logger.info(f"  - {genome}: found in {len(files)} file(s) -> {files}")
    
    return genome_locations


def generate_missing_amr_files(
    genome_locations: Dict[str, List[str]], 
    raw_parquet_files: List[Path], 
    amrf_results_dir: Path
) -> None:
    """
    Generate AMR files for missing genomes by processing relevant parquet files.
    
    Args:
        genome_locations: Dictionary mapping genome_name -> list_of_files_containing_genome
        raw_parquet_files: List of raw parquet file paths
        amrf_results_dir: Directory for AMR results output
    """
    logger.info("=" * 80)
    logger.info("GENERATING MISSING AMR FILES")
    logger.info("=" * 80)
    
    # Build reverse mapping: parquet_file -> list_of_missing_genomes_in_file
    parquet_to_genomes = {}
    
    for genome_name, file_list in genome_locations.items():
        if not file_list:  # Skip genomes not found in any parquet file
            continue
            
        for file_name in file_list:
            if file_name not in parquet_to_genomes:
                parquet_to_genomes[file_name] = []
            parquet_to_genomes[file_name].append(genome_name)
    
    # Convert file names back to full paths
    parquet_file_map = {pf.name: pf for pf in raw_parquet_files}
    files_to_process = []
    
    for file_name, target_genomes in parquet_to_genomes.items():
        if file_name in parquet_file_map:
            files_to_process.append((parquet_file_map[file_name], target_genomes))
        else:
            logger.warning(f"Parquet file not found: {file_name}")
    
    logger.info(f"Will process {len(files_to_process)} parquet files containing missing genomes")
    
    # Track overall statistics
    total_genomes_processed = 0
    total_amr_hits = 0
    total_files_processed = 0
    
    # Process each parquet file with its target genomes
    for parquet_file, target_genomes in tqdm(files_to_process, desc="Processing parquet files"):
        logger.info(f"Processing {parquet_file.name} with {len(target_genomes)} target genomes")
        logger.info(f"Target genomes: {target_genomes[:5]}{'...' if len(target_genomes) > 5 else ''}")
        
        try:
            # Call the AMR processing function with target genomes
            output_dir, stats = process_parquet_file(
                input_file=parquet_file,
                target_genomes=set(target_genomes),  # Pass as set for efficient lookup
                verbose=False,  # Reduce verbosity since we're processing many files
                fasta_dir=None,  # Use default
                output_dir=amrf_results_dir
            )
            
            # Update statistics
            total_genomes_processed += stats['genomes_processed']
            total_amr_hits += stats['amr_hits']
            total_files_processed += 1
            
            logger.info(f"Completed {parquet_file.name}: "
                       f"{stats['genomes_processed']} genomes, "
                       f"{stats['amr_hits']} AMR hits")
                       
        except Exception as e:
            logger.error(f"Error processing {parquet_file.name}: {e}")
            continue
    
    # Final statistics
    logger.info("=" * 80)
    logger.info("AMR GENERATION COMPLETE")
    logger.info("=" * 80)
    logger.info(f"Files processed: {total_files_processed}/{len(files_to_process)}")
    logger.info(f"Total genomes processed: {total_genomes_processed}")
    logger.info(f"Total AMR hits generated: {total_amr_hits}")
    logger.info(f"AMR results saved to: {amrf_results_dir}")


def parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description='Check AMR coverage for all genomes and debug missing files'
    )
    
    parser.add_argument(
        '--species',
        type=str,
        required=True,
        choices=['acinetobacter_baumannii', 'campylobacter_jejuni'],
        help='Bacterial species to check'
    )
    
    parser.add_argument(
        '--raw-sequences-dir',
        type=str,
        default='/home/dca36/rds/hpc-work/data/BacFormer/raw/ast/protein_sequences',
        help='Base directory containing raw protein sequence files'
    )
    
    parser.add_argument(
        '--amrf-results-dir',
        type=str,
        default='/home/dca36/rds/hpc-work/data/amr_finder_plus/data/amrf_results',
        help='Directory containing AMR result files'
    )
    
    
    parser.add_argument(
        '--check-only',
        action='store_true',
        help='Enhanced debugging mode: validate species and locate genomes in parquet files'
    )
    
    
    return parser.parse_args()


def main():
    """Main function."""
    args = parse_args()
    
    # Convert to Path objects
    raw_sequences_base = Path(args.raw_sequences_dir)
    species_raw_dir = raw_sequences_base / args.species
    amrf_results_dir = Path(args.amrf_results_dir)
    
    logger.info("=" * 80)
    logger.info("AMR COVERAGE CHECKER")
    logger.info("=" * 80)
    logger.info(f"Species: {args.species}")
    logger.info(f"Raw sequences directory: {species_raw_dir}")
    logger.info(f"AMR results directory: {amrf_results_dir}")
    logger.info(f"Check only mode: {args.check_only}")
    logger.info("=" * 80)
    
    # Validate directories
    if not species_raw_dir.exists():
        logger.error(f"Species raw directory does not exist: {species_raw_dir}")
        sys.exit(1)
    
    if not amrf_results_dir.exists():
        logger.error(f"AMR results directory does not exist: {amrf_results_dir}")
        sys.exit(1)
    
    # Step 1: Find all raw parquet files
    raw_parquet_files = find_raw_parquet_files(species_raw_dir)
    if not raw_parquet_files:
        logger.error("No raw parquet files found!")
        sys.exit(1)
    
    # Step 2: Extract all genome names
    all_genome_names = get_all_genome_names(raw_parquet_files)
    if not all_genome_names:
        logger.error("No genome names found!")
        sys.exit(1)
    
    # Step 3: Check AMR file availability
    amr_status = check_amr_files_exist(all_genome_names, amrf_results_dir)
    
    missing_genomes = [genome for genome, exists in amr_status.items() if not exists]
    available_genomes = [genome for genome, exists in amr_status.items() if exists]
    
    logger.info("=" * 80)
    logger.info("COVERAGE SUMMARY")
    logger.info("=" * 80)
    logger.info(f"Total genomes: {len(all_genome_names)}")
    logger.info(f"AMR files available: {len(available_genomes)}")
    logger.info(f"AMR files missing: {len(missing_genomes)}")
    logger.info(f"Coverage: {len(available_genomes)/len(all_genome_names)*100:.1f}%")
    
    if missing_genomes:
        logger.warning("Missing AMR files for genomes:")
        for i, genome in enumerate(sorted(missing_genomes)[:10]):  # Show first 10
            logger.warning(f"  {i+1:3d}. {genome}")
        if len(missing_genomes) > 10:
            logger.warning(f"  ... and {len(missing_genomes) - 10} more")
    
    # Enhanced check-only mode with debugging steps
    if args.check_only:
        logger.info("=" * 80)
        logger.info("ENHANCED DEBUGGING MODE")
        logger.info("=" * 80)
        
        # Check AMRF columns in parquet files
        logger.info("Checking AMRF columns in parquet files...")
        file_status = check_parquet_amrf_columns(raw_parquet_files)
        print_amrf_column_examples(file_status, num_examples=5)
        
        if not missing_genomes:
            logger.info("=" * 80)
            logger.info("Check-only debugging complete. All genomes have AMR coverage!")
            logger.info("=" * 80)
            sys.exit(0)
        
        # Step 4a: Validate species in binary_labels.csv
        expected_species_map = {
            'campylobacter_jejuni': 'Campylobacter jejuni',
            'acinetobacter_baumannii': 'Acinetobacter baumannii'
        }
        expected_species = expected_species_map.get(args.species, args.species)
        
        logger.info(f"Step 1: Validating species for missing genomes (expecting '{expected_species}')...")
        species_validation = validate_species_in_binary_labels(missing_genomes, expected_species)
        
        # Check for wrong species
        wrong_species_genomes = [g for g, species in species_validation.items() 
                                if species != expected_species and species != "NOT_FOUND"]
        
        if wrong_species_genomes:
            logger.error("=" * 80)
            logger.error("WRONG SPECIES DETECTED!")
            logger.error("=" * 80)
            logger.error(f"Found {len(wrong_species_genomes)} genomes with wrong species:")
            
            species_counts = {}
            for genome in wrong_species_genomes:
                species = species_validation[genome]
                species_counts[species] = species_counts.get(species, 0) + 1
                if len([g for g in wrong_species_genomes if species_validation[g] == species]) <= 5:
                    logger.error(f"  - {genome}: {species}")
            
            logger.error("Species breakdown of wrong genomes:")
            for species, count in species_counts.items():
                logger.error(f"  - {species}: {count} genomes")
            
            logger.error("Stopping due to wrong species. Fix species assignment first.")
            sys.exit(1)
        
        # Step 4b: Locate genomes in parquet files  
        correct_species_genomes = [g for g, species in species_validation.items() 
                                  if species == expected_species]
        
        if correct_species_genomes:
            logger.info("=" * 80)
            logger.info(f"Step 2: Locating {len(correct_species_genomes)} correct-species genomes in parquet files...")
            genome_locations = locate_genomes_in_parquet_files(correct_species_genomes, raw_parquet_files)
            
            found_in_parquet = [g for g, files in genome_locations.items() if files]
            not_found_in_parquet = [g for g, files in genome_locations.items() if not files]
            
            logger.info("=" * 80)
            logger.info("FINAL DEBUGGING SUMMARY")
            logger.info("=" * 80)
            logger.info(f"Total missing AMR files: {len(missing_genomes)}")
            logger.info(f"Wrong species (excluded): {len(wrong_species_genomes)}")
            logger.info(f"Correct species: {len(correct_species_genomes)}")
            logger.info(f"Found in parquet files: {len(found_in_parquet)}")
            logger.info(f"Not found in parquet files: {len(not_found_in_parquet)}")
            
            if not_found_in_parquet:
                logger.warning("Genomes with correct species but NOT found in parquet files:")
                for genome in not_found_in_parquet[:10]:
                    logger.warning(f"  - {genome}")
                if len(not_found_in_parquet) > 10:
                    logger.warning(f"  ... and {len(not_found_in_parquet) - 10} more")
        
        logger.info("=" * 80)
        logger.info("Check-only debugging complete. Use the information above to understand the missing genomes.")
        sys.exit(0)
    
    # Step 4: Generate missing files if requested (NOT in check-only mode)
    elif not args.check_only and missing_genomes:
        logger.info("=" * 80)
        logger.info("AMR FILE GENERATION MODE")
        logger.info("=" * 80)
        
        # First validate species for missing genomes
        expected_species_map = {
            'campylobacter_jejuni': 'Campylobacter jejuni',
            'acinetobacter_baumannii': 'Acinetobacter baumannii'
        }
        expected_species = expected_species_map.get(args.species, args.species)
        
        logger.info(f"Validating species for {len(missing_genomes)} missing genomes (expecting '{expected_species}')...")
        species_validation = validate_species_in_binary_labels(missing_genomes, expected_species)
        
        # Filter to only correct species genomes
        correct_species_genomes = [g for g, species in species_validation.items() 
                                  if species == expected_species]
        
        if not correct_species_genomes:
            logger.error("No genomes with correct species found for processing")
            sys.exit(1)
            
        logger.info(f"Found {len(correct_species_genomes)} genomes with correct species for processing")
        
        # Locate genomes in parquet files
        logger.info("Locating genomes in parquet files...")
        genome_locations = locate_genomes_in_parquet_files(correct_species_genomes, raw_parquet_files)
        
        # Generate AMR files
        generate_missing_amr_files(genome_locations, raw_parquet_files, amrf_results_dir)
    
    logger.info("=" * 80)
    logger.info("AMR COVERAGE CHECK COMPLETE")
    logger.info("=" * 80)
    
    if not missing_genomes:
        logger.info("All genomes have AMR coverage!")


if __name__ == "__main__":
    main()
