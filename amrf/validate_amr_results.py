#!/usr/bin/env python
"""
Standalone AMR Results Validation Script

This script validates that all genomes from parquet files have corresponding 
AMR Finder Plus result files. Can be run independently to check processing completeness.
"""

import sys
import argparse
import pandas as pd
from pathlib import Path
from collections import defaultdict

def collect_expected_genomes(input_dir):
    """
    Collect all genome names from parquet files in the input directory.
    
    Args:
        input_dir: Path to directory containing parquet files
        
    Returns:
        tuple: (all_expected_genomes_set, file_stats_dict)
    """
    input_dir = Path(input_dir)
    parquet_files = sorted(input_dir.rglob("*.parquet"))
    
    if not parquet_files:
        print(f"Error: No .parquet files found in {input_dir}")
        return set(), {}
    
    all_expected_genomes = set()
    file_stats = {}
    
    print(f"\nScanning {len(parquet_files)} parquet files for genome names...")
    
    for pf in parquet_files:
        try:
            df = pd.read_parquet(pf, engine='pyarrow')
            
            if 'genome_name' in df.columns:
                unique_genomes = set(df['genome_name'].unique())
                genome_count = len(unique_genomes)
                all_expected_genomes.update(unique_genomes)
            else:
                unique_genomes = set()
                genome_count = len(df)  # Fallback
            
            file_stats[pf.name] = {
                'total_rows': len(df),
                'unique_genomes': unique_genomes,
                'genome_count': genome_count
            }
            
            print(f"  {pf.name}: {len(df)} rows, {genome_count} unique genomes")
            
        except Exception as e:
            print(f"  Warning: Could not read {pf.name}: {e}")
            file_stats[pf.name] = {'total_rows': 0, 'unique_genomes': set(), 'genome_count': 0}
    
    return all_expected_genomes, file_stats

def validate_amr_results(all_expected_genomes, results_dir):
    """
    Validate that all expected genomes have corresponding AMR result files.
    
    Args:
        all_expected_genomes: Set of all genome names that should be processed
        results_dir: Path to AMR results directory
        
    Returns:
        dict: Validation results
    """
    results_dir = Path(results_dir)
    
    if not results_dir.exists():
        print(f"Error: Results directory does not exist: {results_dir}")
        return None
    
    # Find all AMR result files
    amr_files = list(results_dir.glob("*_amr_results.tsv"))
    
    # Extract genome names from filenames
    genomes_with_files = set()
    file_sizes = {}
    
    for amr_file in amr_files:
        # Remove the "_amr_results.tsv" suffix to get genome name
        genome_name = amr_file.name.replace("_amr_results.tsv", "")
        genomes_with_files.add(genome_name)
        file_sizes[genome_name] = amr_file.stat().st_size
    
    # Calculate various sets for validation
    missing_files = all_expected_genomes - genomes_with_files
    unexpected_files = genomes_with_files - all_expected_genomes
    
    # Count files by size (to identify empty vs non-empty results)
    empty_files = {genome for genome, size in file_sizes.items() if size == 0}
    non_empty_files = genomes_with_files - empty_files
    
    validation_results = {
        'total_expected_genomes': len(all_expected_genomes),
        'total_amr_files_found': len(amr_files),
        'genomes_with_files': len(genomes_with_files),
        'missing_files': missing_files,
        'unexpected_files': unexpected_files,
        'empty_result_files': empty_files,
        'non_empty_result_files': non_empty_files,
        'file_sizes': file_sizes,
        'validation_success': len(missing_files) == 0
    }
    
    return validation_results

def main():
    parser = argparse.ArgumentParser(
        description='Validate AMR Finder Plus results completeness',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Validate results for acinetobacter_baumannii
  python validate_amr_results.py \\
    --input-dir /path/to/parquets/acinetobacter_baumannii \\
    --results-dir /path/to/amrf_results
  
  # Just check specific results directory against all species
  python validate_amr_results.py \\
    --input-dir /path/to/parquets/all_species \\
    --results-dir /path/to/amrf_results
        """
    )
    
    parser.add_argument(
        '--input-dir',
        type=str,
        required=True,
        help='Directory containing input parquet files'
    )
    parser.add_argument(
        '--results-dir',
        type=str,
        required=True,
        help='Directory containing AMR Finder Plus results'
    )
    parser.add_argument(
        '--show-details',
        action='store_true',
        help='Show detailed lists of missing/unexpected files'
    )
    
    args = parser.parse_args()
    
    print("="*80)
    print("AMR RESULTS VALIDATION")
    print("="*80)
    print(f"Input directory: {args.input_dir}")
    print(f"Results directory: {args.results_dir}")
    
    # Collect expected genomes
    all_expected_genomes, file_stats = collect_expected_genomes(args.input_dir)
    
    if not all_expected_genomes:
        print("No genomes found to validate.")
        sys.exit(1)
    
    # Validate results
    validation_results = validate_amr_results(all_expected_genomes, args.results_dir)
    
    if validation_results is None:
        sys.exit(1)
    
    print(f"\n" + "="*80)
    print("VALIDATION RESULTS")
    print("="*80)
    
    print(f"\nSummary:")
    print(f"  Expected genomes: {validation_results['total_expected_genomes']}")
    print(f"  AMR result files found: {validation_results['total_amr_files_found']}")
    print(f"  Genomes with result files: {validation_results['genomes_with_files']}")
    print(f"  Non-empty result files: {len(validation_results['non_empty_result_files'])}")
    print(f"  Empty result files: {len(validation_results['empty_result_files'])}")
    
    # Report missing files
    if validation_results['missing_files']:
        print(f"\n⚠ MISSING AMR FILES ({len(validation_results['missing_files'])}):")
        missing_list = sorted(list(validation_results['missing_files']))
        
        if args.show_details or len(missing_list) <= 20:
            for genome_name in missing_list:
                print(f"  - {genome_name}")
        else:
            for genome_name in missing_list[:20]:
                print(f"  - {genome_name}")
            print(f"  ... and {len(missing_list) - 20} more (use --show-details to see all)")
    
    # Report unexpected files
    if validation_results['unexpected_files']:
        print(f"\n⚠ UNEXPECTED AMR FILES ({len(validation_results['unexpected_files'])}):")
        unexpected_list = sorted(list(validation_results['unexpected_files']))
        
        if args.show_details or len(unexpected_list) <= 10:
            for genome_name in unexpected_list:
                print(f"  - {genome_name}")
        else:
            for genome_name in unexpected_list[:10]:
                print(f"  - {genome_name}")
            print(f"  ... and {len(unexpected_list) - 10} more (use --show-details to see all)")
    
    # Report empty files
    if validation_results['empty_result_files']:
        print(f"\n⚠ EMPTY RESULT FILES ({len(validation_results['empty_result_files'])}):")
        empty_list = sorted(list(validation_results['empty_result_files']))
        
        if args.show_details or len(empty_list) <= 10:
            for genome_name in empty_list:
                print(f"  - {genome_name}")
        else:
            for genome_name in empty_list[:10]:
                print(f"  - {genome_name}")
            print(f"  ... and {len(empty_list) - 10} more (use --show-details to see all)")
    
    # File statistics by parquet file
    print(f"\nFile Statistics by Parquet File:")
    for filename, stats in file_stats.items():
        print(f"  {filename}: {stats['genome_count']} genomes")
    
    # Final validation status
    print(f"\n" + "="*80)
    if validation_results['validation_success']:
        print("✅ VALIDATION PASSED - All expected genomes have AMR result files")
        completion_rate = len(validation_results['non_empty_result_files']) / validation_results['total_expected_genomes'] * 100
        print(f"📊 Completion rate: {completion_rate:.1f}% ({len(validation_results['non_empty_result_files'])}/{validation_results['total_expected_genomes']} with non-empty results)")
    else:
        print("❌ VALIDATION FAILED - Some genomes are missing AMR result files")
        missing_count = len(validation_results['missing_files'])
        completion_rate = (validation_results['total_expected_genomes'] - missing_count) / validation_results['total_expected_genomes'] * 100
        print(f"📊 Completion rate: {completion_rate:.1f}% ({validation_results['total_expected_genomes'] - missing_count}/{validation_results['total_expected_genomes']} completed)")
    
    print("="*80)
    
    # Exit codes: 0 = success, 1 = validation failed, 2 = error
    if validation_results['validation_success']:
        sys.exit(0)
    else:
        sys.exit(1)

if __name__ == "__main__":
    main()
