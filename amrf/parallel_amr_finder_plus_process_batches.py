#!/usr/bin/env python
"""
Parallel processing script for AMR Finder Plus on HPC
Fixed version with better multiprocessing handling
"""

import sys
import argparse
from pathlib import Path
from multiprocessing import Pool, cpu_count
from datetime import datetime
from tqdm import tqdm

# Add workspace to path so we can import amrf module
WORKSPACE = Path("/home/dca36/workspace/AMR_Finder_Plus")
if str(WORKSPACE) not in sys.path:
    sys.path.insert(0, str(WORKSPACE))


def create_file_chunks(parquet_files, chunks_per_file):
    """
    Read parquet files once and split into dataframe chunks for parallel processing.
    
    Args:
        parquet_files: List of parquet file paths
        chunks_per_file: Number of chunks to split each file into
        
    Returns:
        tuple: (chunks_list, file_stats, all_expected_genomes) where:
            - chunks_list: List of tuples: (parquet_file_path, df_chunk, chunk_id, expected_genomes_in_chunk)
            - file_stats: Dict mapping filename -> {'total_genomes': int, 'unique_genomes': set}
            - all_expected_genomes: Set of all genome names that should be processed
    """
    import pandas as pd
    
    chunks = []
    file_stats = {}
    all_expected_genomes = set()
    
    print(f"\nReading and chunking {len(parquet_files)} files...")
    
    for pf in parquet_files:
        # Read parquet file ONCE
        try:
            df = pd.read_parquet(pf, engine='pyarrow')
            total_rows = len(df)
            
            # Extract genome names and track file statistics
            if 'genome_name' in df.columns:
                unique_genomes = set(df['genome_name'].unique())
                genome_count = len(unique_genomes)
                all_expected_genomes.update(unique_genomes)
            else:
                unique_genomes = set()
                genome_count = total_rows  # Fallback if no genome_name column
            
            file_stats[pf.name] = {
                'total_rows': total_rows,
                'unique_genomes': unique_genomes,
                'genome_count': genome_count
            }
            
            if chunks_per_file == 1 or total_rows <= chunks_per_file:
                # No chunking needed - use entire dataframe
                chunks.append((pf, df, None, unique_genomes))
                if total_rows <= chunks_per_file:
                    print(f"  {pf.name}: {total_rows} rows, {genome_count} genomes (too few to chunk - processing as single unit)")
                else:
                    print(f"  {pf.name}: {total_rows} rows, {genome_count} genomes (no chunking)")
            else:
                # Split dataframe into chunks, but only create non-empty chunks
                chunk_size = total_rows // chunks_per_file
                actual_chunks_created = 0
                for i in range(chunks_per_file):
                    start_idx = i * chunk_size
                    # Last chunk gets any remaining rows
                    end_idx = (i + 1) * chunk_size if i < chunks_per_file - 1 else total_rows
                    
                    # Only create chunk if it would have genomes
                    if start_idx < total_rows:
                        df_chunk = df.iloc[start_idx:end_idx].copy()
                        if len(df_chunk) > 0:  # Double-check chunk is not empty
                            # Get genome names in this specific chunk
                            chunk_genomes = set(df_chunk['genome_name'].unique()) if 'genome_name' in df_chunk.columns else set()
                            chunks.append((pf, df_chunk, i + 1, chunk_genomes))
                            actual_chunks_created += 1
                
                print(f"  {pf.name}: {total_rows} rows, {genome_count} genomes → {actual_chunks_created} chunks of ~{chunk_size} rows")
                
                # Delete the full dataframe to free memory
                del df
                
        except Exception as e:
            print(f"  Warning: Could not read {pf.name}: {e}")
            # Fall back to file-based processing (no pre-loaded chunk)
            file_stats[pf.name] = {'total_rows': 0, 'unique_genomes': set(), 'genome_count': 0}
            chunks.append((pf, None, None, set()))
    
    print(f"\nFile Statistics Summary:")
    total_expected_genomes = len(all_expected_genomes)
    total_files = len(parquet_files)
    total_chunks = len(chunks)
    
    for filename, stats in file_stats.items():
        print(f"  {filename}: {stats['total_rows']} rows, {stats['genome_count']} unique genomes")
    
    print(f"\nOverall Summary:")
    print(f"  Total parquet files: {total_files}")
    print(f"  Total unique genomes across all files: {total_expected_genomes}")
    print(f"  Total processing chunks created: {total_chunks}")
    print(f"  Expected AMR result files: {total_expected_genomes} genome-specific .tsv files\n")
    
    return chunks, file_stats, all_expected_genomes

def process_single_file_wrapper(args):
    """
    Wrapper function that imports and calls the processing function.
    This avoids pickling issues with multiprocessing.
    
    Args:
        args: Tuple of (parquet_file_path, log_file_path, fasta_dir, output_dir, 
                        df_chunk, chunk_id, species, expected_genomes_in_chunk)
        
    Returns:
        Tuple of (filename, success, message, stats, processed_genomes)
    """
    parquet_file, log_file, fasta_dir, output_dir, df_chunk, chunk_id, species, expected_genomes = args
    start_time = datetime.now()
    
    # Create display name
    file_display = Path(parquet_file).name
    if chunk_id is not None:
        file_display = f"{file_display} [chunk {chunk_id}]"
    
    try:
        # Add workspace to path in worker process
        import sys
        from pathlib import Path as PathType
        workspace = PathType("/home/dca36/workspace/AMR_Finder_Plus")
        if str(workspace) not in sys.path:
            sys.path.insert(0, str(workspace))
        
        # Import here to avoid pickling issues
        from amrf.run_amrf_genome_file import process_parquet_file
        
        # Log parquet file info
        chunk_info = f" (chunk {chunk_id})" if chunk_id is not None else ""
        genome_count = len(df_chunk) if df_chunk is not None else "unknown"
        print(f"\nProcessing: {parquet_file.name}{chunk_info} - {genome_count} genomes")
        
        # Call the function with pre-loaded chunk
        result_dir, stats = process_parquet_file(
            parquet_file,
            species=species,
            verbose=True, 
            log_file=log_file,
            fasta_dir=fasta_dir,
            output_dir=output_dir,
            chunk_id=chunk_id,
            df_chunk=df_chunk,  # Pass pre-loaded dataframe chunk
            debug_verbose=True
        )
        
        elapsed = (datetime.now() - start_time).total_seconds()
        
        # Determine which genomes were actually processed by checking if df_chunk exists and has genome_name
        if df_chunk is not None and 'genome_name' in df_chunk.columns:
            processed_genomes = set(df_chunk['genome_name'].unique())
        else:
            processed_genomes = expected_genomes  # Fallback - assume all expected were processed
        
        message = (f"Completed in {elapsed:.1f}s - "
                  f"{stats['genomes_processed']} processed, "
                  f"{stats['amr_hits']} AMR hits")
        
        return (file_display, True, message, stats, processed_genomes)
            
    except Exception as e:
        elapsed = (datetime.now() - start_time).total_seconds()
        import traceback
        error_msg = f"Error after {elapsed:.1f}s: {str(e)}\n{traceback.format_exc()}"
        return (file_display, False, error_msg, None, set())
    
def validate_amr_results(all_expected_genomes, results_dir, successfully_processed_genomes):
    """
    Validate that all expected genomes have corresponding AMR result files.
    
    Args:
        all_expected_genomes: Set of all genome names that should be processed
        results_dir: Path to AMR results directory
        successfully_processed_genomes: Set of genomes that were successfully processed
        
    Returns:
        dict: Validation results with missing files, unexpected files, etc.
    """
    results_dir = Path(results_dir)
    
    # Find all AMR result files
    amr_files = list(results_dir.glob("*_amr_results.tsv"))
    
    # Extract genome names from filenames
    genomes_with_files = set()
    for amr_file in amr_files:
        # Remove the "_amr_results.tsv" suffix to get genome name
        genome_name = amr_file.name.replace("_amr_results.tsv", "")
        genomes_with_files.add(genome_name)
    
    # Calculate various sets for validation
    missing_files = all_expected_genomes - genomes_with_files
    unexpected_files = genomes_with_files - all_expected_genomes
    successfully_processed_with_files = successfully_processed_genomes.intersection(genomes_with_files)
    successfully_processed_missing_files = successfully_processed_genomes - genomes_with_files
    
    validation_results = {
        'total_expected_genomes': len(all_expected_genomes),
        'total_successfully_processed': len(successfully_processed_genomes),
        'total_amr_files_found': len(amr_files),
        'genomes_with_files': len(genomes_with_files),
        'missing_files': missing_files,
        'unexpected_files': unexpected_files,
        'successfully_processed_with_files': successfully_processed_with_files,
        'successfully_processed_missing_files': successfully_processed_missing_files,
        'validation_success': len(missing_files) == 0 and len(successfully_processed_missing_files) == 0
    }
    
    return validation_results

def main():
    # Parse command line arguments
    parser = argparse.ArgumentParser(
        description='Parallel processing of parquet files through AMR Finder Plus',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Process default directory with auto-detected CPUs, splitting each file into 5 chunks
  python parallel_process_batches.py
  
  # Process custom directory with 76 CPUs, 5 chunks per file (optimal for 27 files)
  python parallel_process_batches.py --input-dir /path/to/parquets --species campylobacter_jejuni --cpus 76 --chunks-per-file 5
  
  # Specify custom output directories (no chunking)
  python parallel_process_batches.py --input-dir /path/to/input --species acinetobacter_baumannii \\
      --fasta-dir /path/to/fastas --results-dir /path/to/results --chunks-per-file 1
        """
    )
    
    parser.add_argument(
        '--input-dir', 
        type=str,
        default='/home/dca36/rds/hpc-work/data/BacFormer/raw/ast/protein_sequences/all_species',
        help='Directory containing input parquet files (default: %(default)s)'
    )
    parser.add_argument(
        '--fasta-dir',
        type=str,
        default='/home/dca36/rds/hpc-work/data/amr_finder_plus/data/fasta_files',
        help='Directory for temporary FASTA files (default: %(default)s)'
    )
    parser.add_argument(
        '--results-dir',
        type=str,
        default='/home/dca36/rds/hpc-work/data/amr_finder_plus/data/amrf_results',
        help='Directory for AMR Finder Plus results (default: %(default)s)'
    )
    parser.add_argument(
        '--log-dir',
        type=str,
        default='/home/dca36/rds/hpc-work/data/amr_finder_plus/data/logs',
        help='Directory for log files (default: %(default)s)'
    )
    parser.add_argument(
        '--cpus',
        type=int,
        default=None,
        help=f'Number of CPUs to use (default: auto-detect, currently {cpu_count()})'
    )
    parser.add_argument(
        '--chunks-per-file',
        type=int,
        default=5,
        help='Number of chunks to split each parquet file into for parallel processing (default: 5)'
    )
    parser.add_argument(
        '--species',
        type=str,
        required=True,
        help='Species name (lowercase with underscores, e.g., campylobacter_jejuni, acinetobacter_baumannii)'
    )
    
    args = parser.parse_args()
    
    # Convert to Path objects
    input_dir = Path(args.input_dir)
    fasta_dir = Path(args.fasta_dir)
    results_dir = Path(args.results_dir)
    log_dir = Path(args.log_dir)
    
    # Determine CPU count
    n_cpus = args.cpus if args.cpus else cpu_count()
    
    # Create directories
    fasta_dir.mkdir(parents=True, exist_ok=True)
    results_dir.mkdir(parents=True, exist_ok=True)
    log_dir.mkdir(parents=True, exist_ok=True)
    
    # Find all parquet files (recursively)
    parquet_files = sorted(input_dir.rglob("*.parquet"))
    
    if not parquet_files:
        print(f"Error: No .parquet files found in {input_dir}")
        sys.exit(1)
    
    # Create chunks from parquet files and collect tracking data
    file_chunks, file_stats, all_expected_genomes = create_file_chunks(parquet_files, args.chunks_per_file)
    
    print("="*80)
    print("AMR FINDER PLUS - PARALLEL BATCH PROCESSING")
    print("="*80)
    print(f"\nInput directory: {input_dir}")
    print(f"Species: {args.species}")
    print(f"FASTA directory: {fasta_dir}")
    print(f"Results directory: {results_dir}")
    print(f"Log directory: {log_dir}")
    print(f"Found {len(parquet_files)} parquet files")
    print(f"Chunks per file: {args.chunks_per_file}")
    print(f"Total work items: {len(file_chunks)}")
    print(f"Requested CPUs: {n_cpus}")
    print(f"Parallel workers that will be active: {min(n_cpus, len(file_chunks))}")
    print(f"\nStart time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("="*80)
    
    # Create work items for each chunk
    file_args = []
    for pf, df_chunk, chunk_id, expected_genomes_in_chunk in file_chunks:
        # Create unique log file name
        if chunk_id is not None:
            log_file = log_dir / f"{pf.stem}_chunk{chunk_id}_processing.log"
        else:
            log_file = log_dir / f"{pf.stem}_processing.log"
        
        # Clear existing log file
        if log_file.exists():
            log_file.unlink()
        
        file_args.append((pf, log_file, fasta_dir, results_dir, df_chunk, chunk_id, args.species, expected_genomes_in_chunk))
    
    # Process files in parallel
    start_time = datetime.now()
    results = []
    
    print(f"\nSpawning {n_cpus} worker processes...")
    
    with Pool(processes=n_cpus) as pool:
        print(f"Pool created with {n_cpus} workers")
        print("Starting parallel processing...\n")
        
        # Use imap_unordered for immediate results
        with tqdm(total=len(file_chunks), desc="Processing chunks", unit="chunk") as pbar:
            for result in pool.imap_unordered(process_single_file_wrapper, file_args, chunksize=1):
                results.append(result)
                filename, success, message, stats, processed_genomes = result
                status = "✓" if success else "✗"
                # Print detailed completion message
                print(f"\n{status} {filename}: {message}")
                pbar.set_postfix_str(f"{status} {filename}")
                pbar.update(1)
    
    # Process results and collect tracking data
    elapsed = (datetime.now() - start_time).total_seconds()
    
    print("\n" + "="*80)
    print("PROCESSING COMPLETE - ANALYZING RESULTS")
    print("="*80)
    print(f"\nTotal time: {elapsed/60:.1f} minutes ({elapsed:.1f} seconds)")
    print("\nDetailed Results:")
    
    success_count = 0
    failed_count = 0
    total_genomes_processed = 0
    total_amr_hits = 0
    all_successfully_processed_genomes = set()
    
    for filename, success, message, stats, processed_genomes in results:
        status = "✓" if success else "✗"
        print(f"  {status} {filename}: {message}")
        if success:
            success_count += 1
            if stats:
                total_genomes_processed += stats['genomes_processed']
                total_amr_hits += stats['amr_hits']
            all_successfully_processed_genomes.update(processed_genomes)
        else:
            failed_count += 1
    
    print("\nProcessing Summary:")
    print(f"  Chunks processed: {success_count}/{len(file_chunks)}")
    print(f"  Chunks failed: {failed_count}/{len(file_chunks)}")
    print(f"  Original parquet files: {len(parquet_files)}")
    print(f"  Total genomes processed: {total_genomes_processed}")
    print(f"  Total AMR hits found: {total_amr_hits}")
    if len(file_chunks) > 0:
        print(f"  Average processing time per chunk: {elapsed/len(file_chunks):.1f}s")
    
    # Perform comprehensive validation
    print("\n" + "="*80)
    print("AMR RESULTS VALIDATION")
    print("="*80)
    
    validation_results = validate_amr_results(all_expected_genomes, results_dir, all_successfully_processed_genomes)
    
    print(f"\nExpected vs Actual:")
    print(f"  Expected genomes total: {validation_results['total_expected_genomes']}")
    print(f"  Successfully processed: {validation_results['total_successfully_processed']}")
    print(f"  AMR result files found: {validation_results['total_amr_files_found']}")
    print(f"  Genomes with result files: {validation_results['genomes_with_files']}")
    
    # Report any issues
    if validation_results['missing_files']:
        print(f"\n⚠ MISSING AMR FILES ({len(validation_results['missing_files'])}):")
        for genome_name in sorted(list(validation_results['missing_files'])[:20]):  # Show first 20
            print(f"  - {genome_name}")
        if len(validation_results['missing_files']) > 20:
            print(f"  ... and {len(validation_results['missing_files']) - 20} more")
    
    if validation_results['successfully_processed_missing_files']:
        print(f"\n⚠ SUCCESSFULLY PROCESSED BUT NO AMR FILE ({len(validation_results['successfully_processed_missing_files'])}):")
        for genome_name in sorted(list(validation_results['successfully_processed_missing_files'])[:10]):
            print(f"  - {genome_name}")
        if len(validation_results['successfully_processed_missing_files']) > 10:
            print(f"  ... and {len(validation_results['successfully_processed_missing_files']) - 10} more")
    
    if validation_results['unexpected_files']:
        print(f"\n⚠ UNEXPECTED AMR FILES ({len(validation_results['unexpected_files'])}):")
        for genome_name in sorted(list(validation_results['unexpected_files'])[:10]):
            print(f"  - {genome_name}")
        if len(validation_results['unexpected_files']) > 10:
            print(f"  ... and {len(validation_results['unexpected_files']) - 10} more")
    
    # Final validation status
    print(f"\n" + "="*80)
    if validation_results['validation_success']:
        print("✅ VALIDATION PASSED - All expected genomes have AMR result files")
    else:
        print("❌ VALIDATION FAILED - Some genomes are missing AMR result files")
    
    print(f"\nFile Statistics by Parquet File:")
    for filename, stats in file_stats.items():
        print(f"  {filename}: {stats['genome_count']} genomes")
    
    # List sample output files
    output_files = sorted(results_dir.glob("*_amr_results.tsv"))
    print(f"\nSample AMR result files ({len(output_files)} total):")
    for f in sorted(output_files)[:10]:  # Show first 10
        size_kb = f.stat().st_size / 1024
        print(f"  {f.name} ({size_kb:.1f} KB)")
    if len(output_files) > 10:
        print(f"  ... and {len(output_files) - 10} more files")
    
    # List log files
    log_files = sorted(log_dir.glob("*_processing.log"))
    print(f"\nLog files created ({len(log_files)}):")
    for f in log_files[:10]:  # Show first 10
        print(f"  {f.name}")
    if len(log_files) > 10:
        print(f"  ... and {len(log_files) - 10} more files")
    
    print("="*80)
    
    # Exit with appropriate code
    if failed_count > 0:
        sys.exit(1)
    elif not validation_results['validation_success']:
        sys.exit(2)  # Validation failed
    else:
        sys.exit(0)  # Success

if __name__ == "__main__":
    main()

