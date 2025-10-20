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


def scan_already_processed_genomes(results_dir):
    """
    Scan the AMR results directory and extract genome names that have already been processed.
    
    Args:
        results_dir: Path to directory containing AMR result files (*_amr_results.tsv)
        
    Returns:
        Set of genome names that have already been processed
    """
    results_path = Path(results_dir)
    if not results_path.exists():
        print(f"Results directory does not exist yet: {results_path}")
        return set()
    
    # Find all AMR result files
    result_files = list(results_path.glob("*_amr_results.tsv"))
    
    # Extract genome names from filenames
    processed_genomes = set()
    for result_file in result_files:
        # Remove "_amr_results.tsv" suffix to get genome name
        genome_name = result_file.name.replace("_amr_results.tsv", "")
        processed_genomes.add(genome_name)
    
    print(f"Found {len(processed_genomes)} already-processed genomes in {results_path}")
    return processed_genomes


def create_file_chunks(parquet_files, chunks_per_file):
    """
    Read parquet files once and split into dataframe chunks for parallel processing.
    
    Args:
        parquet_files: List of parquet file paths
        chunks_per_file: Number of chunks to split each file into
        
    Returns:
        List of tuples: (parquet_file_path, df_chunk, chunk_id)
        where df_chunk is the pre-sliced dataframe for this chunk
    """
    import pandas as pd
    
    chunks = []
    print(f"\nReading and chunking {len(parquet_files)} files...")
    
    for pf in parquet_files:
        # Read parquet file ONCE
        try:
            df = pd.read_parquet(pf, engine='pyarrow')
            total_rows = len(df)
            
            if chunks_per_file == 1:
                # No chunking, use entire dataframe
                chunks.append((pf, df, None))
                print(f"  {pf.name}: {total_rows} genomes (no chunking)")
            else:
                # Split dataframe into chunks
                chunk_size = total_rows // chunks_per_file
                for i in range(chunks_per_file):
                    start_idx = i * chunk_size
                    # Last chunk gets any remaining rows
                    end_idx = (i + 1) * chunk_size if i < chunks_per_file - 1 else total_rows
                    
                    # Pre-slice the dataframe
                    df_chunk = df.iloc[start_idx:end_idx].copy()
                    chunks.append((pf, df_chunk, i + 1))
                
                print(f"  {pf.name}: {total_rows} genomes → {chunks_per_file} chunks of ~{chunk_size} genomes")
                
                # Delete the full dataframe to free memory
                del df
                
        except Exception as e:
            print(f"  Warning: Could not read {pf.name}: {e}")
            # Fall back to file-based processing (no pre-loaded chunk)
            chunks.append((pf, None, None))
    
    print(f"Created {len(chunks)} pre-loaded chunks\n")
    return chunks

def process_single_file_wrapper(args):
    """
    Wrapper function that imports and calls the processing function.
    This avoids pickling issues with multiprocessing.
    
    Args:
        args: Tuple of (parquet_file_path, log_file_path, skip_genomes, fasta_dir, output_dir, 
                        df_chunk, chunk_id)
        
    Returns:
        Tuple of (filename, success, message, stats)
    """
    parquet_file, log_file, skip_genomes, fasta_dir, output_dir, df_chunk, chunk_id = args
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
        
        # Call the function with pre-loaded chunk
        result_dir, stats = process_parquet_file(
            parquet_file, 
            verbose=False, 
            log_file=log_file,
            skip_genomes=skip_genomes,
            fasta_dir=fasta_dir,
            output_dir=output_dir,
            chunk_id=chunk_id,
            df_chunk=df_chunk  # Pass pre-loaded dataframe chunk
        )
        
        elapsed = (datetime.now() - start_time).total_seconds()
        
        message = (f"Completed in {elapsed:.1f}s - "
                  f"{stats['genomes_processed']} processed, "
                  f"{stats['genomes_skipped']} skipped, "
                  f"{stats['amr_hits']} AMR hits")
        
        return (file_display, True, message, stats)
            
    except Exception as e:
        elapsed = (datetime.now() - start_time).total_seconds()
        import traceback
        error_msg = f"Error after {elapsed:.1f}s: {str(e)}\n{traceback.format_exc()}"
        return (file_display, False, error_msg, None)

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
  python parallel_process_batches.py --input-dir /path/to/parquets --cpus 76 --chunks-per-file 5
  
  # Specify custom output directories (no chunking)
  python parallel_process_batches.py --input-dir /path/to/input \\
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
        '--no-skip',
        action='store_true',
        help='Do not skip already-processed genomes (process all)'
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
    
    # Scan for already-processed genomes (unless --no-skip specified)
    skip_genomes = set()
    if not args.no_skip:
        skip_genomes = scan_already_processed_genomes(results_dir)
    else:
        print("Skipping check for already-processed genomes (--no-skip flag set)")
    
    # Find all parquet files (recursively)
    parquet_files = sorted(input_dir.rglob("*.parquet"))
    
    if not parquet_files:
        print(f"Error: No .parquet files found in {input_dir}")
        sys.exit(1)
    
    # Create chunks from parquet files
    file_chunks = create_file_chunks(parquet_files, args.chunks_per_file)
    
    print("="*80)
    print("AMR FINDER PLUS - PARALLEL BATCH PROCESSING")
    print("="*80)
    print(f"\nInput directory: {input_dir}")
    print(f"FASTA directory: {fasta_dir}")
    print(f"Results directory: {results_dir}")
    print(f"Log directory: {log_dir}")
    print(f"Found {len(parquet_files)} parquet files")
    print(f"Chunks per file: {args.chunks_per_file}")
    print(f"Total work items: {len(file_chunks)}")
    print(f"Already-processed genomes to skip: {len(skip_genomes)}")
    print(f"Requested CPUs: {n_cpus}")
    print(f"Parallel workers that will be active: {min(n_cpus, len(file_chunks))}")
    print(f"\nStart time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("="*80)
    
    # Create work items for each chunk
    file_args = []
    for pf, df_chunk, chunk_id in file_chunks:
        # Create unique log file name
        if chunk_id is not None:
            log_file = log_dir / f"{pf.stem}_chunk{chunk_id}_processing.log"
        else:
            log_file = log_dir / f"{pf.stem}_processing.log"
        
        # Clear existing log file
        if log_file.exists():
            log_file.unlink()
        
        file_args.append((pf, log_file, skip_genomes, fasta_dir, results_dir, df_chunk, chunk_id))
    
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
                filename, success, message, stats = result
                status = "✓" if success else "✗"
                # Print detailed completion message
                print(f"\n{status} {filename}: {message}")
                pbar.set_postfix_str(f"{status} {filename}")
                pbar.update(1)
    
    # Print results summary
    elapsed = (datetime.now() - start_time).total_seconds()
    
    print("\n" + "="*80)
    print("PROCESSING COMPLETE")
    print("="*80)
    print(f"\nTotal time: {elapsed/60:.1f} minutes ({elapsed:.1f} seconds)")
    print("\nDetailed Results:")
    
    success_count = 0
    failed_count = 0
    total_genomes = 0
    total_amr_hits = 0
    total_skipped = 0
    
    for filename, success, message, stats in results:
        status = "✓" if success else "✗"
        print(f"  {status} {filename}: {message}")
        if success:
            success_count += 1
            if stats:
                total_genomes += stats['genomes_processed']
                total_amr_hits += stats['amr_hits']
                total_skipped += stats.get('genomes_skipped', 0)
        else:
            failed_count += 1
    
    print("\nOverall Summary:")
    print(f"  Chunks processed: {success_count}/{len(file_chunks)}")
    print(f"  Chunks failed: {failed_count}/{len(file_chunks)}")
    print(f"  Original parquet files: {len(parquet_files)}")
    print(f"  Total genomes processed (new): {total_genomes}")
    print(f"  Total genomes skipped (already done): {total_skipped}")
    print(f"  Total AMR hits found (new): {total_amr_hits}")
    if len(file_chunks) > 0:
        print(f"  Average processing time per chunk: {elapsed/len(file_chunks):.1f}s")
    
    # List output files
    output_files = sorted(results_dir.glob("*_amr_results.tsv"))
    print(f"\nTotal AMR result files in output directory: {len(output_files)}")
    if len(output_files) <= 50:
        for f in output_files:
            size_kb = f.stat().st_size / 1024
            print(f"  {f.name} ({size_kb:.1f} KB)")
    else:
        print("  (too many files to list individually)")
    
    # List log files
    log_files = sorted(log_dir.glob("*_processing.log"))
    print(f"\nLog files created ({len(log_files)}):")
    for f in log_files:
        print(f"  {f.name}")
    
    print("="*80)
    
    sys.exit(0 if failed_count == 0 else 1)

if __name__ == "__main__":
    main()

