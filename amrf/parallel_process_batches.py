#!/usr/bin/env python
"""
Parallel processing script for AMR Finder Plus on HPC
Fixed version with better multiprocessing handling
"""

import sys
from pathlib import Path
from multiprocessing import Pool, cpu_count
from datetime import datetime
from tqdm import tqdm

# Add workspace to path so we can import amrf module
WORKSPACE = Path("/home/dca36/workspace/AMR_Finder_Plus")
if str(WORKSPACE) not in sys.path:
    sys.path.insert(0, str(WORKSPACE))

def process_single_file_wrapper(args):
    """
    Wrapper function that imports and calls the processing function.
    This avoids pickling issues with multiprocessing.
    
    Args:
        args: Tuple of (parquet_file_path, log_file_path)
        
    Returns:
        Tuple of (filename, success, message, stats)
    """
    parquet_file, log_file = args
    start_time = datetime.now()
    
    try:
        # Add workspace to path in worker process
        import sys
        from pathlib import Path
        workspace = Path("/home/dca36/workspace/AMR_Finder_Plus")
        if str(workspace) not in sys.path:
            sys.path.insert(0, str(workspace))
        
        # Import here to avoid pickling issues
        from amrf.run_amrf_genome_file import process_parquet_file
        
        # Call the function with log file
        output_dir, stats = process_parquet_file(parquet_file, verbose=False, log_file=log_file)
        
        elapsed = (datetime.now() - start_time).total_seconds()
        
        message = (f"Completed in {elapsed:.1f}s - "
                  f"{stats['genomes_processed']} genomes, "
                  f"{stats['amr_hits']} AMR hits")
        
        return (Path(parquet_file).name, True, message, stats)
            
    except Exception as e:
        elapsed = (datetime.now() - start_time).total_seconds()
        import traceback
        error_msg = f"Error after {elapsed:.1f}s: {str(e)}\n{traceback.format_exc()}"
        return (Path(parquet_file).name, False, error_msg, None)

def main():
    # Configuration
    input_dir = Path("/home/dca36/rds/hpc-work/data/BacFormer/processed/ecoli_batches/")
    #workspace = Path("/home/dca36/workspace/AMR_Finder_Plus")
    data_dir = Path("/home/dca36/rds/hpc-work/data/amr_finder_plus/data")
    output_dir = data_dir / "amrf_results"
    log_dir = data_dir / "logs"
    
    # Create directories
    log_dir.mkdir(parents=True, exist_ok=True)
    
    # Parse command line arguments
    n_cpus = cpu_count()
    
    if len(sys.argv) > 1:
        if sys.argv[1] in ['-h', '--help']:
            print("Usage: python parallel_process_batches_fixed.py [n_cpus]")
            print(f"\n  n_cpus: Number of CPUs to use (optional, default: {cpu_count()})")
            print(f"\nInput directory: {input_dir}")
            print(f"Output directory: {output_dir}")
            sys.exit(0)
        try:
            n_cpus = int(sys.argv[1])
        except ValueError:
            print(f"Error: n_cpus must be an integer, got '{sys.argv[1]}'")
            sys.exit(1)
    
    # Find all parquet files
    parquet_files = sorted(input_dir.glob("*.parquet"))
    
    if not parquet_files:
        print(f"Error: No .parquet files found in {input_dir}")
        sys.exit(1)
    
    print("="*80)
    print("AMR FINDER PLUS - PARALLEL BATCH PROCESSING (FIXED)")
    print("="*80)
    print(f"\nInput directory: {input_dir}")
    print(f"Output directory: {output_dir}")
    print(f"Found {len(parquet_files)} parquet files to process")
    print(f"Requested CPUs: {n_cpus}")
    print(f"Files to process: {len(parquet_files)}")
    print(f"Parallel workers that will be active: {min(n_cpus, len(parquet_files))}")
    print(f"\nStart time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("="*80)
    
    # Create log file paths for each parquet file
    file_args = []
    for pf in parquet_files:
        log_file = log_dir / f"{pf.stem}_processing.log"
        # Clear existing log file
        if log_file.exists():
            log_file.unlink()
        file_args.append((pf, log_file))
    
    # Process files in parallel
    start_time = datetime.now()
    results = []
    
    print(f"\nSpawning {n_cpus} worker processes...")
    
    with Pool(processes=n_cpus) as pool:
        print(f"Pool created with {n_cpus} workers")
        print("Starting parallel processing...\n")
        
        # Use imap_unordered for immediate results
        with tqdm(total=len(parquet_files), desc="Processing files", unit="file") as pbar:
            for result in pool.imap_unordered(process_single_file_wrapper, file_args, chunksize=1):
                results.append(result)
                filename, success, message, _ = result
                status = "✓" if success else "✗"
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
    
    for filename, success, message, stats in results:
        status = "✓" if success else "✗"
        print(f"  {status} {filename}: {message}")
        if success:
            success_count += 1
            if stats:
                total_genomes += stats['genomes_processed']
                total_amr_hits += stats['amr_hits']
        else:
            failed_count += 1
    
    print("\nOverall Summary:")
    print(f"  Files processed: {success_count}/{len(parquet_files)}")
    print(f"  Files failed: {failed_count}/{len(parquet_files)}")
    print(f"  Total genomes processed: {total_genomes}")
    print(f"  Total AMR hits found: {total_amr_hits}")
    if len(parquet_files) > 0:
        print(f"  Average processing time per file: {elapsed/len(parquet_files):.1f}s")
    
    # List output files
    output_files = sorted(output_dir.glob("*_amr_results.tsv"))
    print(f"\nOutput files created ({len(output_files)}):")
    if len(output_files) <= 50:
        for f in output_files:
            size_kb = f.stat().st_size / 1024
            print(f"  {f.name} ({size_kb:.1f} KB)")
    else:
        print(f"  {len(output_files)} genome result files (too many to list individually)")
    
    # List log files
    log_files = sorted(log_dir.glob("*_processing.log"))
    print(f"\nLog files created ({len(log_files)}):")
    for f in log_files:
        print(f"  {f.name}")
    
    print("="*80)
    
    sys.exit(0 if failed_count == 0 else 1)

if __name__ == "__main__":
    main()

