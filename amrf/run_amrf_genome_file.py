import pandas as pd
import subprocess
import os
import sys
from pathlib import Path
from tqdm import tqdm


def process_parquet_file(input_file, n_samples=None, verbose=True, log_file=None):
    """
    Process a parquet file through AMR Finder Plus.
    
    Args:
        input_file (str or Path): Path to input parquet file
        n_samples (int, optional): Number of genomes to process. If None, process all.
        verbose (bool): Whether to print progress information
        log_file (str or Path, optional): Path to log file for output
        
    Returns:
        tuple: (output_dir_path, stats_dict) where stats_dict contains:
            - genomes_processed
            - unique_proteins
            - duplicate_proteins
            - amr_hits
    """
    # Setup logging
    def log(message):
        """Print to console if verbose, always write to log file if provided"""
        if verbose:
            print(message)
        if log_file:
            with open(log_file, 'a') as f:
                f.write(message + '\n')
    
    # Setup directories
    output_dir = Path("/home/dca36/rds/hpc-work/data/amr_finder_plus/data")
    fasta_dir = output_dir / "fasta_files"
    output_dir = output_dir / "amrf_results"
    
    # Create directories
    output_dir.mkdir(exist_ok=True)
    fasta_dir.mkdir(exist_ok=True)
    output_dir.mkdir(exist_ok=True)
    
    input_file = Path(input_file)
    input_basename = input_file.stem  # e.g., "ecoli_batch_000"
    
    # Read the parquet file
    log(f"Reading: {input_file}")
    ecoli_df = pd.read_parquet(input_file, engine='pyarrow')
    
    # Limit to n_samples if specified
    if n_samples is not None:
        ecoli_df = ecoli_df.head(n_samples)
        log(f"⚠ Limited to first {n_samples} genomes for testing")
    
    log(f"Loaded {len(ecoli_df)} genomes")
    log(f"Columns: {ecoli_df.columns.tolist()}")
    
    # Calculate protein and contig statistics
    total_proteins = 0
    total_contigs = 0
    
    for idx, row in ecoli_df.iterrows():
        protein_ids_by_contig = row['protein_id']
        num_contigs = len(protein_ids_by_contig)
        total_contigs += num_contigs
        
        for contig_proteins in protein_ids_by_contig:
            total_proteins += len(contig_proteins)
    
    log("\nDataset statistics:")
    log(f"  Total genomes: {len(ecoli_df)}")
    log(f"  Total contigs: {total_contigs}")
    log(f"  Total proteins: {total_proteins}")
    log(f"  Average contigs per genome: {total_contigs/len(ecoli_df):.2f}")
    log(f"  Average proteins per genome: {total_proteins/len(ecoli_df):.2f}")
    log(f"  Average proteins per contig: {total_proteins/total_contigs:.2f}")
    
    # Check available organisms
    organism_name = "Escherichia"
    log("\nChecking available AMRFinder organisms...")
    org_check = subprocess.run(["amrfinder", "--list_organisms"], capture_output=True, text=True)
    if org_check.returncode == 0:
        if "Escherichia" in org_check.stdout:
            log(f"✓ Using organism: {organism_name}")
        else:
            log("⚠ 'Escherichia' not found, using default")
    else:
        log(f"Using default organism: {organism_name}")
    
    # Global tracking
    global_seen_ids = set()
    global_duplicate_count = 0
    global_unique_count = 0
    total_amr_hits = 0
    genomes_with_hits = 0
    
    log("\n" + "="*80)
    log(f"Processing {len(ecoli_df)} genomes...")
    log("="*80 + "\n")
    
    # Process each genome
    iterator = tqdm(range(len(ecoli_df)), desc="Processing genomes") if verbose else range(len(ecoli_df))
    
    for genome_idx in iterator:
        genome_row = ecoli_df.iloc[genome_idx]
        genome_name = genome_row['genome_name']
        protein_sequences_by_contig = genome_row['protein_sequence']
        protein_ids_by_contig = genome_row['protein_id']
        
        # Write this genome's proteins to FASTA with unique name (includes batch name)
        fasta_file = fasta_dir / f"{input_basename}_{genome_name}_proteins.fasta"
        genome_unique_count = 0
        genome_duplicate_count = 0
        
        with open(fasta_file, 'w') as f:
            for contig_sequences, contig_ids in zip(protein_sequences_by_contig, protein_ids_by_contig):
                for protein_seq, protein_id in zip(contig_sequences, contig_ids):
                    if protein_id in global_seen_ids:
                        genome_duplicate_count += 1
                        global_duplicate_count += 1
                    else:
                        global_seen_ids.add(protein_id)
                        f.write(f">{protein_id}\n{protein_seq}\n")
                        genome_unique_count += 1
                        global_unique_count += 1
        
        # Skip if no unique proteins
        if genome_unique_count == 0:
            os.remove(fasta_file)
            continue
        
        # Run AMR Finder Plus - save directly to final location
        output_file = output_dir / f"{genome_name}_amr_results.tsv"
        cmd = [
            "amrfinder",
            "--protein", str(fasta_file),
            "--organism", organism_name,
            "--output", str(output_file),
            "--plus"
        ]
        
        result = subprocess.run(cmd, capture_output=True, text=True)
        
        if result.returncode == 0 and output_file.exists():
            # Read results to check if there are any hits and clean up columns
            genome_amr_results = pd.read_csv(output_file, sep='\t')
            if len(genome_amr_results) > 0:
                # Rename columns to snake_case to match our convention
                column_mapping = {
                    'Protein id': 'protein_id',
                    'Element symbol': 'element_symbol',
                    'Element name': 'element_name',
                    'Closest reference name': 'closest_reference_name',
                    'HMM description': 'HMM_description'
                }
                genome_amr_results.rename(columns=column_mapping, inplace=True)
                
                # Save without genome_name column (it's in the filename)
                genome_amr_results.to_csv(output_file, sep='\t', index=False)
                
                total_amr_hits += len(genome_amr_results)
                genomes_with_hits += 1
            else:
                # No hits - remove the empty file
                log(f"No hits found for {genome_name}, removing {output_file}")
                os.remove(output_file)
    
    log("\n" + "="*80)
    log("PROCESSING COMPLETE")
    log("="*80)
    log("\nResults summary:")
    log(f"  Genomes processed: {len(ecoli_df)}")
    log(f"  Genomes with AMR hits: {genomes_with_hits}")
    log(f"  Total unique proteins: {global_unique_count}")
    log(f"  Total duplicate proteins (skipped): {global_duplicate_count}")
    log(f"  Total AMR hits found: {total_amr_hits}")
    log(f"\nFASTA files saved to: {fasta_dir}")
    log(f"AMR results saved to: {output_dir}")
    log(f"  ({genomes_with_hits} individual genome result files)")
    
    # Return stats
    stats = {
        'genomes_processed': len(ecoli_df),
        'unique_proteins': global_unique_count,
        'duplicate_proteins': global_duplicate_count,
        'amr_hits': total_amr_hits
    }
    
    return output_dir, stats


def main():
    """Main function for command line usage"""
    # Parse command line arguments
    n_samples = None
    ecoli_file = None
    
    if len(sys.argv) > 1:
        if sys.argv[1] in ['-h', '--help']:
            print("Usage: python run_amrf_genome_file.py [input_file] [n_samples]")
            print("\n  input_file: Path to input parquet file (optional, default: ecoli_batch_000.parquet)")
            print("  n_samples: Number of genomes to process (optional, default: all)")
            print("\nExamples:")
            print("  python run_amrf_genome_file.py                                    # Process all genomes in default file")
            print("  python run_amrf_genome_file.py /path/to/file.parquet              # Process all genomes in specified file")
            print("  python run_amrf_genome_file.py /path/to/file.parquet 5            # Process first 5 genomes")
            sys.exit(0)
        
        # First argument could be file path or n_samples
        try:
            n_samples = int(sys.argv[1])
            ecoli_file = "/home/dca36/rds/hpc-work/data/BacFormer/processed/ecoli_batches/ecoli_batch_000.parquet"
        except ValueError:
            # It's a file path
            ecoli_file = sys.argv[1]
            if len(sys.argv) > 2:
                try:
                    n_samples = int(sys.argv[2])
                except ValueError:
                    print(f"Error: n_samples must be an integer, got '{sys.argv[2]}'")
                    print("Use -h or --help for usage information")
                    sys.exit(1)
    else:
        ecoli_file = "/home/dca36/rds/hpc-work/data/BacFormer/processed/ecoli_batches/ecoli_batch_000.parquet"
    
    # Run processing
    process_parquet_file(ecoli_file, n_samples=n_samples, verbose=True)


if __name__ == "__main__":
    main()

