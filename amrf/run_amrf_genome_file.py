import pandas as pd
import subprocess
import os
import sys
import argparse
from pathlib import Path
from tqdm import tqdm

# Species translation dictionary: maps user species names (lowercase with underscores) 
# to AMR Finder Plus organism names
SPECIES_TRANSLATION = {
    'acinetobacter_baumannii': 'Acinetobacter_baumannii',
    'bordetella_pertussis': 'Bordetella_pertussis',
    'burkholderia_cepacia': 'Burkholderia_cepacia',
    'burkholderia_mallei': 'Burkholderia_mallei',
    'burkholderia_pseudomallei': 'Burkholderia_pseudomallei',
    'campylobacter_jejuni': 'Campylobacter',  # Special case
    'citrobacter_freundii': 'Citrobacter_freundii',
    'clostridioides_difficile': 'Clostridioides_difficile',
    'corynebacterium_diphtheriae': 'Corynebacterium_diphtheriae',
    'enterobacter_asburiae': 'Enterobacter_asburiae',
    'enterobacter_cloacae': 'Enterobacter_cloacae',
    'enterococcus_faecalis': 'Enterococcus_faecalis',
    'enterococcus_faecium': 'Enterococcus_faecium',
    'ecoli': 'Escherichia',  # Special case
    'haemophilus_influenzae': 'Haemophilus_influenzae',
    'klebsiella_oxytoca': 'Klebsiella_oxytoca',
    'klebsiella_pneumoniae': 'Klebsiella_pneumoniae',
    'neisseria_gonorrhoeae': 'Neisseria_gonorrhoeae',
    'neisseria_meningitidis': 'Neisseria_meningitidis',
    'pseudomonas_aeruginosa': 'Pseudomonas_aeruginosa',
    'salmonella': 'Salmonella',
    'serratia_marcescens': 'Serratia_marcescens',
    'staphylococcus_aureus': 'Staphylococcus_aureus',
    'staphylococcus_pseudintermedius': 'Staphylococcus_pseudintermedius',
    'streptococcus_agalactiae': 'Streptococcus_agalactiae',
    'streptococcus_pneumoniae': 'Streptococcus_pneumoniae',
    'streptococcus_pyogenes': 'Streptococcus_pyogenes',
    'vibrio_cholerae': 'Vibrio_cholerae',
    'vibrio_parahaemolyticus': 'Vibrio_parahaemolyticus',
    'vibrio_vulnificus': 'Vibrio_vulnificus',
}


def process_parquet_file(input_file, species=None, n_samples=None, verbose=True, log_file=None, 
                         fasta_dir=None, output_dir=None,
                         start_idx=None, end_idx=None, chunk_id=None, df_chunk=None,
                         target_genomes=None, debug_verbose=False):
    """
    Process a parquet file (or chunk) through AMR Finder Plus.
    
    Args:
        input_file (str or Path): Path to input parquet file (used for naming if df_chunk provided)
        species (str, required): Species name in lowercase with underscores (e.g., 'campylobacter_jejuni')
        n_samples (int, optional): Number of genomes to process. If None, process all.
        verbose (bool): Whether to print progress information
        log_file (str or Path, optional): Path to log file for output
        fasta_dir (Path, optional): Directory for FASTA files. Defaults to standard location.
        output_dir (Path, optional): Directory for AMR results. Defaults to standard location.
        start_idx (int, optional): Start index for processing chunk (row number in parquet)
        end_idx (int, optional): End index for processing chunk (exclusive, row number in parquet)
        chunk_id (int, optional): Chunk identifier for logging
        df_chunk (pd.DataFrame, optional): Pre-loaded dataframe chunk. If provided, skips file reading.
        target_genomes (set, optional): Set of specific genome names to process. If provided, only these genomes will be processed.
        debug_verbose (bool): Whether to print detailed per-genome debugging information
        
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
    
    # Validate and translate species name
    if species is None:
        raise ValueError("species parameter is required. Please provide a species name.")
    
    if species not in SPECIES_TRANSLATION:
        available_species = ', '.join(sorted(SPECIES_TRANSLATION.keys()))
        raise ValueError(
            f"Unknown species '{species}'. Available species are: {available_species}"
        )
    
    amr_finder_species = SPECIES_TRANSLATION[species]
    
    # Validate that species exists in AMR Finder
    log("\nChecking available AMRFinder organisms...")
    org_check = subprocess.run(["amrfinder", "--list_organisms"], capture_output=True, text=True)
    if org_check.returncode != 0:
        raise RuntimeError(
            f"Failed to check AMRFinder organisms. Command returned code {org_check.returncode}.\n"
            f"STDERR: {org_check.stderr}"
        )
    
    if amr_finder_species not in org_check.stdout:
        available_orgs = [line.strip() for line in org_check.stdout.split('\n') if line.strip()]
        available_orgs_str = '\n'.join(available_orgs[:20])  # Show first 20
        if len(available_orgs) > 20:
            available_orgs_str += f"\n... and {len(available_orgs) - 20} more"
        raise ValueError(
            f"AMR Finder organism '{amr_finder_species}' (translated from '{species}') "
            f"not found in AMRFinder.\n\n"
            f"Available organisms:\n{available_orgs_str}"
        )
    
    log(f"✓ Using organism: {amr_finder_species} (from species: {species})")
    
    # Setup directories - use provided paths or defaults
    if fasta_dir is None:
        fasta_dir = Path("/home/dca36/rds/hpc-work/data/amr_finder_plus/data/fasta_files")
    else:
        fasta_dir = Path(fasta_dir)
    
    if output_dir is None:
        output_dir = Path("/home/dca36/rds/hpc-work/data/amr_finder_plus/data/amrf_results")
    else:
        output_dir = Path(output_dir)
    
    # Create directories
    fasta_dir.mkdir(parents=True, exist_ok=True)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    
    input_file = Path(input_file)
    input_basename = input_file.stem  # e.g., "campylobacter_jejuni_batch_000"
    
    # Add chunk identifier to basename if processing a chunk
    if chunk_id is not None:
        input_basename = f"{input_basename}_chunk{chunk_id}"
    
    # Use pre-loaded chunk if provided, otherwise read from file
    chunk_info = f" (chunk {chunk_id})" if chunk_id is not None else ""
    
    if df_chunk is not None:
        # Use pre-loaded dataframe chunk
        log(f"Processing pre-loaded chunk from: {input_file}{chunk_info}")
        genomes_df = df_chunk
    else:
        # Read the parquet file
        log(f"Reading: {input_file}{chunk_info}")
        genomes_df = pd.read_parquet(input_file, engine='pyarrow')
        
        # Handle chunk processing (only needed if not pre-chunked)
        total_genomes = len(genomes_df)
        if start_idx is not None or end_idx is not None:
            start = start_idx if start_idx is not None else 0
            end = end_idx if end_idx is not None else len(genomes_df)
            log(f"Processing chunk: rows {start} to {end} of {total_genomes}")
            genomes_df = genomes_df.iloc[start:end]
    
    # Limit to n_samples if specified
    if n_samples is not None:
        genomes_df = genomes_df.head(n_samples)
        log(f"⚠ Limited to first {n_samples} genomes for testing")
    
    # Check for empty chunk and return early if no genomes to process
    if len(genomes_df) == 0:
        log("⚠ Empty chunk - no genomes to process")
        return output_dir, {
            'genomes_processed': 0,
            'unique_proteins': 0,
            'duplicate_proteins': 0,
            'amr_hits': 0
        }
    
    log(f"Processing {len(genomes_df)} genomes")
    log(f"Columns: {genomes_df.columns.tolist()}")
    
    # Calculate protein and contig statistics
    total_proteins = 0
    total_contigs = 0
    
    for idx, row in genomes_df.iterrows():
        protein_ids_by_contig = row['protein_id']
        num_contigs = len(protein_ids_by_contig)
        total_contigs += num_contigs
        
        for contig_proteins in protein_ids_by_contig:
            total_proteins += len(contig_proteins)
    
    log("\nDataset statistics:")
    log(f"  Total genomes: {len(genomes_df)}")
    log(f"  Total contigs: {total_contigs}")
    log(f"  Total proteins: {total_proteins}")
    log(f"  Average contigs per genome: {total_contigs/len(genomes_df):.2f}")
    log(f"  Average proteins per genome: {total_proteins/len(genomes_df):.2f}")
    log(f"  Average proteins per contig: {total_proteins/total_contigs:.2f}")
    
    # Global tracking
    global_seen_ids = set()
    global_duplicate_count = 0
    global_unique_count = 0
    total_amr_hits = 0
    genomes_with_hits = 0
    
    log("\n" + "="*80)
    log(f"Processing {len(genomes_df)} genomes...")
    if target_genomes:
        log(f"Will only process {len(target_genomes)} target genomes")
    log("="*80 + "\n")
    
    # Process each genome
    iterator = tqdm(range(len(genomes_df)), desc="Processing genomes") if verbose and not debug_verbose else range(len(genomes_df))
    
    for genome_idx in iterator:
        genome_row = genomes_df.iloc[genome_idx]
        genome_name = genome_row['genome_name']
        
        # Skip if not in target genomes (when target filtering is enabled)
        if target_genomes is not None and genome_name not in target_genomes:
            continue
        
        protein_sequences_by_contig = genome_row['protein_sequence']
        protein_ids_by_contig = genome_row['protein_id']
        
        # Count total proteins for this genome
        genome_total_proteins = sum(len(contig_proteins) for contig_proteins in protein_ids_by_contig)
        
        # Debug logging
        if debug_verbose:
            log(f"\n[Genome {genome_idx + 1}/{len(genomes_df)}] Processing: {genome_name}")
            log(f"  Total proteins: {genome_total_proteins}")
        
        # Write this genome's proteins to FASTA with unique name (includes batch name)
        fasta_file = fasta_dir / f"{input_basename}_{genome_name}_proteins.fasta"
        
        if debug_verbose:
            log(f"  FASTA file: {fasta_file}")
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
            if debug_verbose:
                log("  No unique proteins, skipping")
            os.remove(fasta_file)
            continue
        
        if debug_verbose:
            log(f"  Unique proteins: {genome_unique_count} (writing to FASTA)")
            log(f"  Duplicate proteins: {genome_duplicate_count} (skipped)")
        
        # Run AMR Finder Plus - save directly to final location
        output_file = output_dir / f"{genome_name}_amr_results.tsv"
        cmd = [
            "amrfinder",
            "--protein", str(fasta_file),
            "--organism", amr_finder_species,
            "--output", str(output_file),
            "--plus"
        ]
        
        if debug_verbose:
            log(f"  Output file: {output_file}")
            log(f"  Running AMR Finder Plus command: {' '.join(cmd)}")
        
        result = subprocess.run(cmd, capture_output=True, text=True)
        
        if debug_verbose:
            log(f"  AMR Finder Plus return code: {result.returncode}")
        
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
                
                num_hits = len(genome_amr_results)
                total_amr_hits += num_hits
                genomes_with_hits += 1
                
                if debug_verbose:
                    log(f"  AMR hits found: {num_hits}")
                    log("  Status: Complete ✓")
            else:
                # No hits - remove the empty file
                if debug_verbose:
                    log("  AMR hits found: 0")
                    log("  Status: No hits found")
                else:
                    log(f"No hits found for {genome_name}, removing {output_file}")
                os.remove(output_file)
        else:
            if debug_verbose:
                log(f"  ERROR: AMR Finder Plus failed with return code {result.returncode}")
                log(f"  STDERR: {result.stderr}")
                log(f"  STDOUT: {result.stdout}")
    
    log("\n" + "="*80)
    log("PROCESSING COMPLETE")
    log("="*80)
    log("\nResults summary:")
    log(f"  Genomes in file: {len(genomes_df)}")
    log(f"  Genomes processed: {len(genomes_df)}")
    log(f"  Genomes with AMR hits: {genomes_with_hits}")
    log(f"  Total unique proteins: {global_unique_count}")
    log(f"  Total duplicate proteins (skipped): {global_duplicate_count}")
    log(f"  Total AMR hits found: {total_amr_hits}")
    log(f"\nFASTA files saved to: {fasta_dir}")
    log(f"AMR results saved to: {output_dir}")
    log(f"  ({genomes_with_hits} individual genome result files)")
    
    # Return stats
    stats = {
        'genomes_processed': len(genomes_df),
        'unique_proteins': global_unique_count,
        'duplicate_proteins': global_duplicate_count,
        'amr_hits': total_amr_hits
    }
    
    return output_dir, stats


def main():
    """Main function for command line usage"""
    parser = argparse.ArgumentParser(
        description='Process parquet file through AMR Finder Plus',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Process all genomes in file
  python run_amrf_genome_file.py /path/to/file.parquet --species campylobacter_jejuni
  
  # Process first 5 genomes
  python run_amrf_genome_file.py /path/to/file.parquet --species campylobacter_jejuni --n-samples 5
        """
    )
    
    parser.add_argument(
        'input_file',
        type=str,
        help='Path to input parquet file'
    )
    parser.add_argument(
        '--species',
        type=str,
        required=True,
        choices=list(SPECIES_TRANSLATION.keys()),
        help='Species name (lowercase with underscores). Available: ' + ', '.join(sorted(SPECIES_TRANSLATION.keys()))
    )
    parser.add_argument(
        '--n-samples',
        type=int,
        default=None,
        help='Number of genomes to process (optional, default: all)'
    )
    parser.add_argument(
        '--fasta-dir',
        type=str,
        default=None,
        help='Directory for FASTA files (optional)'
    )
    parser.add_argument(
        '--output-dir',
        type=str,
        default=None,
        help='Directory for AMR results (optional)'
    )
    parser.add_argument(
        '--log-file',
        type=str,
        default=None,
        help='Path to log file (optional)'
    )
    
    args = parser.parse_args()
    
    # Run processing
    process_parquet_file(
        args.input_file,
        species=args.species,
        n_samples=args.n_samples,
        verbose=True,
        log_file=args.log_file,
        fasta_dir=args.fasta_dir,
        output_dir=args.output_dir
    )


if __name__ == "__main__":
    main()

