import pandas as pd
import subprocess
import os
from pathlib import Path
from tqdm import tqdm

# Setup directories
workspace = Path("/home/dca36/workspace/AMR_Finder_Plus")
temp_dir = workspace / "temp"
temp_dir.mkdir(exist_ok=True)

# Read the parquet file
ecoli_file = "/home/dca36/rds/hpc-work/data/BacFormer/processed/ecoli_batches/ecoli_batch_000.parquet"
ecoli_df = pd.read_parquet(ecoli_file, engine='pyarrow')

print(f"Loaded {len(ecoli_df)} genomes")
print(f"Columns: {ecoli_df.columns.tolist()}")

# Calculate protein and contig statistics
total_proteins = 0
total_contigs = 0

for idx, row in ecoli_df.iterrows():
    # protein_id and protein_sequence are lists of lists: [[contig1_proteins...], [contig2_proteins...], ...]
    protein_ids_by_contig = row['protein_id']
    
    # Count contigs (first level of list)
    num_contigs = len(protein_ids_by_contig)
    total_contigs += num_contigs
    
    # Count proteins (sum of all proteins across all contigs)
    for contig_proteins in protein_ids_by_contig:
        total_proteins += len(contig_proteins)

print(f"\nDataset statistics:")
print(f"  Total genomes: {len(ecoli_df)}")
print(f"  Total contigs: {total_contigs}")
print(f"  Total proteins: {total_proteins}")
print(f"  Average contigs per genome: {total_contigs/len(ecoli_df):.2f}")
print(f"  Average proteins per genome: {total_proteins/len(ecoli_df):.2f}")
print(f"  Average proteins per contig: {total_proteins/total_contigs:.2f}")

# Check available organisms
print("\nChecking available AMRFinder organisms...")
org_check = subprocess.run(["amrfinder", "--list_organisms"], capture_output=True, text=True)
if org_check.returncode == 0:
    if "Escherichia" in org_check.stdout:
        organism_name = "Escherichia"
        print(f"✓ Using organism: {organism_name}")
    else:
        print("⚠ 'Escherichia' not found, using default")
        organism_name = "Escherichia"
else:
    organism_name = "Escherichia"
    print(f"Using default organism: {organism_name}")

# Output file for all results
final_output_file = temp_dir / "all_genomes_amr_results.tsv"
all_results = []

# Global tracking
global_seen_ids = set()
global_duplicate_count = 0
global_unique_count = 0
total_amr_hits = 0

print(f"\n{'='*80}")
print(f"Processing {len(ecoli_df)} genomes...")
print(f"{'='*80}\n")

# Process each genome
for genome_idx in tqdm(range(len(ecoli_df)), desc="Processing genomes"):
    genome_row = ecoli_df.iloc[genome_idx]
    genome_name = genome_row['genome_name']
    protein_sequences_by_contig = genome_row['protein_sequence']
    protein_ids_by_contig = genome_row['protein_id']
    
    # Write this genome's proteins to FASTA
    fasta_file = temp_dir / f"genome_{genome_idx}_proteins.fasta"
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
    
    # Run AMR Finder Plus
    output_file = temp_dir / f"genome_{genome_idx}_amr.tsv"
    cmd = [
        "amrfinder",
        "--protein", str(fasta_file),
        "--organism", organism_name,
        "--output", str(output_file),
        "--plus"
    ]
    
    result = subprocess.run(cmd, capture_output=True, text=True)
    
    if result.returncode == 0 and output_file.exists():
        # Read results and add genome_name column
        genome_amr_results = pd.read_csv(output_file, sep='\t')
        if len(genome_amr_results) > 0:
            genome_amr_results['genome_name'] = genome_name
            all_results.append(genome_amr_results)
            total_amr_hits += len(genome_amr_results)
    
    # Clean up temporary files
    if fasta_file.exists():
        os.remove(fasta_file)
    if output_file.exists():
        os.remove(output_file)

# Combine all results
if all_results:
    combined_results = pd.concat(all_results, ignore_index=True)
    # Reorder columns to put genome_name first
    cols = combined_results.columns.tolist()
    cols.remove('genome_name')
    combined_results = combined_results[['genome_name'] + cols]
    
    # Save to file
    combined_results.to_csv(final_output_file, sep='\t', index=False)
    
    print(f"\n{'='*80}")
    print(f"PROCESSING COMPLETE")
    print(f"{'='*80}")
    print(f"\nResults summary:")
    print(f"  Genomes processed: {len(ecoli_df)}")
    print(f"  Total unique proteins: {global_unique_count}")
    print(f"  Total duplicate proteins (skipped): {global_duplicate_count}")
    print(f"  Total AMR hits found: {total_amr_hits}")
    print(f"\nResults saved to: {final_output_file}")
    print(f"\nFirst 10 hits preview:")
    print(combined_results.head(10).to_string(index=False))
else:
    print(f"\n{'='*80}")
    print(f"No AMR hits found across any genomes")
    print(f"{'='*80}")

