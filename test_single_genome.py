import pandas as pd
import subprocess
import os
from pathlib import Path

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

# Get first genome
first_genome = ecoli_df.iloc[0]
genome_name = first_genome['genome_name']
protein_sequences_by_contig = first_genome['protein_sequence']  # List of lists
protein_ids_by_contig = first_genome['protein_id']  # List of lists

# Count proteins in this genome
num_proteins_in_genome = sum(len(contig_proteins) for contig_proteins in protein_ids_by_contig)
num_contigs_in_genome = len(protein_ids_by_contig)

print(f"\nTesting with genome: {genome_name}")
print(f"  Contigs in genome: {num_contigs_in_genome}")
print(f"  Proteins in genome: {num_proteins_in_genome}")

# Write all proteins from all contigs to FASTA file (skip duplicates)
fasta_file = temp_dir / "genome_proteins.fasta"
seen_ids = set()
duplicate_count = 0
unique_count = 0

with open(fasta_file, 'w') as f:
    for contig_idx, (contig_sequences, contig_ids) in enumerate(zip(protein_sequences_by_contig, protein_ids_by_contig)):
        for protein_seq, protein_id in zip(contig_sequences, contig_ids):
            if protein_id in seen_ids:
                duplicate_count += 1
            else:
                seen_ids.add(protein_id)
                f.write(f">{protein_id}\n{protein_seq}\n")
                unique_count += 1

print(f"\nFASTA file statistics:")
print(f"  Total proteins in genome: {num_proteins_in_genome}")
print(f"  Unique protein IDs: {unique_count}")
print(f"  Duplicate protein IDs (skipped): {duplicate_count}")
print(f"  Wrote to: {fasta_file}")

# Check available organisms
print("\nChecking available AMRFinder organisms...")
org_check = subprocess.run(["amrfinder", "--list_organisms"], capture_output=True, text=True)
if org_check.returncode == 0:
    print("Available organisms:")
    print(org_check.stdout)
    
    # Look for Escherichia in the list
    if "Escherichia" in org_check.stdout:
        organism_name = "Escherichia"
        print(f"\n✓ Using organism: {organism_name}")
    else:
        print("\n⚠ 'Escherichia' not found in list - check output above for correct name")
        organism_name = "Escherichia"  # Default, but may need to change
else:
    organism_name = "Escherichia"
    print(f"Could not list organisms, using default: {organism_name}")

# Run AMR Finder Plus
output_file = temp_dir / "amr_output.tsv"
cmd = [
    "amrfinder",
    "--protein", str(fasta_file),
    "--organism", organism_name,
    "--output", str(output_file),
    "--plus"
]

print(f"\nRunning AMRFinder command: {' '.join(cmd)}")
result = subprocess.run(cmd, capture_output=True, text=True)

if result.returncode == 0:
    print("✓ AMR Finder Plus completed successfully!")
    
    # Read and display results
    if output_file.exists():
        amr_results = pd.read_csv(output_file, sep='\t')
        print(f"\nAMR Results: {len(amr_results)} hits found")
        print(f"Columns: {amr_results.columns.tolist()}")
        print("\n" + "="*100)
        print(amr_results.to_string(index=False))
        print("="*100)
    else:
        print("No AMR hits found for this genome")
else:
    print(f"✗ AMR Finder Plus failed:")
    print(f"STDERR: {result.stderr}")
    print(f"STDOUT: {result.stdout}")

print(f"\nFiles created in: {temp_dir}")

