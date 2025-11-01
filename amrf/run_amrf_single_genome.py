import pandas as pd
import subprocess
from pathlib import Path

# Setup directories
workspace = Path("/home/dca36/workspace/AMR_Finder_Plus")
temp_dir = workspace / "temp"
temp_dir.mkdir(exist_ok=True)

# Read the parquet file
parquet_file = "/home/dca36/rds/hpc-work/data/BacFormer/raw/ast/protein_sequences/campylobacter_jejuni/genomes_chunk_0.parquet"
df = pd.read_parquet(parquet_file, engine='pyarrow')

print(f"Loaded {len(df)} genomes")
print(f"Columns: {df.columns.tolist()}")

# Calculate protein and contig statistics
total_proteins = 0
total_contigs = 0

for idx, row in df.iterrows():
    # protein_id and protein_sequence are lists of lists: [[contig1_proteins...], [contig2_proteins...], ...]
    protein_ids_by_contig = row['protein_id']
    
    # Count contigs (first level of list)
    num_contigs = len(protein_ids_by_contig)
    total_contigs += num_contigs
    
    # Count proteins (sum of all proteins across all contigs)
    for contig_proteins in protein_ids_by_contig:
        total_proteins += len(contig_proteins)

print(f"\nDataset statistics:")
print(f"  Total genomes: {len(df)}")
print(f"  Total contigs: {total_contigs}")
print(f"  Total proteins: {total_proteins}")
print(f"  Average contigs per genome: {total_contigs/len(df):.2f}")
print(f"  Average proteins per genome: {total_proteins/len(df):.2f}")
print(f"  Average proteins per contig: {total_proteins/total_contigs:.2f}")

# Get first genome
first_genome = df.iloc[0]
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
    
    # Look for Campylobacter in the list
    organism_name = None
    if "Campylobacter" in org_check.stdout:
        #Possible organisms: Acinetobacter_baumannii, Bordetella_pertussis, Burkholderia_cepacia, Burkholderia_mallei, Burkholderia_pseudomallei, Campylobacter, Citrobacter_freundii, Clostridioides_difficile, Corynebacterium_diphtheriae, Enterobacter_asburiae, Enterobacter_cloacae, Enterococcus_faecalis, Enterococcus_faecium, Escherichia, Haemophilus_influenzae, Klebsiella_oxytoca, Klebsiella_pneumoniae, Neisseria_gonorrhoeae, Neisseria_meningitidis, Pseudomonas_aeruginosa, Salmonella, Serratia_marcescens, Staphylococcus_aureus, Staphylococcus_pseudintermedius, Streptococcus_agalactiae, Streptococcus_pneumoniae, Streptococcus_pyogenes, Vibrio_cholerae, Vibrio_parahaemolyticus, Vibrio_vulnificus
        organism_name = "Campylobacter"
        # Print the full output to see exact format
        print("\nSearching for Campylobacter in organism list...")
        for line in org_check.stdout.split('\n'):
            if 'Campylobacter' in line:
                print(f"Found: {line}")
        organism_name = "Campylobacter"
        print(f"\n✓ Using organism: {organism_name}")
    else:
        print("\n⚠ 'Campylobacter' not found in list - check output above for correct name")
        organism_name = None  # Will try without organism parameter
else:
    organism_name = None
    print(f"Could not list organisms, will try without organism parameter")

# Run AMR Finder Plus
output_file = temp_dir / "amr_output.tsv"

# Build command - try with organism first, fall back to without if needed
cmd = [
    "amrfinder",
    "--protein", str(fasta_file),
    "--output", str(output_file),
    "--plus"
]

if organism_name:
    cmd.extend(["--organism", organism_name])
else:
    print("\n⚠ Note: Running without --organism parameter")

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

