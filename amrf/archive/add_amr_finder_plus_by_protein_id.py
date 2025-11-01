#!/usr/bin/env python3
"""
Add AMR Finder Plus data to embeddings parquet file by matching protein IDs.

This script enriches the embeddings parquet file with AMR Finder Plus results
by matching protein IDs between the parquet file and corresponding TSV files.
"""

import pandas as pd
import numpy as np
from pathlib import Path
from tqdm import tqdm
import sys


def main():
    # File paths
    receiver_file = Path("/home/dca36/rds/hpc-work/data/BacFormer/models/captum/cluster_embeddings_by_attribution/drug_resistant_high_attribution_embeddings.parquet")
    amrf_dir = Path("/home/dca36/rds/hpc-work/data/amr_finder_plus/data/amrf_results")
    
    print(f"Loading parquet file: {receiver_file}")
    
    # Load the main parquet file
    try:
        df_embeddings = pd.read_parquet(receiver_file)
        print(f"Loaded parquet with {len(df_embeddings)} rows and {len(df_embeddings.columns)} columns")
    except Exception as e:
        print(f"Error loading parquet file: {e}")
        sys.exit(1)
    
    # Check if required columns exist
    if 'genome_name' not in df_embeddings.columns:
        print("Error: 'genome_name' column not found in parquet file")
        sys.exit(1)
    if 'protein_id' not in df_embeddings.columns:
        print("Error: 'protein_id' column not found in parquet file")
        sys.exit(1)
    
    # Get unique genome names
    unique_genomes = df_embeddings['genome_name'].unique()
    print(f"Found {len(unique_genomes)} unique genome names to process")
    
    # Initialize tracking variables
    missing_files = []
    missing_protein_ids = []
    successful_matches = 0
    total_embeddings_rows = len(df_embeddings)
    
    # Initialize new columns for AMR data
    amr_columns = ['element_symbol', 'element_name', 'type', 'subtype', 'class', 'subclass']
    for col in amr_columns:
        df_embeddings[col] = np.nan
    
    # Process each genome
    print("Processing genomes...")
    for genome_name in tqdm(unique_genomes, desc="Processing genomes"):
        # Construct TSV file path
        tsv_file = amrf_dir / f"{genome_name}_amr_results.tsv"
        
        if not tsv_file.exists():
            missing_files.append(genome_name)
            continue
        
        try:
            # Load the TSV file
            df_amr = pd.read_csv(tsv_file, sep='\t')
            
            # Check if required columns exist in TSV
            required_tsv_columns = ['protein_id', 'element_symbol', 'element_name', 'Type', 'Subtype', 'Class', 'Subclass']
            missing_cols = [col for col in required_tsv_columns if col not in df_amr.columns]
            if missing_cols:
                print(f"Warning: Missing columns in {tsv_file}: {missing_cols}")
                continue
            
            # Select and rename columns
            df_amr_selected = df_amr[required_tsv_columns].copy()
            df_amr_selected = df_amr_selected.drop_duplicates(subset=['protein_id'], keep='first')  # Handle duplicates
            df_amr_selected = df_amr_selected.rename(columns={
                'Type': 'type',
                'Subtype': 'subtype', 
                'Class': 'class',
                'Subclass': 'subclass'
            })
            
            # Get protein_ids for this genome from embeddings
            genome_mask = df_embeddings['genome_name'] == genome_name
            genome_protein_ids = df_embeddings.loc[genome_mask, 'protein_id'].unique()
            
            # Find matches and missing protein_ids
            matched_protein_ids = set(df_amr_selected['protein_id']).intersection(set(genome_protein_ids))
            unmatched_protein_ids = set(genome_protein_ids) - matched_protein_ids
            
            # Track missing protein_ids
            for pid in unmatched_protein_ids:
                missing_protein_ids.append((genome_name, pid))
            
            # Merge AMR data for matched protein_ids
            for protein_id in matched_protein_ids:
                # Get AMR data for this protein_id
                amr_data = df_amr_selected[df_amr_selected['protein_id'] == protein_id].iloc[0]
                
                # Update embeddings dataframe
                mask = (df_embeddings['genome_name'] == genome_name) & (df_embeddings['protein_id'] == protein_id)
                for col in amr_columns:
                    df_embeddings.loc[mask, col] = amr_data[col]
                
                successful_matches += mask.sum()
            
        except Exception as e:
            print(f"Error processing {tsv_file}: {e}")
            continue
    
    # Print summary statistics
    print("\n" + "="*60)
    print("PROCESSING SUMMARY")
    print("="*60)
    print(f"Total embeddings rows: {total_embeddings_rows}")
    print(f"Unique genomes processed: {len(unique_genomes)}")
    print(f"Missing TSV files: {len(missing_files)}")
    print(f"Successful protein_id matches: {successful_matches}")
    print(f"Missing protein_id matches: {len(missing_protein_ids)}")
    print(f"Match rate: {successful_matches/total_embeddings_rows*100:.2f}%")
    
    # Print missing files
    if missing_files:
        print(f"\nMissing TSV files for genomes:")
        for genome in missing_files[:10]:  # Show first 10
            print(f"  - {genome}")
        if len(missing_files) > 10:
            print(f"  ... and {len(missing_files) - 10} more")
    
    # Print some missing protein_ids
    if missing_protein_ids:
        print(f"\nSample missing protein_ids:")
        for genome, protein_id in missing_protein_ids[:10]:  # Show first 10
            print(f"  - Genome: {genome}, Protein ID: {protein_id}")
        if len(missing_protein_ids) > 10:
            print(f"  ... and {len(missing_protein_ids) - 10} more")
    
    # Save the enriched parquet file
    print(f"\nSaving enriched data back to: {receiver_file}")
    try:
        df_embeddings.to_parquet(receiver_file, index=False)
        print("Successfully saved enriched parquet file!")
    except Exception as e:
        print(f"Error saving parquet file: {e}")
        sys.exit(1)
    
    print("\nScript completed successfully!")


if __name__ == "__main__":
    main()
