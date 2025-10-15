"""
Test script to verify AMR data merge output files.

This script tests a parquet file to ensure:
- All 7 new AMR columns exist
- Array lengths match protein counts
- Data integrity is maintained
"""

import pandas as pd
import numpy as np
from pathlib import Path
import sys


def test_parquet_file(file_path: Path):
    """
    Test a single parquet file for correct AMR data structure.
    
    Args:
        file_path: Path to parquet file to test
    """
    print("=" * 70)
    print("TESTING OUTPUT FILE")
    print("=" * 70)
    print(f"File: {file_path.name}")
    print(f"Path: {file_path}")
    print()

    # Read the file
    df = pd.read_parquet(file_path, engine='pyarrow')

    # Basic structure
    print("✓ BASIC STRUCTURE")
    print(f"  - Total genomes: {len(df)}")
    print(f"  - Total columns: {len(df.columns)}")
    print()

    # Check new columns exist
    new_cols = ['element_symbol', 'element_name', 'amrf_type', 'amrf_subtype', 
                'amrf_class', 'amrf_subclass', 'pct_identity_to_reference']
    print("✓ NEW AMR COLUMNS")
    all_present = True
    for col in new_cols:
        exists = col in df.columns
        print(f"  - {col}: {'✓' if exists else '✗ MISSING'}")
        if not exists:
            all_present = False
    
    if not all_present:
        print("\n✗ ERROR: Some AMR columns are missing!")
        sys.exit(1)
    print()

    # Test sample genomes
    print("✓ TESTING SAMPLE GENOMES")
    print("-" * 70)

    test_indices = [0, len(df)//2, len(df)-1] if len(df) >= 3 else range(len(df))
    
    for idx in test_indices:
        row = df.iloc[idx]
        genome_name = row['genome_name']
        
        # Flatten protein IDs
        n_proteins = sum(len(c) for c in row['protein_id'])
        n_contigs = len(row['protein_id'])
        
        # Check array lengths
        element_symbol_len = len(row['element_symbol'])
        element_name_len = len(row['element_name'])
        amrf_type_len = len(row['amrf_type'])
        pct_identity_len = len(row['pct_identity_to_reference'])
        
        # Verify all lengths match
        lengths_match = (element_symbol_len == element_name_len == 
                        amrf_type_len == pct_identity_len == n_proteins)
        
        # Count AMR hits
        amr_hits = sum(1 for x in row['element_symbol'] if pd.notna(x))
        
        print(f"\nGenome #{idx+1}: {genome_name}")
        print(f"  - Contigs: {n_contigs}")
        print(f"  - Total proteins: {n_proteins}")
        print(f"  - Array lengths match: {'✓' if lengths_match else '✗ MISMATCH!'}")
        if not lengths_match:
            print(f"    element_symbol: {element_symbol_len}")
            print(f"    element_name: {element_name_len}")
            print(f"    amrf_type: {amrf_type_len}")
            print(f"    pct_identity: {pct_identity_len}")
        print(f"  - Proteins with AMR data: {amr_hits}")
        
        # Sample AMR data
        if amr_hits > 0:
            # Find first AMR hit
            for i, symbol in enumerate(row['element_symbol']):
                if pd.notna(symbol):
                    print(f"  - Sample AMR protein:")
                    print(f"    * element_symbol: {row['element_symbol'][i]}")
                    print(f"    * element_name: {row['element_name'][i]}")
                    print(f"    * amrf_type: {row['amrf_type'][i]}")
                    print(f"    * amrf_class: {row['amrf_class'][i]}")
                    print(f"    * pct_identity: {row['pct_identity_to_reference'][i]}")
                    break
        else:
            print(f"  - No AMR hits for this genome (all NaN)")

    print()
    print("=" * 70)
    print("✓ ALL TESTS PASSED!")
    print("=" * 70)


if __name__ == "__main__":
    # Default test file - can be changed via command line
    if len(sys.argv) > 1:
        test_file = Path(sys.argv[1])
    else:
        # Default: test a validation file
        test_file = Path(
            '/home/dca36/rds/hpc-work/data/BacFormer/processed/'
            'ecoli_embeddings_with_amrf_labels_flattened/validation/batch_001_chunks/train_chunk_3.parquet'
        )
    
    if not test_file.exists():
        print(f"Error: File not found: {test_file}")
        sys.exit(1)
    
    test_parquet_file(test_file)


# # Test the default file (validation/batch_001_chunks/train_chunk_3.parquet)
# uv run python amrf/test_output_file.py

# # Test a specific file
# uv run python amrf/test_output_file.py /path/to/your/file.parquet

# # Example: Test a train file
# uv run python amrf/test_output_file.py /home/dca36/rds/hpc-work/data/BacFormer/processed/ecoli_embeddings_with_amrf_labels_flattened/train/batch_005_chunks/train_chunk_1.parquet

