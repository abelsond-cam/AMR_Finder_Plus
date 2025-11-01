#!/usr/bin/env python3
"""
Script to check parquet file for corruption issues.
"""
import os
import sys

try:
    import pandas as pd
    HAS_PANDAS = True
except ImportError:
    HAS_PANDAS = False

try:
    import pyarrow.parquet as pq
    import pyarrow as pa
    HAS_PYARROW = True
except ImportError:
    HAS_PYARROW = False

def check_parquet_file(filepath):
    print(f"Checking file: {filepath}")
    print(f"File exists: {os.path.exists(filepath)}")
    
    if os.path.exists(filepath):
        size = os.path.getsize(filepath)
        print(f"File size: {size} bytes")
    
    # Try reading with pyarrow first
    if HAS_PYARROW:
        print("\n--- Trying PyArrow ---")
        try:
            parquet_file = pq.ParquetFile(filepath)
            print(f"Schema: {parquet_file.schema}")
            print(f"Number of row groups: {parquet_file.num_row_groups}")
            
            metadata = parquet_file.metadata
            print(f"Number of rows: {metadata.num_rows}")
            print(f"Number of columns: {len(parquet_file.schema)}")
            
            # Try to read first row group
            print("\nFirst row group schema:")
            print(parquet_file.read_row_group(0).schema)
            
            # Try to read as table
            table = parquet_file.read()
            print(f"\nRead table successfully with shape: {table.shape}")
            print(f"Columns: {table.column_names}")
            
        except Exception as e:
            print(f"ERROR with PyArrow: {type(e).__name__}: {e}")
    
    # Try reading with pandas
    if HAS_PANDAS:
        print("\n--- Trying Pandas ---")
        try:
            df = pd.read_parquet(filepath)
            print(f"Read successfully with shape: {df.shape}")
            print(f"Columns: {df.columns.tolist()}")
            print(f"Memory usage: {df.memory_usage(deep=True).sum() / 1024**2:.2f} MB")
            
            # Show first few rows
            print("\nFirst 3 rows:")
            print(df.head(3))
            
        except Exception as e:
            print(f"ERROR with Pandas: {type(e).__name__}: {e}")

if __name__ == "__main__":
    filepath = "/home/dca36/rds/hpc-work/data/BacFormer/processed/ast_esm_embeddings/campylobacter_jejuni/evaluate/genomes_chunk_0_chunks/train_chunk_12.parquet"
    check_parquet_file(filepath)

