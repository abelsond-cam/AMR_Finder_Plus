import pandas as pd

# Load the parquet file
df = pd.read_parquet('/home/dca36/rds/hpc-work/data/BacFormer/processed/ecoli_embeddings_with_ast_labels_flattened/train/batch_000_chunks/train_chunk_1.parquet')

print('Columns:', df.columns.tolist())
print('Shape:', df.shape)
print('\nFirst row sample:')
for col in df.columns[:10]:
    val = df.iloc[0][col] if pd.notna(df.iloc[0][col]) else "NaN"
    print(f'{col}: {str(val)[:100]}...' if len(str(val)) > 100 else f'{col}: {val}')

# Look for sequence-related columns
sequence_cols = [col for col in df.columns if 'seq' in col.lower() or 'protein' in col.lower() or 'gene' in col.lower()]
print(f'\nSequence-related columns: {sequence_cols}')
