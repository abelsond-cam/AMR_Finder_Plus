import pandas as pd

# Load the parquet file with sequences
df = pd.read_parquet('/home/dca36/rds/hpc-work/data/BacFormer/processed/ecoli_batches/ecoli_batch_000.parquet')

print('Columns:', df.columns.tolist())
print('Shape:', df.shape)
print('\nFirst row data:')
print(f'genome_name: {df.iloc[0]["genome_name"]}')
print(f'protein_id: {df.iloc[0]["protein_id"]}')
print(f'product: {df.iloc[0]["product"]}')
print(f'protein_sequence (first 100 chars): {df.iloc[0]["protein_sequence"][:100]}...')
print(f'protein_sequence length: {len(df.iloc[0]["protein_sequence"])}')

# Check if sequences look like valid protein sequences
seq = df.iloc[0]["protein_sequence"]
print(f'\nSequence type: {type(seq)}')
if hasattr(seq, '__len__') and len(seq) > 0:
    if isinstance(seq, str):
        print(f'Sequence characters present: {set(seq[:50])}')
    else:
        print(f'First element of sequence array: {seq[0] if len(seq) > 0 else "empty"}')
        if len(seq) > 0 and isinstance(seq[0], str):
            print(f'Characters in first element: {set(seq[0][:50])}')
