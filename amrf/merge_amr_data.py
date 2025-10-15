"""
Core module for merging AMR Finder Plus results into parquet files.

This module processes individual parquet files containing E. coli genome embeddings
and merges AMR data from corresponding TSV files.
"""

import pandas as pd
import numpy as np
from pathlib import Path
import logging
from typing import List, Dict, Optional

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


def flatten_protein_ids(protein_id_nested: np.ndarray, genome_name: str = "") -> List[str]:
    """
    Flatten nested protein_id structure from list of arrays to single list.
    
    Args:
        protein_id_nested: Nested structure like [array(['prot1', 'prot2']), array(['prot3', 'prot4'])]
        genome_name: Optional genome name for debugging
    
    Returns:
        Flat list of all protein IDs
    """
    flat_list = []
    n_contigs = len(protein_id_nested)
    
    for contig_proteins in protein_id_nested:
        flat_list.extend(contig_proteins.tolist())
    
    n_proteins = len(flat_list)
    n_unique_proteins = len(set(flat_list))
    
    if genome_name:
        logger.info(f"  {genome_name}:")
        logger.info(f"    - Number of contigs: {n_contigs}")
        logger.info(f"    - Total proteins (after flattening): {n_proteins}")
        logger.info(f"    - Unique protein IDs: {n_unique_proteins}")
        if n_unique_proteins < n_proteins:
            logger.warning(f"    - WARNING: {n_proteins - n_unique_proteins} duplicate protein IDs found!")
    
    return flat_list


def load_amr_data(genome_name: str, amrf_results_dir: Path) -> Optional[pd.DataFrame]:
    """
    Load AMR results TSV file for a given genome.
    
    Args:
        genome_name: Name of the genome (e.g., 'GCA_000474825.1')
        amrf_results_dir: Directory containing AMR result files
    
    Returns:
        DataFrame with AMR results or None if file doesn't exist
    """
    amr_file_path = amrf_results_dir / f"{genome_name}_amr_results.tsv"
    
    if not amr_file_path.exists():
        logger.warning(f"AMR file not found for {genome_name}: {amr_file_path}")
        return None
    
    try:
        amr_df = pd.read_csv(amr_file_path, sep='\t')
        logger.debug(f"Loaded AMR data for {genome_name}: {len(amr_df)} proteins")
        return amr_df
    except Exception as e:
        logger.error(f"Error loading AMR file {amr_file_path}: {e}")
        return None


def map_amr_to_proteins(
    protein_ids_flat: List[str],
    amr_df: Optional[pd.DataFrame],
    genome_name: str = ""
) -> Dict[str, np.ndarray]:
    """
    Map AMR data to flattened protein list.
    
    Args:
        protein_ids_flat: Flat list of protein IDs
        amr_df: DataFrame with AMR results (or None if no AMR data)
        genome_name: Optional genome name for debugging
    
    Returns:
        Dictionary with 7 arrays (one for each AMR column)
    """
    n_proteins = len(protein_ids_flat)
    
    # Initialize arrays with NaN
    amr_data = {
        'element_symbol': np.full(n_proteins, np.nan, dtype=object),
        'element_name': np.full(n_proteins, np.nan, dtype=object),
        'amrf_type': np.full(n_proteins, np.nan, dtype=object),
        'amrf_subtype': np.full(n_proteins, np.nan, dtype=object),
        'amrf_class': np.full(n_proteins, np.nan, dtype=object),
        'amrf_subclass': np.full(n_proteins, np.nan, dtype=object),
        'pct_identity_to_reference': np.full(n_proteins, np.nan, dtype=float)
    }
    
    if amr_df is None or len(amr_df) == 0:
        logger.debug(f"No AMR data to map for {n_proteins} proteins")
        return amr_data
    
    # Debug AMR file info
    n_amr_rows = len(amr_df)
    n_unique_amr_proteins = len(amr_df['protein_id'].unique())
    
    if genome_name:
        logger.info(f"    - AMR file rows: {n_amr_rows}")
        logger.info(f"    - Unique protein IDs in AMR file: {n_unique_amr_proteins}")
        if n_unique_amr_proteins < n_amr_rows:
            logger.warning(f"    - WARNING: {n_amr_rows - n_unique_amr_proteins} duplicate protein IDs in AMR file!")
    
    # Create lookup dictionary from AMR data
    amr_lookup = {}
    for _, row in amr_df.iterrows():
        protein_id = row['protein_id']
        amr_lookup[protein_id] = {
            'element_symbol': row.get('element_symbol', np.nan),
            'element_name': row.get('element_name', np.nan),
            'amrf_type': row.get('Type', np.nan),
            'amrf_subtype': row.get('Subtype', np.nan),
            'amrf_class': row.get('Class', np.nan),
            'amrf_subclass': row.get('Subclass', np.nan),
            'pct_identity_to_reference': row.get('% Identity to reference', np.nan)
        }
    
    # Map AMR data to protein positions
    matches = 0
    for i, protein_id in enumerate(protein_ids_flat):
        if protein_id in amr_lookup:
            matches += 1
            for key in amr_data.keys():
                amr_data[key][i] = amr_lookup[protein_id][key]
    
    if genome_name:
        logger.info(f"    - Matched proteins from AMR file: {matches}")
    
    return amr_data


def validate_amr_data(
    protein_ids_flat: List[str],
    amr_data: Dict[str, np.ndarray],
    amr_df: Optional[pd.DataFrame],
    genome_name: str
) -> bool:
    """
    Validate that AMR data was correctly mapped.
    
    Args:
        protein_ids_flat: Flat list of protein IDs
        amr_data: Dictionary of AMR arrays
        amr_df: Original AMR DataFrame (or None)
        genome_name: Genome name for logging
    
    Returns:
        True if validation passes, False otherwise
    """
    n_proteins = len(protein_ids_flat)
    
    # Check 1: All arrays have same length as protein list
    for col_name, arr in amr_data.items():
        if len(arr) != n_proteins:
            logger.error(
                f"Validation failed for {genome_name}: "
                f"{col_name} has length {len(arr)}, expected {n_proteins}"
            )
            return False
    
    # Check 2: Count of non-NA values matches unique protein IDs in AMR file
    # Note: AMR files may have duplicate protein IDs (same protein ID with multiple hits)
    non_na_count = np.sum(~pd.isna(amr_data['element_symbol']))
    expected_count = len(amr_df['protein_id'].unique()) if amr_df is not None else 0
    
    if non_na_count != expected_count:
        logger.error(
            f"Validation failed for {genome_name}: "
            f"Found {non_na_count} non-NA values, expected {expected_count} unique proteins from AMR file"
        )
        return False
    
    logger.info(f"    - Validation passed: {non_na_count} AMR proteins matched")
    return True


def process_single_parquet(
    parquet_path: Path,
    amrf_results_dir: Path,
    output_dir: Path
) -> bool:
    """
    Process a single parquet file and merge AMR data.
    
    Args:
        parquet_path: Path to input parquet file
        amrf_results_dir: Directory containing AMR result TSV files
        output_dir: Root directory for output files
    
    Returns:
        True if successful, False otherwise
    """
    try:
        logger.info(f"Processing: {parquet_path}")
        
        # Read parquet file with pyarrow engine
        df = pd.read_parquet(parquet_path, engine='pyarrow')
        logger.info(f"Loaded {len(df)} genomes from {parquet_path.name}")
        
        # Pre-create new columns - collect data for all rows
        new_columns = [
            'element_symbol', 'element_name', 'amrf_type', 'amrf_subtype',
            'amrf_class', 'amrf_subclass', 'pct_identity_to_reference'
        ]
        
        # Collect all AMR data for each row
        all_amr_data = {col: [] for col in new_columns}
        
        # Process each genome row
        for idx, row in df.iterrows():
            genome_name = row['genome_name']
            protein_ids_nested = row['protein_id']
            
            # Step 1: Flatten protein IDs
            protein_ids_flat = flatten_protein_ids(protein_ids_nested, genome_name)
            
            # Step 2: Load AMR data for this genome
            amr_df = load_amr_data(genome_name, amrf_results_dir)
            
            # Step 3: Map AMR data to protein positions
            amr_data = map_amr_to_proteins(protein_ids_flat, amr_df, genome_name)
            
            # Step 4: Validate
            if not validate_amr_data(protein_ids_flat, amr_data, amr_df, genome_name):
                logger.error(f"Validation failed for {genome_name} in {parquet_path}")
                return False
            
            # Step 5: Collect arrays for this row
            for col_name, arr in amr_data.items():
                all_amr_data[col_name].append(arr)
        
        # Assign all collected data to dataframe columns
        for col_name, data_list in all_amr_data.items():
            df[col_name] = data_list
        
        # Create output path maintaining directory structure
        # Get relative path from input base directory
        relative_path = parquet_path.relative_to(
            Path('/home/dca36/rds/hpc-work/data/BacFormer/processed/ecoli_embeddings_with_ast_labels_flattened')
        )
        output_path = output_dir / relative_path
        
        # Create output directory if needed
        output_path.parent.mkdir(parents=True, exist_ok=True)
        
        # Write output parquet
        df.to_parquet(output_path, engine='pyarrow', index=False)
        logger.info(f"Successfully wrote output to: {output_path}")
        
        return True
        
    except Exception as e:
        logger.error(f"Error processing {parquet_path}: {e}", exc_info=True)
        return False


if __name__ == "__main__":
    # Test with a single file
    test_parquet = Path(
        '/home/dca36/rds/hpc-work/data/BacFormer/processed/'
        'ecoli_embeddings_with_ast_labels_flattened/test/batch_020_chunks/train_chunk_1.parquet'
    )
    amrf_dir = Path('/home/dca36/rds/hpc-work/data/amr_finder_plus/data/amrf_results')
    output_dir = Path('/home/dca36/rds/hpc-work/data/BacFormer/processed/ecoli_embeddings_with_amrf_labels_flattened')
    
    success = process_single_parquet(test_parquet, amrf_dir, output_dir)
    print(f"Test run {'succeeded' if success else 'failed'}")

