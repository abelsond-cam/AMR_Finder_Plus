# AMR Finder Plus Integration Plan

## Project Overview
Integrate AMR Finder Plus to annotate ~5000 E. coli protein sequences with antimicrobial resistance information. The goal is to add AMR annotations to existing protein embeddings while maintaining the contig and MAG structure.

## Data Structure Understanding

### Source Files
- **Protein Sequences**: `/home/dca36/rds/hpc-work/data/BacFormer/processed/ecoli_batches/ecoli_batch_000.parquet`
  - Contains: `genome_name`, `contig_name`, `protein_id`, `protein_sequence`, `taxid`, `locus_tag`, `start`, `end`, `product`
  - Shape: (250, 9) - this is just one batch, need to process all batches
  - `protein_sequence` column contains a list of lists [[xxxx yyyy zzzz][aaa bbb cccc]] with each element in the internal list being an amino acid sequence as a numpy array.  Each elemnt of top list corresponsds to a contig in the genome.  Each internal element on a contig list is a protein with a amino acids.  protein_id has the same structure, but with id's for each protein. We will flatten the protein sequence into a fasta file and flatten the protein_id list into a flat list of ids to use after the run of amrfinderplus. 

- **Target File (Test)**: `/home/dca36/rds/hpc-work/data/BacFormer/processed/ecoli_embeddings_with_ast_labels_flattened/train/batch_000_chunks/train_chunk_1.parquet`
  - Contains: protein embeddings + AST labels + metadata
  - Shape: (50, 47) 
  - Key linking field: `protein_id`

### Data Mapping Strategy
- Use `protein_id` to match sequences from ecoli_batches to embeddings files
- Preserve all existing data structure and add new AMR annotation columns

## AMR Finder Plus Installation

### Prerequisites
```bash
# Install conda/miniconda if not available
cd /home/dca36/workspace/AMR_Finder_Plus - DONE
```

### Installation Steps
```bash
# Create dedicated conda environment
conda create -n my-bioinformatics -c bioconda ncbi-amrfinderplus python=3.9 pandas pyarrow - DONE
conda activate my-bioinformatics - DONE

# Update AMR Finder Plus database
amrfinder -u - DONE

# Test installation
amrfinder --help
```

## AMR Finder Plus Output Columns

Based on AMR Finder Plus documentation, we'll capture these key fields as separate parquet columns:

1. **AMR_gene_symbol** - Gene symbol (e.g., "blaOXA-1")
2. **AMR_sequence_name** - Protein/sequence identifier  
3. **AMR_scope** - Scope (e.g., "core", "plus")
4. **AMR_element_type** - Type (e.g., "AMR", "VIRULENCE", "STRESS")
5. **AMR_element_subtype** - Subtype (e.g., "AMR", "VIRULENCE")
6. **AMR_class** - Resistance class (e.g., "BETA-LACTAM")
7. **AMR_subclass** - Resistance subclass
8. **AMR_method** - Detection method (e.g., "EXACTX", "BLASTX")
9. **AMR_target_length** - Target sequence length
10. **AMR_reference_length** - Reference sequence length
11. **AMR_percent_coverage** - Percent coverage of reference
12. **AMR_percent_identity** - Percent identity to reference
13. **AMR_alignment_length** - Length of alignment
14. **AMR_accession** - Reference accession number

For proteins with no AMR hits, all columns will contain `NaN`.

## Implementation Plan

### Phase 1: Setup and Testing (Single Sample)
```bash
cd /home/dca36/workspace/AMR_Finder_Plus - DONE
```

1. **Create test script** (`test_single_protein.py`):
   - Read one protein sequence from ecoli_batch_000.parquet
   - Write to temporary FASTA file
   - Run AMRFinderPlus via subprocess
   - Parse output and display results

2. **Test AMR Finder Plus**:
   ```bash
   conda activate amrfinder
   uv run python test_single_protein.py
   ```

### Phase 2: Full Implementation

1. **Create main processing script** (`process_amr_annotations.py`):
   - Read all protein sequences from ecoli batches
   - Process each sequence through AMRFinderPlus
   - Create mapping dictionary: protein_id -> AMR results
   - Handle batch processing for memory efficiency

2. **Create annotation script** (`add_amr_to_embeddings.py`):
   - Read target parquet file (embeddings with AST labels)
   - Match protein_ids with AMR results
   - Add 14 new AMR columns
   - Save updated parquet file

### Phase 3: Batch Processing

1. **Process all ecoli batches**:
   ```bash
   # Find all ecoli batch files
   find /home/dca36/rds/hpc-work/data/BacFormer/processed/ecoli_batches/ -name "*.parquet"
   ```

2. **Apply to all embedding files**:
   ```bash
   # Find all embedding files that need annotation
   find /home/dca36/rds/hpc-work/data/BacFormer/processed/ecoli_embeddings_with_ast_labels_flattened/ -name "*.parquet"
   ```

## File Organization

```
/home/dca36/workspace/AMR_Finder_Plus/
├── README_AMR_F_PLAN.md                 # This file
├── environment.yml                      # Conda environment spec
├── test_single_protein.py              # Phase 1: Single protein test
├── process_amr_annotations.py          # Phase 2: Batch AMR processing
├── add_amr_to_embeddings.py            # Phase 2: Add columns to parquet
├── run_full_pipeline.py                # Phase 3: Complete pipeline
├── temp/                               # Temporary FASTA files
│   ├── input_sequences.fasta
│   └── amr_output.tsv
└── logs/                               # Processing logs
    └── amr_processing.log
```

## Expected Workflow

### For Testing (Single Sample):
```bash
conda activate amrfinder
cd /home/dca36/workspace/AMR_Finder_Plus
uv run python test_single_protein.py
```

### For Full Processing:
```bash
conda activate amrfinder
cd /home/dca36/workspace/AMR_Finder_Plus

# Step 1: Process all protein sequences through AMR Finder Plus
uv run python process_amr_annotations.py

# Step 2: Add AMR annotations to embedding files
uv run python add_amr_to_embeddings.py

# Or run complete pipeline:
uv run python run_full_pipeline.py
```

## Key Technical Details

### AMR Finder Plus Command Structure:
```bash
amrfinder --protein input_sequences.fasta --organism Escherichia --output amr_results.tsv --plus
```

### Important Parameters:
- `--protein`: Input protein sequences in FASTA format
- `--organism Escherichia`: Specify E. coli for species-specific database
- `--plus`: Include additional resistance genes beyond core set
- `--output`: Output TSV file with results

### Error Handling:
- Handle sequences that fail AMR analysis
- Manage temporary file cleanup
- Log processing progress and errors
- Validate output before updating parquet files

## Performance Considerations

- **Memory**: Process proteins in batches to avoid memory issues
- **Speed**: AMR Finder Plus can be slow; consider parallel processing
- **Storage**: Temporary FASTA files will be created and cleaned up
- **Validation**: Verify results before overwriting original data

## Next Steps

1. Set up conda environment in new workspace
2. Install AMR Finder Plus and dependencies
3. Implement and test single protein processing
4. Scale to batch processing
5. Validate results on sample data
6. Apply to full dataset

## Troubleshooting

- **Installation issues**: Check bioconda channel configuration
- **Database errors**: Ensure `amrfinder -u` completes successfully
- **Memory problems**: Reduce batch sizes in processing
- **File permissions**: Ensure write access to output directories

---

*Generated: $(date)*
*Workspace: /home/dca36/workspace/AMR_Finder_Plus*
