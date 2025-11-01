#!/bin/bash
#SBATCH --job-name=amrf_acinetobacter_baumannii
#SBATCH --partition=icelake-himem 
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=76
#SBATCH --time=06:00:00
#SBATCH --output=amrf_acinetobacter_baumannii_%j.out
#SBATCH --error=amrf_acinetobacter_baumannii_%j.err
#SBATCH --account=FLOTO-SL3-CPU

# AMR Coverage Generation HPC Processing Script
# Generates missing AMR files for target species genomes

echo "=========================================="
echo "MERGE AMR DATA - HPC PROCESSING"
echo "=========================================="
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $SLURM_NODELIST"
echo "CPUs: $SLURM_CPUS_PER_TASK"
echo "Start time: $(date)"
echo "=========================================="

# Load conda environment
module purge
module load miniconda/3
source activate my-bioinformatics

# Verify AMRFinder is available
echo "Checking AMRFinder installation..."
amrfinder --version || { echo "Error: AMRFinder not found"; exit 1; }

# Navigate to workspace
cd /home/dca36/workspace/AMR_Finder_Plus

# Generate missing AMR files for target species genomes
echo ""
echo "Starting AMR generation for acinetobacter_baumannii..."
echo ""

# Uses parallel_amr_finder_plus_process_batches.py to process all parquet files
uv sync # Ensures all recently added dependencies are updated
cd /home/dca36/workspace/AMR_Finder_Plus 
uv run python amrf/parallel_amr_finder_plus_process_batches.py \
    --input-dir /home/dca36/rds/hpc-work/data/BacFormer/raw/ast/protein_sequences/acinetobacter_baumannii \
    --species acinetobacter_baumannii \
    --fasta-dir /home/dca36/rds/hpc-work/data/amr_finder_plus/data/fasta_files \
    --results-dir /home/dca36/rds/hpc-work/data/amr_finder_plus/data/amrf_results

# uv run python amrf/run_amr_merge_multiprocess.py \
#     --input-dir /home/dca36/rds/hpc-work/data/BacFormer/processed/ast_esm_embeddings/campylobacter_jejuni \
#     --output-dir /home/dca36/rds/hpc-work/data/BacFormer/processed/ast_esm_embeddings_with_amrf/campylobacter_jejuni \
#     --amrf-dir /home/dca36/rds/hpc-work/data/amr_finder_plus/data/amrf_results

# Check exit status
EXIT_CODE=$?

echo ""
echo "=========================================="
echo "Job completed with exit code: $EXIT_CODE"
echo "End time: $(date)"
echo "=========================================="

exit $EXIT_CODE

# sbatch run_on_hpc.sh
# squeue -u dca36
# scancel -u dca36
