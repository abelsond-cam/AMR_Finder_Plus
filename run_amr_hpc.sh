#!/bin/bash
#SBATCH --job-name=amr_finder_parallel_chunks
#SBATCH --partition=icelake-himem 
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=76
#SBATCH --time=05:00:00
#SBATCH --output=amr_finder_parallel_chunks_%j.out
#SBATCH --error=amr_finder_parallel_chunks_%j.err
#SBATCH --account=FLOTO-SL2-CPU

# AMR Finder Plus HPC Processing Script
# Processes all E. coli batch parquet files in parallel

echo "=========================================="
echo "AMR FINDER PLUS - HPC PROCESSING"
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

# Run parallel processing script using all available CPUs
echo ""
echo "Starting parallel processing with $SLURM_CPUS_PER_TASK CPUs..."
echo ""

# Run AMR Finder Plus on parquet files
# Default: processes /home/dca36/rds/hpc-work/data/BacFormer/raw/ast/protein_sequences/all_species
# Outputs to: /home/dca36/rds/hpc-work/data/amr_finder_plus/data/amrf_results
# Each parquet file is split into 5 chunks for optimal CPU utilization (27 files × 5 = 135 chunks for 76 CPUs)
uv sync # Ensures all recently added dependencies are updated)
uv run python amrf/parallel_process_batches.py --cpus $SLURM_CPUS_PER_TASK --chunks-per-file 5

# To customize directories or chunking, use:
# uv run python amrf/parallel_process_batches.py \
#     --input-dir /path/to/parquets \
#     --fasta-dir /path/to/fastas \
#     --results-dir /path/to/results \
#     --cpus $SLURM_CPUS_PER_TASK \
#     --chunks-per-file 3

# Check exit status
EXIT_CODE=$?

echo ""
echo "=========================================="
echo "Job completed with exit code: $EXIT_CODE"
echo "End time: $(date)"
echo "=========================================="

exit $EXIT_CODE

# sbatch run_amr_hpc.sh
# squeue -u dca36
# scancel -u dca36
