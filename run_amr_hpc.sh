#!/bin/bash
#SBATCH --job-name=amr_finder_merge
#SBATCH --partition=icelake-himem 
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=76
#SBATCH --time=00:15:00
#SBATCH --output=amr_finder_merge_%j.out
#SBATCH --error=amr_finder_merge_%j.err
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

#python amrf/parallel_process_batches.py $SLURM_CPUS_PER_TASK
uv run python amrf/run_amr_merge_multiprocess.py

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
