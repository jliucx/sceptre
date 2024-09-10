#!/bin/bash

#SBATCH --job-name=run_sceptre
#SBATCH --cpus-per-task=5
#SBATCH --ntasks-per-node=1
#SBATCH --mem=200GB
#SBATCH --time=10:00:00
#SBATCH --partition=xuanyao-hm
#SBATCH --qos=xuanyao
#SBATCH --account=pi-xuanyao
#SBATCH --output=job-info.out
#SBATCH --error=job-info.err

module load R/4.1.0
Rscript try.R

echo "Job completed successfully"

