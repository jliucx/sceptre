#!/bin/bash

#SBATCH --job-name=nature_data
#SBATCH --cpus-per-task=5
#SBATCH --ntasks-per-node=1
#SBATCH --mem=100GB
#SBATCH --time=0:30:00
#SBATCH --partition=xuanyao-hm
#SBATCH --qos=xuanyao
#SBATCH --account=pi-xuanyao
#SBATCH --output=rcc_out/seruat_%j.out
#SBATCH --error=rcc_err/seruat_%j.err

module load R/4.1.0
Rscript rcc_jobs/t_test.R

echo "Job completed successfully"

