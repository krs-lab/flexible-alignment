#!/bin/bash
#SBATCH --job-name=colabfold_multifasta_msa
##SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=200G
#SBATCH --time=48:00:00
#SBATCH --output=slurm_%j.log
#SBATCH --error=slurm_%j.err

set -e
set -x

chmod +x ./msa_biomix.sh
./msa_biomix.sh
