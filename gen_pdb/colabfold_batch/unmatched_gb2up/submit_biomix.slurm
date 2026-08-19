#!/bin/bash
#SBATCH --job-name=colabfold_multifasta
#SBATCH --partition=gpu
#SBATCH --qos=gpu
#SBATCH --gres=gpu:1
#SBATCH --nodelist=biomix64,biomix66,biomix61,biomix62,biomix63,biomix65,biomix67,biomix68,biomix70
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --output=slurm_%j.log
#SBATCH --error=slurm_%j.err

set -e

chmod +x ./batch_biomix.sh
./batch_biomix.sh
