#!/bin/bash

#SBATCH --account=biostat
#SBATCH --partition=all-20c96g
#SBATCH --time=48:00:00

#SBATCH --job-name=greenleaf
#SBATCH --mem-per-cpu=50gb

#SBATCH --error=/home/users/%u/kzlinlab/projects/veloUncertainty/out/kevin/Writeup15/slurm/slurm_%j.err
#SBATCH --output=/home/users/%u/kzlinlab/projects/veloUncertainty/out/kevin/Writeup15/slurm/slurm_%j.out

source ~/miniconda3/etc/profile.d/conda.sh

conda activate scVI
conda env list

python --version
conda list
python Writeup15_greenleaf.py > /home/users/kzlin/kzlinlab/projects/veloUncertainty/out/kevin/Writeup15/slurm/$SLURM_JOB_NAME-$SLURM_JOB_ID.txt