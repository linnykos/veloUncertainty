#!/bin/bash

#SBATCH --account=biostat
#SBATCH --partition=all-12c128g
#SBATCH --time=48:00:00

#SBATCH --job-name=sctour_shuffle-all
#SBATCH --mem-per-cpu=100gb

#SBATCH --error=/home/users/%u/kzlinlab/projects/veloUncertainty/out/kevin/Writeup9/slurm/slurm_%j.err
#SBATCH --output=/home/users/%u/kzlinlab/projects/veloUncertainty/out/kevin/Writeup9/slurm/slurm_%j.out

source ~/miniconda3/etc/profile.d/conda.sh

conda activate sctour
conda env list

python --version
conda list
python Writeup9_sctour_larry_shuffle-all.py > /home/users/kzlin/kzlinlab/projects/veloUncertainty/out/kevin/Writeup9/slurm/$SLURM_JOB_NAME-$SLURM_JOB_ID.txt