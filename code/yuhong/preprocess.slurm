#!/bin/bash

#SBATCH --account=biostat
#SBATCH --partition=all-12c128g
#SBATCH --time=48:00:00

#SBATCH --job-name=veloSkeleton
#SBATCH --mem-per-cpu=50gb

#SBATCH --error=/home/users/%u/kzlinlab/projects/veloUncertainty/out/yuhong/Writeup1/slurm/slurm_%j.err
#SBATCH --output=/home/users/%u/kzlinlab/projects/veloUncertainty/out/yuhong/Writeup1/slurm/slurm_%j.out

conda init
conda activate scVelo
python --version
conda list
python preprocess.py > /home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/Writeup1/slurm/$SLURM_JOB_NAME-$SLURM_JOB_ID.txt