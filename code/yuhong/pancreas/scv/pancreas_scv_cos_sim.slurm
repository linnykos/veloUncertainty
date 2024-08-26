#!/bin/bash

#SBATCH --account=biostat
#SBATCH --partition=all-12c128g
#SBATCH --time=48:00:00

#SBATCH --job-name=panscv
#SBATCH --mem-per-cpu=50gb

#SBATCH --error=/home/users/%u/kzlinlab/projects/veloUncertainty/out/yuhong/Writeup1/slurm_scv/slurm_%j.err
#SBATCH --output=/home/users/%u/kzlinlab/projects/veloUncertainty/out/yuhong/Writeup1/slurm_scv/slurm_%j.out

source ~/miniconda3/etc/profile.d/conda.sh

conda activate scVelo
conda env list

python --version
conda list
python /home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/code/yuhong/pancreas/scv/pancreas_scv_cos_sim.py > /home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/Writeup1/slurm_scv/$SLURM_JOB_NAME-$SLURM_JOB_ID.txt

