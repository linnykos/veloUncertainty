#!/bin/bash

#SBATCH --account=biostat
#SBATCH --partition=all-20c96g
#SBATCH --time=48:00:00

#SBATCH --job-name=scv_larry-sub
#SBATCH --mem-per-cpu=15gb

#SBATCH --error=/home/users/%u/kzlinlab/projects/veloUncertainty/out/kevin/Writeup9/slurm/slurm_%j.err
#SBATCH --output=/home/users/%u/kzlinlab/projects/veloUncertainty/out/kevin/Writeup9/slurm/slurm_%j.out

source ~/miniconda3/etc/profile.d/conda.sh

conda activate scvelo
conda env list

cd /home/users/kzlin/kzlinlab/projects/veloUncertainty/out/kevin/Writeup9/tmp1

python --version
conda list
python /home/users/kzlin/kzlinlab/projects/veloUncertainty/git/veloUncertainty_kevin/code/kevin/Writeup9_sctour/Writeup9_scvelo_larry-subset.py > /home/users/kzlin/kzlinlab/projects/veloUncertainty/out/kevin/Writeup9/slurm/$SLURM_JOB_NAME-$SLURM_JOB_ID.txt

