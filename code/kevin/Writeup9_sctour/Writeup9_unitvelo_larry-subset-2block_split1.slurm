#!/bin/bash

#SBATCH --account=biostat
#SBATCH --partition=all-20c96g
#SBATCH --time=48:00:00

#SBATCH --job-name=unitvelo_larry-subset-2b_s1
#SBATCH --mem-per-cpu=15gb

#SBATCH --error=/home/users/%u/kzlinlab/projects/veloUncertainty/out/kevin/Writeup9/slurm/slurm_%j.err
#SBATCH --output=/home/users/%u/kzlinlab/projects/veloUncertainty/out/kevin/Writeup9/slurm/slurm_%j.out

source ~/miniconda3/etc/profile.d/conda.sh

conda activate unitvelo
conda env list

python --version
conda list
python Writeup9_unitvelo_larry-subset-2block_split1.py > /home/users/kzlin/kzlinlab/projects/veloUncertainty/out/kevin/Writeup9/slurm/$SLURM_JOB_NAME-$SLURM_JOB_ID.txt

