#!/bin/bash

#SBATCH --account=biostat
#SBATCH --partition=all-20c96g
#SBATCH --time=48:00:00

#SBATCH --job-name=unitvelo_larry_s1
#SBATCH --mem-per-cpu=20gb

#SBATCH --error=/home/users/%u/kzlinlab/projects/veloUncertainty/out/kevin/Writeup11/slurm/slurm_%j.err
#SBATCH --output=/home/users/%u/kzlinlab/projects/veloUncertainty/out/kevin/Writeup11/slurm/slurm_%j.out

source ~/miniconda3/etc/profile.d/conda.sh

conda activate unitvelo
conda env list

cd /home/users/kzlin/kzlinlab/projects/veloUncertainty/out/kevin/Writeup11/tmp2

python --version
conda list
python /home/users/kzlin/kzlinlab/projects/veloUncertainty/git/veloUncertainty_kevin/code/kevin/Writeup11_full-block-shuffle/Writeup11_unitvelo_larry_full-block-shuffle_split1.py > /home/users/kzlin/kzlinlab/projects/veloUncertainty/out/kevin/Writeup11/slurm/$SLURM_JOB_NAME-$SLURM_JOB_ID.txt

