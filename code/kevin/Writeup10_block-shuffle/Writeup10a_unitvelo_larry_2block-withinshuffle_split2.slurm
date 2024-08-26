#!/bin/bash

#SBATCH --account=biostat
#SBATCH --partition=all-20c96g
#SBATCH --time=48:00:00

#SBATCH --job-name=unitvelo_larry_2block-ws_s2
#SBATCH --mem-per-cpu=20gb

#SBATCH --error=/home/users/%u/kzlinlab/projects/veloUncertainty/out/kevin/Writeup10/slurm/slurm_%j.err
#SBATCH --output=/home/users/%u/kzlinlab/projects/veloUncertainty/out/kevin/Writeup10/slurm/slurm_%j.out

source ~/miniconda3/etc/profile.d/conda.sh

conda activate unitvelo
conda env list

cd /home/users/kzlin/kzlinlab/projects/veloUncertainty/out/kevin/Writeup10/tmp3

python --version
conda list
python /home/users/kzlin/kzlinlab/projects/veloUncertainty/git/veloUncertainty_kevin/code/kevin/Writeup10_block-shuffle/Writeup10a_unitvelo_larry_2block-withinshuffle_split2.py > /home/users/kzlin/kzlinlab/projects/veloUncertainty/out/kevin/Writeup10/slurm/$SLURM_JOB_NAME-$SLURM_JOB_ID.txt

