#!/bin/bash

#SBATCH --account=biostat
#SBATCH --partition=all-12c128g
#SBATCH --time=48:00:00

#SBATCH --job-name=eryutv20
#SBATCH --mem-per-cpu=50gb

#SBATCH --error=/home/users/%u/kzlinlab/projects/veloUncertainty/out/yuhong/Writeup3_eryutv/slurm_%j.err
#SBATCH --output=/home/users/%u/kzlinlab/projects/veloUncertainty/out/yuhong/Writeup3_eryutv/slurm_%j.out

source ~/miniconda3/etc/profile.d/conda.sh

conda activate utvClone
conda env list

python --version
conda list
python /home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/code/yuhong/erythroid/1_countsplit_and_umap/tmp_ery_utv_code/tmp_ery_utv_umap_320.py > /home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/Writeup3_eryutv/$SLURM_JOB_NAME-$SLURM_JOB_ID.txt