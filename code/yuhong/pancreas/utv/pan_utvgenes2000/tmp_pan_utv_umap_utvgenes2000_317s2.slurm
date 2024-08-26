#!/bin/bash

#SBATCH --account=biostat
#SBATCH --partition=all-12c128g
#SBATCH --time=48:00:00

#SBATCH --job-name=panutv7s22000
#SBATCH --mem-per-cpu=50gb

#SBATCH --error=/home/users/%u/kzlinlab/projects/veloUncertainty/out/yuhong/Writeup5_panutv_utvgenes/slurm_%j.err
#SBATCH --output=/home/users/%u/kzlinlab/projects/veloUncertainty/out/yuhong/Writeup5_panutv_utvgenes/slurm_%j.out

source ~/miniconda3/etc/profile.d/conda.sh

conda activate utvClone
conda env list

python --version
conda list
python /home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/code/yuhong/pancreas/utv/pan_utvgenes2000/tmp_pan_utv_umap_utvgenes2000_317s2.py > /home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/Writeup5_panutv_utvgenes/$SLURM_JOB_NAME-$SLURM_JOB_ID.txt
