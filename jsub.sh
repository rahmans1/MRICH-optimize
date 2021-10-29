#!/bin/bash
#SBATCH --account=rrg-wdconinc
#SBATCH --job-name=mrich
#SBATCH --time=24:00:00 
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=4G
#SBATCH --output==%x-%j.out
cd ~/projects/rrg-wdconinc/rahmans/MRICH-optimization
source  OPTENV/bin/activate
module load python/3.9.6
python optimize.py > results.txt
