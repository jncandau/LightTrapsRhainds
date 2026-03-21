#!/bin/bash
#SBATCH --job-name=Hierarchical_fit
#SBATCH --time=48:00:00
#SBATCH --ntasks=1
#SBATCH --mem=256G
#SBATCH --cpus-per-task=4   # Stan chains are single-threaded by default
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jean.noel.candau@gmail.com
 
module load r/4.5.0
 
Rscript moths_hierarchical_fit.R 
 
exit
