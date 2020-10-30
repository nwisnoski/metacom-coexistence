#!/bin/bash -l

#SBATCH --account=coexistence
#SBATCH --time=2:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --cpus-per-task=1
#SBATCH --mem=64GB
#SBATCH --mail-type=ALL
#SBATCH --mail-user=nathan.wisnoski@uwyo.edu
#SBATCH --job-name=make_plots

cd /gscratch/nwisnosk/GitHub/metacom-coexistence

module load gcc/7.3.0 r/3.6.1

R CMD BATCH --no-restore --no-save  make_plots.R make_plots.log
