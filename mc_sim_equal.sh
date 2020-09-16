#!/bin/bash -l

#SBATCH --account=coexistence
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --cpus-per-task=1
#SBATCH --mem=64GB
#SBATCH --mail-type=ALL
#SBATCH --mail-user=nathan.wisnoski@uwyo.edu
#SBATCH --job-name=metacom_equal

cd /home/nwisnosk/GitHub/metacom-coexistence

module load gcc/7.3.0 r/3.6.1

R CMD BATCH --no-restore --no-save --quiet mc_sim_equal.R mc_sim_equal.log
