#!/bin/bash
#SBATCH -p main
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --mem 200gb 
#SBATCH -t 24:00 
#SBATCH -J project_name
#SBATCH --output=project_name.o

hostname
env
ulimit -v unlimited 
ulimit -a

source $HOME/.bashrc
module load R/4.2.0
module list
date
Rscript project_name.R 
date