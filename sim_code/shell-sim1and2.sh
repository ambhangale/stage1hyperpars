#!/bin/bash

#SBATCH -J sim1and2
#SBATCH -e sim1and2.SERR
#SBATCH -o sim1and2.SOUT
#SBATCH -N 1
#SBATCH --cpus-per-task 32
#SBACTH -t 3-20:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=aditibhangale@gmail.com

cd $SLURM_SUBMIT_DIR

module load R/4.2.2
export MKL_NUM_THREADS=1

Rscript --vanilla runsim-sim1to3.R
