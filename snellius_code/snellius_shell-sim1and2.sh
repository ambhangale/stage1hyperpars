"$TMP#!/bin/bash

#SBATCH -J sim1and2
#SBATCH -e sim1and2.SERR
#SBATCH -o sim1and2.SOUT
#SBATCH -p rome
#SBATCH -N 1
#SBATCH -n 32
#SBATCH -t 3-20:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=aditibhangale@gmail.com

cd "$TMPDIR"

module load 2022
module load R/4.2.1-foss-2022a
export MKL_NUM_THREADS=1

cp $HOME/SR-SEM/stage1hyperpars/functions-sim1to3.R "$TMPDIR"

Rscript --vanilla functions-sim1to3.R

cp "$TMPDIR"/*.rds $HOME/SR-SEM/stage1hyperpars/
