#!/bin/bash
#SBATCH --job-name=calibration
#SBATCH --partition=standard
#SBATCH -c 16
#SBATCH --mem=16G

# change to work directory
cd $SLURM_SUBMIT_DIR

# export R_LIBS with the desired version of calibrar
export R_LIBS=/home/$USER/R/x86_64-pc-linux-gnu-library/4.4

# Run R script in parallel mode
time Rscript calibrar_minimal.R --ncores=6 &>> mycalibration.log
