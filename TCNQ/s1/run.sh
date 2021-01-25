#!/bin/bash
#SBATCH -n 64                
#SBATCH --ntasks-per-node=64 
#SBATCH -p medium            
#SBATCH --job-name=SlabCalc

cd $SLURM_SUBMIT_DIR

mpirun -np $SLURM_NTASKS pw.x < TCNQ-s1.in 1> TCNQ-s1.out 2> err
