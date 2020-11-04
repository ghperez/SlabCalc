#!/bin/bash
#SBATCH -n 128                
#SBATCH --ntasks-per-node=64 
#SBATCH -p medium            
#SBATCH --job-name=SlabCalc

cd $SLURM_SUBMIT_DIR

export PYTHONPATH=/home/gabriel.perez/miniconda3/envs/molSimplify/bin/python

python main.py $SLURM_NTASKS 1> calc.out 2> calc.err

