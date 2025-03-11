#!/bin/sh
#SBATCH --partition=longjobs
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --time=14-00:00:00
#SBATCH --job-name=test_job
#SBATCH -o %x_%j.out      # File to which STDOUT will be written
#SBATCH -e %x_%j.err      # File to which STDERR will be written
#SBATCH --mail-type=ALL
#SBATCH --mail-user=<usuario>@eafit.edu.co


## se cargan los módulos y el entorno
module load python/3.9_miniconda-4.9.2

source activate entorno


## se corren los programas 
python prog.py

