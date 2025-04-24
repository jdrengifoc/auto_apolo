#!/bin/bash
#SBATCH --partition=longjobs
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --time=04-00:00:00
#SBATCH --job-name=job_app=usBanks
#SBATCH --output=abc4sfa/3b_gibbs_usbanks//apolo_info/job_app=usBanks.out
#SBATCH --error=abc4sfa/3b_gibbs_usbanks//apolo_info/job_app=usBanks.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mmanzurg@eafit.edu.co

echo "Wed Apr 23 01:14:03 -05 2025: Cargando módulo python/3.6.0_miniconda-4.3.11_gcc-11.2.0"
module load python/3.6.0_miniconda-4.3.11_gcc-11.2.0
echo "Wed Apr 23 01:14:03 -05 2025: activando venv abc4sfa"
source activate abc4sfa

echo "Wed Apr 23 01:14:03 -05 2025: Ejecutando Rscript abc4sfa/3b_gibbs_usbanks//jobs/main_app=usBanks.R"
Rscript abc4sfa/3b_gibbs_usbanks//jobs/main_app=usBanks.R
