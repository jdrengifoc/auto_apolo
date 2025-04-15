#!/bin/bash
#SBATCH --partition=longjobs
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --time=01-00:00:00
#SBATCH --job-name=job_app=usBanks_time_limit_seconds=900
#SBATCH --output=abc4sfa/4_enrealidad_3_gibbs_usbank//apolo_info/job_app=usBanks_time_limit_seconds=900.out
#SBATCH --error=abc4sfa/4_enrealidad_3_gibbs_usbank//apolo_info/job_app=usBanks_time_limit_seconds=900.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mmanzurg@eafit.edu.co

echo "Tue Apr  8 19:55:11 -05 2025: Cargando módulo python/3.6.0_miniconda-4.3.11_gcc-11.2.0"
module load python/3.6.0_miniconda-4.3.11_gcc-11.2.0
echo "Tue Apr  8 19:55:11 -05 2025: activando venv abc4sfa"
source activate abc4sfa

echo "Tue Apr  8 19:55:11 -05 2025: Ejecutando Rscript abc4sfa/4_enrealidad_3_gibbs_usbank//jobs/main_app=usBanks_time_limit_seconds=900.R"
Rscript abc4sfa/4_enrealidad_3_gibbs_usbank//jobs/main_app=usBanks_time_limit_seconds=900.R
