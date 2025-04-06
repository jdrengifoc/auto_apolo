#!/bin/bash
#SBATCH --partition=longjobs
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=07-00:00:00
#SBATCH --job-name=job_app_name=usBanks_time_limit_minutes=900
#SBATCH --output=abc4sfa/2a_empirical_applications_sa/apolo_info/job_app_name=usBanks_time_limit_minutes=900.out
#SBATCH --error=abc4sfa/2a_empirical_applications_sa/apolo_info/job_app_name=usBanks_time_limit_minutes=900.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mmanzurg@eafit.edu.co

echo "Tue Apr  1 23:04:41 -05 2025: Cargando módulo python/3.6.0_miniconda-4.3.11_gcc-11.2.0"
module load python/3.6.0_miniconda-4.3.11_gcc-11.2.0
echo "Tue Apr  1 23:04:41 -05 2025: activando venv abc4sfa"
source activate abc4sfa

echo "Tue Apr  1 23:04:41 -05 2025: Ejecutando Rscript abc4sfa/2a_empirical_applications_sa/jobs/main_app_name=usBanks_time_limit_minutes=900.R"
Rscript abc4sfa/2a_empirical_applications_sa/jobs/main_app_name=usBanks_time_limit_minutes=900.R
