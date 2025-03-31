#!/bin/bash
#SBATCH --partition=longjobs
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=07-00:00:00
#SBATCH --job-name=job_app_name=spainDairy_time_limit_minutes=20000
#SBATCH --output=abc4sfa/2_empirical_applications/apolo_info/job_app_name=spainDairy_time_limit_minutes=20000.out
#SBATCH --error=abc4sfa/2_empirical_applications/apolo_info/job_app_name=spainDairy_time_limit_minutes=20000.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mmanzurg@eafit.edu.co

echo "Tue Mar 25 22:38:35 -05 2025: Cargando módulo python/3.6.0_miniconda-4.3.11_gcc-11.2.0"
module load python/3.6.0_miniconda-4.3.11_gcc-11.2.0
echo "Tue Mar 25 22:38:35 -05 2025: activando venv abc4sfa"
source activate abc4sfa

echo "Tue Mar 25 22:38:35 -05 2025: Ejecutando Rscript abc4sfa/2_empirical_applications/jobs/main_app_name=spainDairy_time_limit_minutes=20000.R"
Rscript abc4sfa/2_empirical_applications/jobs/main_app_name=spainDairy_time_limit_minutes=20000.R
