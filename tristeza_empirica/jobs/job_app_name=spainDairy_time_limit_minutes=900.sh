#!/bin/bash
#SBATCH --partition=longjobs
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=07-00:00:00
#SBATCH --job-name=job_app_name=spainDairy_time_limit_minutes=900
#SBATCH --output=tristeza_empirica/apolo_info/job_app_name=spainDairy_time_limit_minutes=900.out
#SBATCH --error=tristeza_empirica/apolo_info/job_app_name=spainDairy_time_limit_minutes=900.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mmanzurg@eafit.edu.co

echo "Sat Mar 22 10:43:40 -05 2025: Cargando módulo python/3.6.0_miniconda-4.3.11_gcc-11.2.0"
module load python/3.6.0_miniconda-4.3.11_gcc-11.2.0
echo "Sat Mar 22 10:43:40 -05 2025: activando venv abc4sfa"
source activate abc4sfa

echo "Sat Mar 22 10:43:40 -05 2025: Ejecutando Rscript tristeza_empirica/jobs/main_app_name=spainDairy_time_limit_minutes=900.R"
Rscript tristeza_empirica/jobs/main_app_name=spainDairy_time_limit_minutes=900.R
