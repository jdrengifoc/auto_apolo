#!/bin/bash
#SBATCH --partition=longjobs
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=07-00:00:00
#SBATCH --job-name=job_app_name=usBanks
#SBATCH --output=tristeza_empirica/apolo_info/job_app_name=usBanks.out
#SBATCH --error=tristeza_empirica/apolo_info/job_app_name=usBanks.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mmanzurg@eafit.edu.co

echo "Wed Mar 19 19:08:37 -05 2025: Cargando módulo python/3.6.0_miniconda-4.3.11_gcc-11.2.0"
module load python/3.6.0_miniconda-4.3.11_gcc-11.2.0
echo "Wed Mar 19 19:08:37 -05 2025: activando venv abc4sfa"
source activate abc4sfa

echo "Wed Mar 19 19:08:37 -05 2025: Ejecutando Rscript tristeza_empirica/jobs/main_app_name=usBanks.R"
Rscript tristeza_empirica/jobs/main_app_name=usBanks.R
