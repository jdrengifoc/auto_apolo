#!/bin/bash
#SBATCH --partition=longjobs
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=07-00:00:00
#SBATCH --job-name=job_app_name=usElectricity
#SBATCH --output=tristeza_empirica//apolo_info/job_app_name=usElectricity.out
#SBATCH --error=tristeza_empirica//apolo_info/job_app_name=usElectricity.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mmanzurg@eafit.edu.co

echo "Thu Mar 13 22:41:18 -05 2025: Cargando módulo python/3.6.0_miniconda-4.3.11_gcc-11.2.0"
module load python/3.6.0_miniconda-4.3.11_gcc-11.2.0
echo "Thu Mar 13 22:41:18 -05 2025: activando venv abc4sfa"
source activate abc4sfa

echo "Thu Mar 13 22:41:18 -05 2025: Ejecutando Rscript tristeza_empirica//jobs/main_app_name=usElectricity.R"
Rscript tristeza_empirica//jobs/main_app_name=usElectricity.R
