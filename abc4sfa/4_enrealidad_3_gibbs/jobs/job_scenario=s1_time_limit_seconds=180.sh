#!/bin/bash
#SBATCH --partition=longjobs
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --time=07-00:00:00
#SBATCH --job-name=job_scenario=s1_time_limit_seconds=180
#SBATCH --output=abc4sfa/4_enrealidad_3_gibbs/apolo_info/job_scenario=s1_time_limit_seconds=180.out
#SBATCH --error=abc4sfa/4_enrealidad_3_gibbs/apolo_info/job_scenario=s1_time_limit_seconds=180.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mmanzurg@eafit.edu.co

echo "Mon Apr  7 23:47:34 -05 2025: Cargando módulo python/3.6.0_miniconda-4.3.11_gcc-11.2.0"
module load python/3.6.0_miniconda-4.3.11_gcc-11.2.0
echo "Mon Apr  7 23:47:34 -05 2025: activando venv abc4sfa"
source activate abc4sfa

echo "Mon Apr  7 23:47:34 -05 2025: Ejecutando Rscript abc4sfa/4_enrealidad_3_gibbs/jobs/main_scenario=s1_time_limit_seconds=180.R"
Rscript abc4sfa/4_enrealidad_3_gibbs/jobs/main_scenario=s1_time_limit_seconds=180.R
