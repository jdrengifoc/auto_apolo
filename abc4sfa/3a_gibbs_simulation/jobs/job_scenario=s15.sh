#!/bin/bash
#SBATCH --partition=longjobs
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --time=02-00:00:00
#SBATCH --job-name=job_scenario=s15
#SBATCH --output=abc4sfa/3a_gibbs_simulation//apolo_info/job_scenario=s15.out
#SBATCH --error=abc4sfa/3a_gibbs_simulation//apolo_info/job_scenario=s15.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mmanzurg@eafit.edu.co

echo "Thu Apr 24 11:02:09 -05 2025: Cargando módulo python/3.6.0_miniconda-4.3.11_gcc-11.2.0"
module load python/3.6.0_miniconda-4.3.11_gcc-11.2.0
echo "Thu Apr 24 11:02:09 -05 2025: activando venv abc4sfa"
source activate abc4sfa

echo "Thu Apr 24 11:02:09 -05 2025: Ejecutando Rscript abc4sfa/3a_gibbs_simulation//jobs/main_scenario=s15.R"
Rscript abc4sfa/3a_gibbs_simulation//jobs/main_scenario=s15.R
