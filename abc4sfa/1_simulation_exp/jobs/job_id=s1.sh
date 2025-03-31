#!/bin/bash
#SBATCH --partition=longjobs
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --time=08-00:00:00
#SBATCH --job-name=job_id=s1
#SBATCH --output=abc4sfa/1_simulation_hn/apolo_info/job_id=s1.out
#SBATCH --error=abc4sfa/1_simulation_hn/apolo_info/job_id=s1.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mmanzurg@eafit.edu.co

echo "Tue Mar 25 22:38:54 -05 2025: Cargando módulo python/3.6.0_miniconda-4.3.11_gcc-11.2.0"
module load python/3.6.0_miniconda-4.3.11_gcc-11.2.0
echo "Tue Mar 25 22:38:54 -05 2025: activando venv abc4sfa"
source activate abc4sfa

echo "Tue Mar 25 22:38:54 -05 2025: Ejecutando Rscript abc4sfa/1_simulation_hn/jobs/main_id=s1.R"
Rscript abc4sfa/1_simulation_hn/jobs/main_id=s1.R
