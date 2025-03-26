#!/bin/bash
#SBATCH --partition=longjobs
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --time=08-00:00:00
#SBATCH --job-name=job_id=s1
#SBATCH --output=negacion_empirica/apolo_info/job_id=s1.out
#SBATCH --error=negacion_empirica/apolo_info/job_id=s1.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mmanzurg@eafit.edu.co

echo "Wed Mar 19 18:59:16 -05 2025: Cargando módulo python/3.6.0_miniconda-4.3.11_gcc-11.2.0"
module load python/3.6.0_miniconda-4.3.11_gcc-11.2.0
echo "Wed Mar 19 18:59:16 -05 2025: activando venv abc4sfa"
source activate abc4sfa

echo "Wed Mar 19 18:59:16 -05 2025: Ejecutando Rscript negacion_empirica/jobs/main_id=s1.R"
Rscript negacion_empirica/jobs/main_id=s1.R
