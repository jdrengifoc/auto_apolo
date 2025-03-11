#!/bin/bash

# Verifica que se proporcione un argumento
if [ "$#" -ne 1 ]; then
    echo "Uso: $0 <ruta_al_proyecto>"
    exit 1
fi

PROYECTO_DIR="experiment_1"
CONFIG_FILE="$PROYECTO_DIR/params.config"
TEMPLATE_R="$PROYECTO_DIR/main_template.R"
JOBS_DIR="$PROYECTO_DIR/aa_JOBS_DIR"
LOGS_DIR="$PROYECTO_DIR/apolo_info"

# Verifica que existan los archivos necesarios
if [ ! -f "$CONFIG_FILE" ] || [ ! -f "$TEMPLATE_R" ]; then
    echo "ERROR: No se encontraron params.config o main_template.R en $PROYECTO_DIR"
    exit 1
fi

# Crea las carpetas si no existen
mkdir -p "$JOBS_DIR"
mkdir -p "$LOGS_DIR"

# Carga los parámetros desde el archivo de configuración
source "$CONFIG_FILE"

# Itera sobre las combinaciones de parámetros
for ((i = 0; i < ${#param1[@]}; i++)); do
    VAL1="${param1[$i]}"
    VAL2="${param2[$i]}"

    # Nombres de archivos
    JOB_FILE="$JOBS_DIR/job_${VAL1}_${VAL2}.sh"
    R_FILE="$JOBS_DIR/main_${VAL1}_${VAL2}.R"

    # Genera el script SLURM
    cat > "$JOB_FILE" <<EOF
#!/bin/bash
#SBATCH --partition=longjobs
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --time=14-00:00:00
#SBATCH --job-name=job_${VAL1}_${VAL2}
#SBATCH --output=$LOGS_DIR/job_${VAL1}_${VAL2}.out
#SBATCH --error=$LOGS_DIR/job_${VAL1}_${VAL2}.err

module load python/3.6.0_miniconda-4.3.11_gcc-11.2.0
source activate abc4sfa

Rscript $R_FILE
EOF

    chmod +x "$JOB_FILE"

    # Genera el script R reemplazando los valores
    sed "s/__param1__/$VAL1/g; s/__param2__/$VAL2/g" "$TEMPLATE_R" > "$R_FILE"

    echo "Generado: $JOB_FILE y $R_FILE"
done