#!/bin/bash

# Verifica que se proporcione un argumento
if [ "$#" -ne 1 ]; then
    echo "Uso: $0 <ruta_al_proyecto>"
    exit 1
fi

PROYECTO_DIR="$1"
CONFIG_FILE="$PROYECTO_DIR/params.config"
TEMPLATE_R="$PROYECTO_DIR/main_template.R"

# Verifica que existan los archivos necesarios
if [ ! -f "$CONFIG_FILE" ] || [ ! -f "$TEMPLATE_R" ]; then
    echo "ERROR: No se encontraron params.config o main_template.R en $PROYECTO_DIR"
    exit 1
fi

# Carga los parámetros desde el archivo de configuración
source "$CONFIG_FILE"

# Crea las carpetas si no existen
mkdir -p "$PROYECTO_DIR/$aa_JOBS_folder"
mkdir -p "$PROYECTO_DIR/$aa_LOGS_folder"

# Convierte los arrays declarados en el .config a variables locales
PARAM1_LIST=("${aa_param1[@]}")
PARAM2_LIST=("${aa_param2[@]}")

# Verifica si hay valores en los arrays
if [ ${#PARAM1_LIST[@]} -eq 0 ] || [ ${#PARAM2_LIST[@]} -eq 0 ]; then
  echo "Error: Los parámetros están vacíos"
  exit 1
fi

# Genera los jobs dependiendo del modo
if [ "$aa_USE_ALL_COMBINATIONS" = true ]; then
    # Genera todas las combinaciones posibles de PARAM1 y PARAM2
    for p1 in "${PARAM1_LIST[@]}"; do
        for p2 in "${PARAM2_LIST[@]}"; do
            JOB_FILE="$PROYECTO_DIR/$aa_JOBS_folder/job_${p1}_${p2}.sh"
            R_FILE="$PROYECTO_DIR/$aa_JOBS_folder/main_${p1}_${p2}.R"

            # Genera el script SLURM
            cat > "$JOB_FILE" <<EOF
#!/bin/bash
#SBATCH --partition=$aa_PARTITION
#SBATCH --nodes=$aa_NODES
#SBATCH --ntasks=$aa_TASKS
#SBATCH --time=$aa_TIME
#SBATCH --job-name=job_${p1}_${p2}
#SBATCH --output=$PROYECTO_DIR/$aa_LOGS_folder/job_${p1}_${p2}.out
#SBATCH --error=$PROYECTO_DIR/$aa_LOGS_folder/job_${p1}_${p2}.err

module load $aa_PYTHON_MODULE
source activate $aa_CONDA_ENV

Rscript $R_FILE
EOF
            chmod +x "$JOB_FILE"

            # Genera el script R reemplazando los valores
            sed "s/__param1__/$p1/g; s/__param2__/$p2/g" "$TEMPLATE_R" > "$R_FILE"

            echo "Generado: $JOB_FILE y $R_FILE"
        done
    done
else
    # Genera solo pares ordenados (índices coincidentes)
    for ((i=0; i<${#PARAM1_LIST[@]}; i++)); do
        p1="${PARAM1_LIST[i]}"
        p2="${PARAM2_LIST[i]}"

        JOB_FILE="$PROYECTO_DIR/$aa_JOBS_folder/job_${p1}_${p2}.sh"
        R_FILE="$PROYECTO_DIR/$aa_JOBS_folder/main_${p1}_${p2}.R"

        # Genera el script SLURM
        cat > "$JOB_FILE" <<EOF
#!/bin/bash
#SBATCH --partition=$aa_PARTITION
#SBATCH --nodes=$aa_NODES
#SBATCH --ntasks=$aa_TASKS
#SBATCH --time=$aa_TIME
#SBATCH --job-name=job_${p1}_${p2}
#SBATCH --output=$PROYECTO_DIR/$aa_LOGS_folder/job_${p1}_${p2}.out
#SBATCH --error=$PROYECTO_DIR/$aa_LOGS_folder/job_${p1}_${p2}.err

module load $aa_PYTHON_MODULE
source activate $aa_CONDA_ENV

Rscript $R_FILE
EOF
        chmod +x "$JOB_FILE"

        # Genera el script R reemplazando los valores
        sed "s/__param1__/$p1/g; s/__param2__/$p2/g" "$TEMPLATE_R" > "$R_FILE"

        echo "Generado: $JOB_FILE y $R_FILE"
    done
fi