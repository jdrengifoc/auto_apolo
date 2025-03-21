#!/bin/bash

# Verify that the project directory is passed as an argument
if [ "$#" -ne 1 ]; then
    echo "Uso: $0 <ruta_al_proyecto>"
    exit 1
fi
PROYECTO_DIR="$1"
# Verify that the folder exists
if [ ! -d "$PROYECTO_DIR" ]; then
    echo "ERROR: La carpeta '$PROYECTO_DIR' no existe."
    exit 1
fi


# Find the first `.config` and `.R` files in the project folder.
CONFIG_FILE=$(find "$PROYECTO_DIR" -maxdepth 1 -type f -name "*.config" | head -n 1)
TEMPLATE_R=$(find "$PROYECTO_DIR" -maxdepth 1 -type f -name "*.R" | head -n 1)

# Check if files were found.
if [ -z "$CONFIG_FILE" ] || [ -z "$TEMPLATE_R" ]; then
    echo "ERROR: No se encontró un archivo .config o .R en $PROYECTO_DIR"
    exit 1
fi
echo "✅ CONFIG_FILE encontrado: $CONFIG_FILE"
echo "✅ TEMPLATE_R encontrado: $TEMPLATE_R"

# Load the parameters from the configuration file.
source "$CONFIG_FILE"

# Crea las carpetas si no existen
mkdir -p "$PROYECTO_DIR/$aa_JOBS_folder"
mkdir -p "$PROYECTO_DIR/$aa_LOGS_folder"

# Genera las combinaciones de parámetros usando el script externo
COMBINACIONES=($(./utils/gen_params_combinations.sh "$PROYECTO_DIR"))

# Verifica si se generaron combinaciones
if [ ${#COMBINACIONES[@]} -eq 0 ]; then
    echo "ERROR: No se generaron combinaciones de parámetros."
    exit 1
fi

param_names=($(echo "${!aa_PARAMS[@]}" | tr ' ' '\n' | sort))

# Genera los scripts de trabajo
for combo in "${COMBINACIONES[@]}"; do
    IFS=',' read -r -a values <<< "$combo"

    # Build file name using parameter names
    filename=""
    for ((j=0; j<${#param_names[@]}; j++)); do
        filename+="${param_names[j]}=${values[j]}_"
    done
    filename="${filename%_}"  # Remove trailing underscore

    JOB_FILE="$PROYECTO_DIR/$aa_JOBS_folder/job_${filename}.sh"
    R_FILE="$PROYECTO_DIR/$aa_JOBS_folder/main_${filename}.R"

    # Genera el script SLURM
    cat > "$JOB_FILE" <<EOF
#!/bin/bash
#SBATCH --partition=$aa_PARTITION
#SBATCH --nodes=$aa_NODES
#SBATCH --ntasks=$aa_TASKS
#SBATCH --time=$aa_TIME
#SBATCH --job-name=job_${filename}
#SBATCH --output=$PROYECTO_DIR/$aa_LOGS_folder/job_${filename}.out
#SBATCH --error=$PROYECTO_DIR/$aa_LOGS_folder/job_${filename}.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=$aa_EMAIL

echo "$(date): Cargando módulo $aa_PYTHON_MODULE"
module load $aa_PYTHON_MODULE
echo "$(date): activando venv $aa_CONDA_ENV"
source activate $aa_CONDA_ENV

echo "$(date): Ejecutando Rscript $R_FILE"
Rscript $R_FILE
EOF
    chmod +x "$JOB_FILE"

    # Genera el script R reemplazando los valores en la plantilla
    R_CONTENT=$(<"$TEMPLATE_R")
    for ((i=0; i<${#values[@]}; i++)); do
        R_CONTENT=$(echo "$R_CONTENT" | sed "s/__${param_names[i]}__/${values[i]}/g")
    done
    echo "$R_CONTENT" > "$R_FILE"
    
    echo "✅ Generado: $JOB_FILE y $R_FILE"
done